#include <vector>
#include <unordered_map>

#include "core/graph.hpp"
#include "fma-common/string_formatter.h"
#include "fma-common/check_date.h"
#include "fma-common/logging.h"
#include "toolkits/config_common.h"
#include "fma-common/license.h"

#define MAX_VID 1099511627775ul
double* new_m;    //记录节点的二度邻居的总链接数


size_t parse_edge(const char *p, const char * end, EdgeUnit<double> & e) {
  const char * orig = p;
  int64_t t = 0;
  size_t r = 0;
  if (*p == '#') {
    while (*p != '\n') p++;
    e.src = MAX_VID;
    e.dst = MAX_VID;
    return p - orig;
  }
  r = fma_common::TextParserUtils::ParseInt64(p, end, t);
  e.src = t;
  p += r;
  while (p != end && (*p == ' ' || *p == '\t')) p++;
  r = fma_common::TextParserUtils::ParseInt64(p, end, t);
  e.dst = t;
  p += r;
  double w = 1;
  //r = fma_common::TextParserUtils::ParseDigit(p, end, w);
  e.edge_data = w;
  if (e.src == e.dst) e.edge_data /= 2;
  //p += r;
  return p - orig;
}

bool filter_edge(EdgeUnit<double> & e) {
  return e.src != MAX_VID && e.dst != MAX_VID;
}

class Louvain_graph {
  Graph<double>* my_graph;
  double m;
  VertexId* label;
  Bitmap* active;
  double* e_tot;   //e_tot means degree of each community
  double* k;      //k means edges weight sum of each vertex
  VertexId comm_num;
  double Q;
  double LQ;

public:
  Louvain_graph(Graph<double>* _graph){
   my_graph = _graph;
  }

  VertexId update_comm_num(){
    comm_num = my_graph->stream_vertices<VertexId>(
      [&](VertexId dst){
        if(e_tot[dst] == 0){
          return 0;
        }
        return 1;
      },
      active
    );
    printf("community num is %lu\n",comm_num);
    return comm_num;
  }

  VertexId get_comm_num(){
    return comm_num;
  }

  double update_Q() {
    Q = my_graph->stream_vertices<double> (
      [&] (VertexId v) {
        double q = 0.0;
        for (auto e : my_graph->out_edges(v)) {
          VertexId nbr = e.neighbour;
          if (label[v] == label[nbr]) q += e.edge_data;
        }
        q -= 1.0 * k[v] * e_tot[label[v]] / (2 * m);
        return q;
      },
      active
    ) / (2.0 * m);
    printf("Q =  %f\n",Q);
    return Q;
  }

  double update_LQ() {
    LQ = my_graph->stream_vertices<double> (
      [&] (VertexId v) {
        double temp_m = new_m[v];
        // double temp_m = m;
        // printf("vertexId = %lu,temp_m = %f\n",v,new_m[v]);
        double q = 0.0;
        for (auto e : my_graph->out_edges(v)) {
          VertexId nbr = e.neighbour;
          if (label[v] == label[nbr]) q += e.edge_data;
        }
        q -= 1.0 * k[v] * e_tot[label[v]] / (2 * temp_m);
        // printf("vertexId = %lu,LQ = %f\n",v,q);
        return q / (2.0 * temp_m);
      },
      active
    );
    printf("LQ = %f\n",LQ);
    return LQ;
  }

  double get_Q() {
    return Q;
  }

  void update_e_tot() {
    my_graph->fill_vertex_array(e_tot, (double)0.0);
    my_graph->stream_vertices<VertexId>(
      [&] (VertexId v) {
        write_add(&e_tot[label[v]], k[v]);
        return 0;
      },
      active
    );
  }

  void update_all() {
    update_e_tot();
    update_comm_num();
    update_Q();
    update_LQ();
    update_cp();
  }

  void update_cp(){
    double A = 0.0;
    double B = 0.0;
    VertexId* count = my_graph->alloc_vertex_array<VertexId>();
    my_graph -> fill_vertex_array(count, (VertexId)0);
    my_graph->stream_vertices<char>(
      [&](VertexId v){
        write_add(&count[label[v]],(VertexId)1);
        return 0;
      },
      active
    );
    B = my_graph->stream_vertices<VertexId>(
      [&](VertexId v){
        if(count[v] != 0){
          return count[v]*(count[v]-1);
        }
        return (VertexId)0;
      },
      active
    );
    A = my_graph->stream_vertices<VertexId>(
      [&](VertexId v){
        for(auto e : my_graph->out_edges(v)){
          VertexId s = e.neighbour;
          if(label[v] == label[s]){
            return 1;
          }
        }
        return 0;
      },
      active
    );
    printf("graph condition probability = %f\n",A / B);
  }

  void init(){
    if(my_graph == NULL){
      exit(0);
    }
    active = my_graph->alloc_vertex_bitmap();
    label = my_graph->alloc_vertex_array<VertexId>();
    k = my_graph->alloc_vertex_array<double>();
    e_tot = my_graph->alloc_vertex_array<double>();

    active->fill();
    my_graph->fill_vertex_array(k, (double)0.0);
    my_graph->fill_vertex_array(e_tot, (double)0.0);

    m = my_graph->stream_vertices<double>(
      [&](VertexId dst){
        label[dst] = dst;
        for(auto e : my_graph->out_edges(dst)){
          k[dst] += e.edge_data;
        }
        e_tot[dst] = k[dst];
        return k[dst];
      },
      active
    ) / 2;
    printf("m = %f\n",m);
  }

  void update_new_m(){
    double init_time = 0.0;
    init_time -= get_time();
    my_graph->stream_vertices<char>(
      [&](VertexId v){
        double temp_m = 0.0;
        std::unordered_map<VertexId,int> temp;
        for (auto e : my_graph->out_edges(v)) {
          VertexId nbr = e.neighbour;
          temp[nbr] = 0;
          for(auto m : my_graph->out_edges(nbr)){
            VertexId nbr2 = m.neighbour;
            temp[nbr2] = 0;
          }
        }
        for(auto v : temp){
          temp_m += my_graph->out_degree[v.first];
        }
        new_m[v] = temp_m;
        return 0;
      },
      active
    );
    printf("exec_time = %lf(s),finding second neighbour over\n",init_time + get_time());
  }
  //init over
  VertexId async_louvain(){
    VertexId active_vertices = my_graph->stream_vertices<VertexId>(
      [&](VertexId dst){
        std::unordered_map<VertexId, double> count;
        for (auto e : my_graph->out_edges(dst)) {
          if (dst == e.neighbour) continue;
          VertexId src = e.neighbour;
          auto it = count.find(label[src]);
          if (it == count.end()) {
            count[label[src]] = e.edge_data;
          }else{
            it->second += e.edge_data;
          }
        }
        VertexId best_c = label[dst];
        // double delta_in = 0.0;
        //diffrent init of delta_in
        VertexId old_label = label[dst];

        double k_in_out = 0.0;
        //k_in_out 表示待移动节点和当前社区的连接权重和
        if (count.find(old_label) != count.end()) {
          k_in_out = count[old_label];
        }
        double delta_in = k[dst] * (e_tot[old_label]-k[dst]) - 2.0 * k_in_out * m;
        delta_in = -delta_in;
        //delta_in init over

        for (auto & ele : count){
          if(label[dst] != ele.first){
            double k_i_in = ele.second;
            double delta_out = k_i_in*2.0*m - k[dst]*e_tot[ele.first];
            if(delta_out > delta_in){
              best_c = ele.first;
              delta_in = delta_out;
            }else
            if(delta_out == delta_in && ele.first < best_c){
              best_c = ele.first;
            }
          }
        }
        if(best_c != label[dst]){
          write_add(&e_tot[label[dst]],-k[dst]);
          my_graph->lock_vertex(dst);
          label[dst] = best_c;
          write_add(&e_tot[best_c],k[dst]);
          my_graph->unlock_vertex(dst);
          return 1;
        }
        return 0;
      },
      active
    );
    return active_vertices;
  }
  //one thread
  VertexId single_louvain(){
    VertexId active_vertices = 0;
    for(VertexId dst = 0;dst < my_graph->get_num_vertices();dst++){
        std::unordered_map<VertexId, double> count;
        for (auto e : my_graph->out_edges(dst)) {
          if (dst == e.neighbour) continue;
          VertexId src = e.neighbour;
          auto it = count.find(label[src]);
          if (it == count.end()) {
            count[label[src]] = e.edge_data;
          }else{
            it->second += e.edge_data;
          }
        }
        VertexId best_c = label[dst];
        VertexId old_label = label[dst];

        double k_in_out = 0.0;
        if (count.find(old_label) != count.end()) {
          k_in_out = count[old_label];
        }
        double delta_in = k[dst] * (e_tot[old_label]-k[dst]) - 2.0 * k_in_out * m;
        delta_in = -delta_in;

        for (auto & ele : count){
          if(label[dst] != ele.first){
            double k_i_in = ele.second;
            double delta_out = k_i_in*2.0*m - k[dst]*e_tot[ele.first];
            if(delta_out > delta_in){
              best_c = ele.first;
              delta_in = delta_out;
            }else
            if(delta_out == delta_in && ele.first < best_c){
              best_c = ele.first;
            }
          }
        }
        if(best_c != label[dst]){
          write_add(&e_tot[label[dst]],-k[dst]);
          my_graph->lock_vertex(dst);
          label[dst] = best_c;
          write_add(&e_tot[best_c],k[dst]);
          my_graph->unlock_vertex(dst);
          active_vertices++;
        }
      }
    return active_vertices;
  }

  //first iter of Louvain
  void Louvain_propagate(){
    printf("enter Louvain_propagate\n");
    VertexId active_vertices = my_graph->get_num_vertices();
    // LOG()<<"active_vertices = "<<active_vertices;
    VertexId old_vertices = active_vertices;
    int cyc = 0;
    int comm_cyc = 0;
    while(active_vertices > 0){
      old_vertices = active_vertices;
      active_vertices = single_louvain();
      printf("active_vertices(%d) = %lu\n", cyc, active_vertices);
        cyc++;
        if(old_vertices == active_vertices){
          comm_cyc++;
        }else{
          comm_cyc = 0;
        }

        if(comm_cyc == 5){
          printf("Something Wrong\n");
          break;
        }
    }
    update_comm_num();
    update_Q();
  }

  void update_by_subgraph(){
    //label,my_graph,
    long int * sub_index = my_graph->alloc_vertex_array<long int>();
    VertexId num_sub_vertices = 0;      //vertices of new graph
    for (VertexId i = 0; i < my_graph->get_num_vertices();i ++) {
      if (e_tot[i] == 0) {
        sub_index[i] = -1;
      } else {
        sub_index[i] = num_sub_vertices;
        num_sub_vertices++;
      }
    }

    std::vector<std::unordered_map<VertexId, double> > sub_edges(num_sub_vertices);

    for(VertexId i = 0;i < my_graph->get_num_vertices();i ++){
      if(sub_index[label[i]] == -1) continue;
      for(auto e : my_graph->out_edges(i)){
        //only save one direction out_edges
        if(i <= e.neighbour){
          double e_double = e.edge_data;
          if(i == e.neighbour){
            e_double /= 2;
          }
          auto iter = sub_edges[sub_index[label[i]]].find(sub_index[label[e.neighbour]]);
          if(iter == sub_edges[sub_index[label[i]]].end()){
            sub_edges[sub_index[label[i]]][sub_index[label[e.neighbour]]] = e_double;
          }else{
            iter->second += e_double;
          }
        }
      }
    }

    VertexId sub_edges_num = 0;
    for(VertexId i = 0;i < sub_edges.size();i++){
      sub_edges_num += sub_edges[i].size();
    }

    EdgeUnit<double> * sub_edge_array = new EdgeUnit<double>[sub_edges_num];
    VertexId index = 0;
    for(VertexId sub_node = 0;sub_node < num_sub_vertices;sub_node ++){
      for(auto ele : sub_edges[sub_node]){
        assert(index < sub_edges_num);
        sub_edge_array[index].src = sub_node;
        sub_edge_array[index].dst = ele.first;
        sub_edge_array[index].edge_data = ele.second;
        index++;
      }
    }
    sub_edges.clear();

    Graph<double> * sub_graph;
    sub_graph = new Graph<double>();
    sub_graph->load_from_array(*sub_edge_array, num_sub_vertices, sub_edges_num, true);
    delete sub_edge_array;
    printf("subgraph build over\n");

    VertexId* sub_to_parent_map = sub_graph -> alloc_vertex_array<VertexId>();

    sub_graph->fill_vertex_array(sub_to_parent_map, my_graph->get_num_vertices());
    my_graph->stream_vertices<VertexId> (
      [&] (VertexId v) {
        if (sub_index[v] < 0) {
          return 0;
        }
        sub_to_parent_map[sub_index[v]] = v;
        return 1;
      },
      active
    );

    Louvain_graph* sub_my_graph = new Louvain_graph(sub_graph);
    sub_my_graph->init();
    sub_my_graph->Louvain_propagate();

    if (sub_my_graph->get_comm_num() == sub_graph->get_num_vertices()) {
      return;
    } else {
      sub_my_graph->update_by_subgraph();
    }

    //将当前sub_graph的点映射回去
    my_graph->stream_vertices<VertexId> (
      [&] (VertexId v) {
        if (sub_index[label[v]] <= 0) {
          return 0;
        }
        VertexId sub_index_v_comm = sub_my_graph->label[sub_index[label[v]]];
        label[v] = sub_to_parent_map[sub_index_v_comm];
        return 0;
      },
      active
    );
  }
  void output(std::string output_dir) {
    int num_threads = omp_get_num_procs();
    fma_common::OutputFmaStream * fout = new fma_common::OutputFmaStream [num_threads];
    for (int t_i = 0; t_i < num_threads; t_i++) {
      std::string filename = fma_common::StringFormatter::Format("{}/part-r-{}", output_dir, t_i);
      fout[t_i].Open(filename, 64 << 20);
    }

    #pragma omp parallel for
    for (VertexId v_i = 0; v_i < my_graph->get_num_vertices(); v_i++) {
        if (my_graph->out_degree[v_i] > 0) {
            int t_i = omp_get_thread_num();
            std::string line = fma_common::StringFormatter::Format("{} {}\n", v_i, label[v_i]);
            fout[t_i].Write(line.c_str(), line.size());
        }
    }
    for (int t_i = 0; t_i < num_threads; t_i++) {
      fout[t_i].Close();
    }
    LOG() << "writing to files done!";
  }
};

int main(int argc, char ** argv) {
  int input_format = 1;
  std::string input_dir = "";
  VertexId vertices = 0;
  std::string output_dir = "";
  std::string license_file = "";

  fma_common::Configuration config;
  toolkits_common::config_common(config, license_file, input_format, input_dir, vertices, true, output_dir);

  config.Parse(argc, argv);
  config.ExitAfterHelp();
  config.Finalize();

  double exec_time = 0;
  exec_time -= get_time();

  Graph<double> * graph;
  graph = new Graph<double>();
  if (input_format == 0) {
    graph->load(input_dir, vertices, true);
  } else {
    graph->load_txt_undirected(input_dir, vertices, parse_edge, filter_edge);
  }
  exec_time += get_time();
  printf("prepare time = %lf(s)\n", exec_time);
  exec_time = 0;
  exec_time -= get_time();
  Louvain_graph * l_graph = new Louvain_graph(graph);
  l_graph->init();
  new_m = graph->alloc_vertex_array<double>();
  l_graph->update_new_m();
  //init测试完成
  l_graph->Louvain_propagate();
  l_graph->update_by_subgraph();
  l_graph->update_all();
  exec_time += get_time();
  printf("exec_time = %lf(s)\n", exec_time);
  l_graph->output(output_dir);
}
