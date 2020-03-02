#include <vector>
#include <unordered_map>

#include "core/graph.hpp"
#include "fma-common/string_formatter.h"
#include "fma-common/check_date.h"
#include "fma-common/logging.h"
#include "toolkits/config_common.h"
#include "fma-common/license.h"

#define MAX_VID 1099511627775ul


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
  double* k;      //k means degree of each vertex

public:
  Louvain_graph(Graph<double>* _graph){
   my_graph = *graph;
  }

  void init(){
   if(my_graph == NULL){
     exit(0);
   }
   active = my_graph->alloc_vertex_bitmap();
   label = my_graph->alloc_vertex_array<VertexId>();
   active->fill();
    graph->fill_vertex_array(k, (double)0.0);
    graph->fill_vertex_array(e_tot, (double)0.0);

   m = my_graph->stream_vertex<VertexId>(
      [&](VertexId dst){
        label[dst] = dst;
        for(auto e : my_graph->out_edges(dst)){
          VertexId src = dst.neighbour;
          k[dst] += src.edge_data;
        }
        e_tot[dst] = k[dst];
        return k[dst]
      },
      active
    ) / 2;
    LOG()<<"m = "<<m;
  }
  //init over

  VertexId louvain_first(){
    VertexId active_vertices = stream_vertex(

      std::unordered_map<VertexId,double> count;

      [&](VertexId dst){
        for(auto e : my_graph->out_edges(dst)){
          VertexId src = e.neighbour;
          if(label[src] != label[dst]){
            auto it = count.find(label[src]);
            if(it != count.end()){
              it->second += src.edge_data;
            }else{
              it->second = src.edge_data;
            }
          }
        }
        VertexId best_c = label[src];
        double delta_in = 0.0;

        for (auto & ele : count){
          if(ele != ele.first){
            double delta_out = count[ele]*2.0*m - k[dst]*e_tot[ele.second];
            if(delta_out > delta_in){
              best_c = ele.second;
              delta_in = delta_out;
            }else
            if(delta_out == delta_in && ele.second < best_c){
              best_c = ele.second;
            }
          }
        }
        if(best_c != label[src]){
          write_add(&e_tot[label[src]],-k[src]);
          my_graph->lock_vertex(src);
          label[src] = best_c;
          my_graph->unlock_vertex(src);
          return 1;
        }
        return 0;
      },
      active
    );
    return active_vertices;
  }

}

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

  Graph<double> * graph;
  graph = new Graph<double>();
  if (input_format == 0) {
    graph->load(input_dir, vertices, true);
  } else {
    graph->load_txt_undirected(input_dir, vertices, parse_edge, filter_edge);
  }
  Louvain_graph * l_graph = new LouvainGraph(graph, active_threshold);
  l_graph->init();
}