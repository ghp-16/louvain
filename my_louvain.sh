for i in {1..5}
do
    echo "<====     my_louvain  livejournal     ====>"
    ./bin/my_louvain --input_dir /data/disk4/clk/testdata/livejournal-4.85M/txt-no-edge/ --output_dir ~/GeminiLite/output/

done

for i in {1..5}
do
    echo "<====     my_louvain  facebook     ====>"
    ./bin/my_louvain --input_dir /data/disk4/clk/testdata/facebook-4039/txt_no_edge/ --output_dir ~/GeminiLite/output/
done
for i in {1..5}
do
    echo "<====     my_louvain  twitter-half     ====>"
    ./bin/my_louvain --input_dir /data/disk4/clk/testdata/twitter/twitter-txt/ --output_dir ~/GeminiLite/output/
done

