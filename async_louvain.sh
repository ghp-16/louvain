for i in {1..5}
do
    echo "<====     async_louvain  livejournal     ====>"
    ./bin/async_louvain --input_dir /data/disk4/clk/testdata/livejournal-4.85M/txt-no-edge/ --output_dir ~/GeminiLite/output/
done

for i in {1..5}
do
    echo "<====     async_louvain  facebook     ====>"
    ./bin/my_louvain --input_dir /data/disk4/clk/testdata/facebook-4039/txt_no_edge/ --output_dir ~/GeminiLite/output/
done
./bin/async_louvain --input_dir /home/ghp/twitter/ --output_dir /home/ghp/GeminiLite/out/
for i in {1..5}
do
    echo "<====     async_louvain  twitter-half     ====>"
    /bin/my_louvain --input_dir /data/disk4/clk/testdata/twitter/twitter-txt/ --output_dir ~/GeminiLite/output/
done

