for i in {1..5}
do
    echo "<====     my_louvain  livejournal     ====>"
    ./bin/my_louvain --input_dir /data/disk1/clk/livejournal-txt/ --output_dir /home/ghp/GeminiLite/out/
done

for i in {1..5}
do
    echo "<====     my_louvain  facebook     ====>"
    ./bin/my_louvain --input_dir /data/disk1/clk/facebook-txt/ --output_dir /home/ghp/GeminiLite/out/
done
./bin/async_louvain --input_dir /home/ghp/twitter/ --output_dir /home/ghp/GeminiLite/out/
for i in {1..5}
do
    echo "<====     my_louvain  twitter-half     ====>"
    ./bin/my_louvain --input_dir /home/ghp/twitter/ --output_dir /home/ghp/GeminiLite/out/
done

