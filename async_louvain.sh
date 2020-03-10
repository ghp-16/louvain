for i in {1..5}
do
    echo "<====     async_louvain  livejournal     ====>"
    ./bin/async_louvain --input_dir /data/disk1/clk/livejournal-txt/ --output_dir /home/ghp/GeminiLite/out/
done

for i in {1..5}
do
    echo "<====     async_louvain  facebook     ====>"
    ./bin/async_louvain --input_dir /data/disk1/clk/facebook-txt/ --output_dir /home/ghp/GeminiLite/out/
done
./bin/async_louvain --input_dir /home/ghp/twitter/ --output_dir /home/ghp/GeminiLite/out/
for i in {1..5}
do
    echo "<====     async_louvain  twitter-half     ====>"
    ./bin/async_louvain --input_dir /home/ghp/twitter/ --output_dir /home/ghp/GeminiLite/out/
done

