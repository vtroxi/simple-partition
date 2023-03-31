use simple_partition::{Graph, InitialPartitioning, PartitioningConfig};
use std::fs::File;
use std::time;

fn main() {
    env_logger::init();
    let mut graph =
        Graph::deserialize_metis(&mut File::open("graphs/4elt.graph").unwrap()).unwrap();
    let t1 = time::Instant::now();
    graph.partition(
        &PartitioningConfig {
            initial_partitioning: InitialPartitioning::Bfs,
            ..Default::default()
        },
        4,
    );
    println!("time: {}ms", t1.elapsed().as_millis());
    let size_info = graph.partition_sizes();
    println!(
        "partitions: {}, with sizes from {} to {} with an average of {})",
        graph.count_partitions(),
        size_info.0,
        size_info.1,
        size_info.2
    );
}
