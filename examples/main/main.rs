use simple_partition::{
    AnnealingPartitioningConfig, BisectionPartitioningConfig, Graph, InitialPartitioningMethod,
};
use std::fs::File;
use std::time;

fn main() {
    env_logger::init();
    let mut graph =
        Graph::deserialize_metis(&mut File::open("graphs/4elt.graph").unwrap()).unwrap();
    let t1 = time::Instant::now();
    // graph.partition_annealing(
    //     &AnnealingPartitioningConfig {
    //         initial_partitioning: InitialPartitioningMethod::Bfs,
    //         ..Default::default()
    //     },
    //     4,
    // );
    graph.partition_bisection(&BisectionPartitioningConfig {
        initial_partitioning: InitialPartitioningMethod::Bfs,
        ..Default::default()
    });
    println!("time: {}ms", t1.elapsed().as_millis());
    let infos = graph.partition_info();
    println!(
        "partitions: {}, with sizes from {} to {} with an average of {})",
        graph.count_partitions(),
        infos.0,
        infos.1,
        infos.2
    );
    let sizes = graph.partition_sizes();
    for c in 0..=graph.max_partition_color() {
        let s = sizes[c as usize];
        println!(
            "partition {c} with size {s} has {} parts",
            graph.count_partition_parts(c, s)
        );
    }
}
