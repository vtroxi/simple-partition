use std::time;
use simple_partition::{Graph, PartitioningConfig};

fn main() {
	let t1 = time::Instant::now();
	let mut graph = Graph::deserialize_metis("graphs/4elt.graph").unwrap();
	//dbg!(&graph.vertices);
	graph.partition(&PartitioningConfig::default(), 4);
	println!("time: {}ms", t1.elapsed().as_millis());
}