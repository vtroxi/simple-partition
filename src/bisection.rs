use crate::{divide_round_up, Graph, InitialPartitioningMethod};

pub struct BisectionPartitioningConfig {
    /// The seed for the random number generator.
    pub rng_seed: u64,
    /// What initial partitioning method to use.
    pub initial_partitioning: InitialPartitioningMethod,
}

impl Default for BisectionPartitioningConfig {
    fn default() -> Self {
        Self {
            rng_seed: 1234,
            initial_partitioning: InitialPartitioningMethod::Modulo,
        }
    }
}

impl Graph {
    /// Splits the graph into two parts while minimizing the edge cut cost.
    pub fn partition_bisection(&mut self, config: &BisectionPartitioningConfig) {
        todo!()
    }
}
