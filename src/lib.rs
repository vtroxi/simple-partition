mod annealing;
mod graph;
mod partition_util;

pub use annealing::*;
pub use graph::*;
pub use partition_util::*;

pub(crate) fn divide_round_up(a: u32, b: u32) -> u32 {
    (a + b - 1) / b
}
