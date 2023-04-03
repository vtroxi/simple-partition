use crate::{divide_round_up, Graph};
use rand::rngs::StdRng;
use rand::Rng;
use std::collections::VecDeque;

#[derive(Copy, Clone, Eq, PartialEq, Debug)]
pub enum InitialPartitioningMethod {
    /// Vertices are initialized to partitions based on their index in the graph (idx % partitions).
    Modulo,
    /// Vertices are initialized to random partitions.
    Random,
    /// Vertex patches get initialized based on locality using breadth first search.
    Bfs,
}

impl Graph {
    pub fn partition_initial(
        &mut self,
        method: InitialPartitioningMethod,
        rng: &mut StdRng,
        partitions: u32,
    ) {
        match method {
            InitialPartitioningMethod::Modulo => {
                for i in 0..self.vertices.len() {
                    self.vertices[i].color = i as u32 % partitions;
                }
            }
            InitialPartitioningMethod::Random => {
                for i in 0..self.vertices.len() {
                    self.vertices[i].color = rng.gen_range(0..partitions);
                }
            }
            InitialPartitioningMethod::Bfs => {
                for v in self.vertices.iter_mut() {
                    v.color = u32::MAX;
                }

                let target_size = divide_round_up(self.vertices.len() as u32, partitions);
                //let target_size = self.vertices.len() as u32 / partitions;

                let mut color = 0u32;
                let mut partition_size = 0u32;
                let mut queue = VecDeque::new();

                let mut visited = 0;
                while visited < self.vertices.len() {
                    queue.clear();
                    let start = self
                        .vertices
                        .iter()
                        .position(|v| v.color == u32::MAX)
                        .unwrap();
                    queue.push_back(start as u32);

                    while let Some(vx) = queue.pop_front() {
                        if self.vertices[vx as usize].color == u32::MAX {
                            self.vertices[vx as usize].color = color;

                            visited += 1;
                            partition_size += 1;
                            if partition_size >= target_size {
                                color += 1;
                                partition_size = 0;
                                break;
                            }

                            for e in self.vertices[vx as usize].edges.iter() {
                                if self.vertices[e.dst as usize].color == u32::MAX {
                                    queue.push_back(e.dst);
                                }
                            }
                        }
                    }
                }
                log::trace!("bsf partitions {}", color + 1);
            }
        }
    }

    pub(crate) fn swap_colors(&mut self, va: u32, vb: u32) {
        let tmp = self.vertices[va as usize].color;
        self.vertices[va as usize].color = self.vertices[vb as usize].color;
        self.vertices[vb as usize].color = tmp;
    }
}
