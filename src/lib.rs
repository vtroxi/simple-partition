// Graph partitioning using simulated annealing based on Ja-be-Ja:
// https://www.diva-portal.org/smash/get/diva2:1043244/FULLTEXT01.pdf

use anyhow::{Context, Result};
use rand::rngs::StdRng;
use rand::{Rng, SeedableRng};
use std::collections::VecDeque;
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::Path;
use std::time;

#[derive(Copy, Clone, Debug)]
pub struct GraphEdge {
    pub dst: u32,
    pub weight: u32,
}

#[derive(Clone, Debug)]
pub struct GraphVertex {
    pub edges: Vec<GraphEdge>,
    pub color: u32,
}

#[derive(Clone)]
pub struct Graph {
    pub vertices: Vec<GraphVertex>,
}

impl Graph {
    pub const MAX_EDGES_PER_VERTEX: u32 = 256;

    pub fn deserialize_metis<P: AsRef<Path>>(path: P) -> Result<Self> {
        let reader = BufReader::new(File::open(path)?);

        let mut lines = reader.lines().map(|l| l.unwrap());

        // Parse the header line
        let header = lines.next().context("could not get header line")?;
        let header_parts = header.split_ascii_whitespace().collect::<Vec<_>>();
        let vertex_count = header_parts[0]
            .parse::<usize>()
            .context("could not parse vertex count")?;
        //let edge_count = header_parts[1].parse::<usize>().context("could not parse edge count")?;

        // Initialize the graph
        let mut graph = Graph {
            vertices: vec![
                GraphVertex {
                    edges: vec![],
                    color: u32::MAX
                };
                vertex_count
            ],
        };

        // Parse the edges
        let mut src = 0;
        for line in lines {
            if line.starts_with("%") || line.starts_with("#") {
                continue;
            }
            let parts: Vec<&str> = line.split_ascii_whitespace().collect();
            for &dst_str in parts.iter() {
                let dst = dst_str.parse::<u32>()? - 1;
                graph.vertices[src as usize].edges.push(GraphEdge { dst, weight: 1 });
                if graph.vertices[src as usize].edges.len() as u32 > Self::MAX_EDGES_PER_VERTEX {
                    panic!("vertex has too many edges");
                }
            }
            src += 1;
        }
        Ok(graph)
    }
}

#[derive(Copy, Clone, Eq, PartialEq, Debug)]
pub enum NeighbourSelection {
    Local,
    Hybrid,
}

#[derive(Copy, Clone, Eq, PartialEq, Debug)]
pub enum InitialPartitioning {
    Modulo,
    Random,
    Bfs,
}

pub struct PartitioningConfig {
    pub rng_seed: u64,
    pub initial_partitioning: InitialPartitioning,
    pub max_iterations: u32,
    pub temperature_start: f32,
    pub temperature_delta: f32,
    pub alpha: f32,
    pub neighbour_sample_size: u32,
    pub random_sample_size: u32,
    pub neighbour_selection: NeighbourSelection,
}

impl Default for PartitioningConfig {
    fn default() -> Self {
        Self {
            rng_seed: 1234,
            initial_partitioning: InitialPartitioning::Bfs, // Bfs even beats the paper with 949 vs 1424 (4elt)
            max_iterations: 1000,
            temperature_start: 2.0,
            temperature_delta: 0.003,
            alpha: 2.0,
            neighbour_sample_size: 3,
            random_sample_size: 6,
            neighbour_selection: NeighbourSelection::Hybrid,
        }
    }
}

impl Graph {
    pub fn partition(&mut self, config: &PartitioningConfig, partitions: u32) {
        let mut rng = StdRng::seed_from_u64(config.rng_seed);

        match config.initial_partitioning {
            InitialPartitioning::Modulo => {
                for i in 0..self.vertices.len() {
                    self.vertices[i].color = i as u32 % partitions;
                }
            }
            InitialPartitioning::Random => {
                for i in 0..self.vertices.len() {
                    self.vertices[i].color = rng.gen_range(0..partitions);
                }
            }
            InitialPartitioning::Bfs => { // TODO maybe use this for initial partitioning: https://iccl.inf.tu-dresden.de/w/images/4/4f/Gvdb-2016.pdf
                for v in self.vertices.iter_mut() {
                    v.color = u32::MAX;
                }

                let avg_size = self.vertices.len() as u32 / partitions;

                let mut color = 0u32;
                let mut partition_size = 0u32;
                let mut queue = VecDeque::new();

                let mut visited = 0;
                while visited < self.vertices.len() {
                    queue.clear();
                    let start = self.vertices.iter().position(|v| v.color == u32::MAX).unwrap();
                    queue.push_back(start as u32);

                    while let Some(vx) = queue.pop_front() {
                        if self.vertices[vx as usize].color == u32::MAX {
                            self.vertices[vx as usize].color = color;

                            visited += 1;
                            partition_size += 1;
                            if partition_size >= avg_size {
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

                println!("bsf partitions {}", color+1);
            }
        }

        println!("initial edge cut: {}", self.calculate_edge_cut());

        let mut swaps = 0;

        let mut temperature = config.temperature_start;
        let mut candidate_buf = [0u32; Self::MAX_EDGES_PER_VERTEX as usize];

        for iteration in 0..config.max_iterations {
            for vx in 0..self.vertices.len() as u32 {
                if self.vertices[vx as usize].edges.is_empty() {
                    continue;
                }

                let mut candidate_count =
                    self.get_n_random_neighbours(&mut rng, config.neighbour_sample_size, vx, &mut candidate_buf[..]);
                //dbg!(&candidate_buf[..candidate_count as usize]);
                let mut candidate =
                    self.select_candidate(config, temperature, vx, &candidate_buf[..candidate_count as usize]);

                if candidate.is_none() && config.neighbour_selection == NeighbourSelection::Hybrid {
                    candidate_count =
                        self.get_n_random_vertices(&mut rng, config.random_sample_size, vx, &mut candidate_buf[..]);
                    candidate =
                        self.select_candidate(config, temperature, vx, &candidate_buf[..candidate_count as usize]);
                }

                if let Some(cx) = candidate {
                    self.swap_colors(vx, cx);
                    swaps += 1;
                }
            }

            // Simulated annealing.
            temperature = (temperature - config.temperature_delta).max(1.0);

            let edge_cut = self.calculate_edge_cut();
            println!("iteration: {iteration}, edge cut: {edge_cut}, swaps: {swaps}");
        }
    }

    /// Returns the sum of edge weights of neighbours having the given color.
    pub fn get_degree(&self, vx: u32, color: u32) -> u32 {
        let mut degree = 0;
        for e in self.vertices[vx as usize].edges.iter() {
            if self.vertices[e.dst as usize].color == color {
                degree += e.weight;
            }
        }
        degree
    }

    pub fn swap_colors(&mut self, va: u32, vb: u32) {
        let tmp = self.vertices[va as usize].color;
        self.vertices[va as usize].color = self.vertices[vb as usize].color;
        self.vertices[vb as usize].color = tmp;
    }

    pub fn calculate_edge_cut(&self) -> u32 {
        let mut edge_cut = 0;
        for v in self.vertices.iter() {
            for e in v.edges.iter() {
                let n = &self.vertices[e.dst as usize];
                if n.color != v.color {
                    edge_cut += e.weight;
                }
            }
        }
        edge_cut / 2
    }

    // pub fn count_partitions(&self) -> u32 {
    //     self.vertices.iter().filter(|v.pa|)
    // }

    /// Returns n random neighbours of the given vertex. If more are requested than present, duplicates may occur.
    /// The returned array is only filled with n elements!
    fn get_n_random_neighbours(&self, rng: &mut StdRng, n: u32, vx: u32, dst: &mut [u32]) -> u32 {
        let mut count = 0;
        let v = &self.vertices[vx as usize];
        let v_edge_count = v.edges.len() as u32;

        // Less edges than requested. Return all neighbours.
        if v_edge_count <= n {
            for e in v.edges.iter() {
                dst[count] = e.dst;
                count += 1;
            }
            return count as u32;
        }

        // More neighbours are present than requested. Select neighbours randomly.
        while count < n as usize {
            let r = rng.gen_range(0..v_edge_count);
            if !dst[..count].contains(&r) {
                dst[count] = r;
                count += 1;
            }
        }
        count as u32
    }

    /// Returns n random vertices of the graph excluding the given vertex.
    /// The returned array is only filled with n elements!
    fn get_n_random_vertices(&self, rng: &mut StdRng, n: u32, vx: u32, dst: &mut [u32]) -> u32 {
        let mut count = 0;
        while count < (n as usize).min(self.vertices.len()) {
            let r = rng.gen_range(0..self.vertices.len() as u32);
            if r != vx && !dst[..count].contains(&r) {
                dst[count] = r;
                count += 1;
            }
        }
        count as u32
    }

    fn select_candidate(
        &self,
        config: &PartitioningConfig,
        temperature: f32,
        vx: u32,
        candidates: &[u32],
    ) -> Option<u32> {
        let v = &self.vertices[vx as usize];
        let v_old_degree = self.get_degree(vx, v.color);

        let mut max_degree_sum = 0.0;
        let mut max_degree_candidate = None;

        for &cx in candidates.iter() {
            let c = &self.vertices[cx as usize];
            if c.color != v.color {
                let c_old_degree = self.get_degree(cx, c.color);

                let v_new_degree = self.get_degree(vx, c.color);
                let c_new_degree = self.get_degree(cx, v.color);

                let old_degree_sum =
                    (v_old_degree as f32).powf(config.alpha) + (c_old_degree as f32).powf(config.alpha);
                let new_degree_sum =
                    (v_new_degree as f32).powf(config.alpha) + (c_new_degree as f32).powf(config.alpha);

                // Higher temperatures make changes more likely.
                if new_degree_sum > max_degree_sum && new_degree_sum * temperature > old_degree_sum {
                    max_degree_sum = new_degree_sum;
                    max_degree_candidate = Some(cx);
                }
            }
        }
        max_degree_candidate
    }
}
