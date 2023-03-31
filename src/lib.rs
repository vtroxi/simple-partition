// Graph partitioning using simulated annealing based on Ja-be-Ja:
// https://www.diva-portal.org/smash/get/diva2:1043244/FULLTEXT01.pdf

use anyhow::{Context, Result};
use rand::rngs::StdRng;
use rand::{Rng, SeedableRng};
use std::collections::{HashMap, HashSet, VecDeque};
use std::io::{BufRead, BufReader, Read};
use std::ops::{Range, RangeInclusive};

fn divide_round_up(a: u32, b: u32) -> u32 {
    (a + b - 1) / b
}

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
    pub const MAX_EDGES_PER_VERTEX: u32 = 2048;

    /// Deserializes an unweighted graph in the format used by the METIS library.
    pub fn deserialize_metis<R: Read>(r: &mut R) -> Result<Self> {
        let reader = BufReader::new(r);
        let mut lines = reader.lines().map(|l| l.unwrap());

        // Parse the header
        let header = lines.next().context("could not get header line")?;
        let header_parts = header.split_ascii_whitespace().collect::<Vec<_>>();
        let vertex_count = header_parts[0]
            .parse::<usize>()
            .context("could not parse vertex count")?;
        //let edge_count = header_parts[1].parse::<usize>().context("could not parse edge count")?;

        let mut graph = Graph {
            vertices: vec![
                GraphVertex {
                    edges: vec![],
                    color: u32::MAX
                };
                vertex_count
            ],
        };

        let mut src = 0;
        for line in lines {
            if line.starts_with('%') || line.starts_with('#') {
                continue;
            }
            let parts: Vec<&str> = line.split_ascii_whitespace().collect();
            for &dst_str in parts.iter() {
                let dst = dst_str.parse::<u32>()? - 1;
                graph.vertices[src as usize]
                    .edges
                    .push(GraphEdge { dst, weight: 1 });
                if graph.vertices[src as usize].edges.len() as u32 > Self::MAX_EDGES_PER_VERTEX {
                    panic!(
                        "vertex has too many edges (the limit is {})",
                        Self::MAX_EDGES_PER_VERTEX
                    );
                }
            }
            src += 1;
        }
        Ok(graph)
    }
}

#[derive(Copy, Clone, Eq, PartialEq, Debug)]
pub enum NeighbourSelection {
    /// Just direct neighbours of vertices are selected as candidates for swapping.
    Local,
    /// If direct neighbours cannot improve the partitioning, a set number of random vertices is sampled.
    Hybrid,
}

#[derive(Copy, Clone, Eq, PartialEq, Debug)]
pub enum InitialPartitioning {
    /// Vertices are initialized to partitions based on their index in the graph (idx % partitions).
    Modulo,
    /// Vertices are initialized to random partitions.
    Random,
    /// Vertex patches get initialized based on locality using breadth first search.
    Bfs,
}

pub struct PartitioningConfig {
    /// The seed for the random number generator.
    pub rng_seed: u64,
    /// What initial partitioning method to use.
    pub initial_partitioning: InitialPartitioning,
    /// The maximum amount of iterations the algorithm is allowed to run.
    pub max_iterations: u32,
    /// If this is set to Some(n), the algorithm will stop if no changes have occurred in the last n iterations (the algorithm has converged).
    pub stop_early_n: Option<u32>,
    /// The initial temperature for the simulated annealing.
    pub temperature_start: f32,
    /// How much the temperature decreases every round.
    pub temperature_delta: f32,
    /// The alpha parameter used to decide if a swap of two vertices is beneficial (see paper).
    pub alpha: f32,
    /// How many random neighbours should be sampled for each vertex.
    pub neighbour_sample_size: u32,
    /// How many random vertices should be sampled if no swaps with neighbours can improve the partitioning.
    pub random_sample_size: u32,
    /// The way that neighbours are selected.
    pub neighbour_selection: NeighbourSelection,
}

impl Default for PartitioningConfig {
    fn default() -> Self {
        Self {
            rng_seed: 1234,
            initial_partitioning: InitialPartitioning::Random,
            max_iterations: 500,
            stop_early_n: Some(4),
            temperature_start: 2.0,
            temperature_delta: 0.003,
            alpha: 2.0,
            neighbour_sample_size: 2,
            random_sample_size: 16,
            neighbour_selection: NeighbourSelection::Hybrid,
        }
    }
}

impl Graph {
    // General functions and metrics

    /// Counts the number of distinct partitions in the graph.
    pub fn count_partitions(&self) -> u32 {
        self.vertices
            .iter()
            .map(|v| v.color)
            .collect::<HashSet<_>>()
            .len() as u32
    }

    /// Calculates the min, max and average partition size.
    pub fn partition_sizes(&self) -> (u32, u32, u32) {
        let mut counts = vec![0u32; self.count_partitions() as usize];
        for v in self.vertices.iter() {
            counts[v.color as usize] += 1;
        }
        let min = *counts.iter().min().unwrap_or(&0);
        let max = *counts.iter().max().unwrap_or(&0);
        let avg = counts.iter().sum::<u32>() / counts.len() as u32;
        (min, max, avg)
    }

    /// Returns the summed weights of all cut edges.
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

    // Partitioning

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
            InitialPartitioning::Bfs => {
                // TODO maybe use this for initial partitioning: https://iccl.inf.tu-dresden.de/w/images/4/4f/Gvdb-2016.pdf
                for v in self.vertices.iter_mut() {
                    v.color = u32::MAX;
                }

                let target_size = divide_round_up(self.vertices.len() as u32, partitions);

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

        log::trace!("initial edge cut: {}", self.calculate_edge_cut());

        let mut stop_indicators = vec![0; config.stop_early_n.unwrap_or(1) as usize];

        let mut swaps = 0;

        let mut temperature = config.temperature_start;
        let mut candidate_buf = vec![0u32; Self::MAX_EDGES_PER_VERTEX as usize];

        for iteration in 0..config.max_iterations {
            for vx in 0..self.vertices.len() as u32 {
                if self.vertices[vx as usize].edges.is_empty() {
                    continue;
                }

                let mut candidate_count = self.get_n_random_neighbours(
                    &mut rng,
                    config.neighbour_sample_size,
                    vx,
                    &mut candidate_buf[..],
                );
                //dbg!(&candidate_buf[..candidate_count as usize]);
                let mut candidate = self.select_candidate(
                    config,
                    temperature,
                    vx,
                    &candidate_buf[..candidate_count as usize],
                );

                if candidate.is_none() && config.neighbour_selection == NeighbourSelection::Hybrid {
                    candidate_count = self.get_n_random_vertices(
                        &mut rng,
                        config.random_sample_size,
                        vx,
                        &mut candidate_buf[..],
                    );
                    candidate = self.select_candidate(
                        config,
                        temperature,
                        vx,
                        &candidate_buf[..candidate_count as usize],
                    );
                }

                if let Some(cx) = candidate {
                    self.swap_colors(vx, cx);
                    swaps += 1;
                }
            }

            // Simulated annealing.
            temperature = (temperature - config.temperature_delta).max(1.0);

            let edge_cut = self.calculate_edge_cut();
            log::trace!("iteration: {iteration}, edge cut: {edge_cut}, swaps: {swaps}");

            if config.stop_early_n.is_some() {
                stop_indicators.rotate_right(1);
                stop_indicators[0] = swaps;
                // Stop if no swaps occurred in the last n rounds (all indicators are equal).
                if stop_indicators.iter().min() == stop_indicators.iter().max() {
                    log::trace!("stopping early");
                    break;
                }
            }
        }
    }

    fn swap_colors(&mut self, va: u32, vb: u32) {
        let tmp = self.vertices[va as usize].color;
        self.vertices[va as usize].color = self.vertices[vb as usize].color;
        self.vertices[vb as usize].color = tmp;
    }

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

                let old_degree_sum = (v_old_degree as f32).powf(config.alpha)
                    + (c_old_degree as f32).powf(config.alpha);
                let new_degree_sum = (v_new_degree as f32).powf(config.alpha)
                    + (c_new_degree as f32).powf(config.alpha);

                // Higher temperatures make changes more likely.
                if new_degree_sum > max_degree_sum && new_degree_sum * temperature > old_degree_sum
                {
                    max_degree_sum = new_degree_sum;
                    max_degree_candidate = Some(cx);
                }
            }
        }
        max_degree_candidate
    }
}
