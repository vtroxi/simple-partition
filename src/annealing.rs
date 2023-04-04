use crate::{Graph, InitialPartitioningMethod};
use rand::rngs::StdRng;
use rand::{Rng, SeedableRng};

#[derive(Copy, Clone, Eq, PartialEq, Debug)]
pub enum NeighbourSelection {
    /// Just direct neighbours of vertices are selected as candidates for swapping.
    Local,
    /// If direct neighbours cannot improve the partitioning, a set number of random vertices is sampled.
    Hybrid,
}

pub struct AnnealingPartitioningConfig {
    /// The seed for the random number generator.
    pub rng_seed: u64,
    /// What initial partitioning method to use.
    pub initial_partitioning: InitialPartitioningMethod,
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

impl Default for AnnealingPartitioningConfig {
    fn default() -> Self {
        Self {
            rng_seed: 1234,
            initial_partitioning: InitialPartitioningMethod::Random,
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
    /// Splits the graph into a given number of partitions while minimizing the edge cut cost.
    /// This algorithm uses simulated annealing based on the [Ja-be-Ja](https://www.diva-portal.org/smash/get/diva2:1043244/FULLTEXT01.pdf) paper.
    pub fn partition_annealing(&mut self, config: &AnnealingPartitioningConfig, partitions: u32) {
        let mut rng = StdRng::seed_from_u64(config.rng_seed);

        self.partition_initial(config.initial_partitioning, &mut rng, partitions);
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
        config: &AnnealingPartitioningConfig,
        temperature: f32,
        vx: u32,
        candidates: &[u32],
    ) -> Option<u32> {
        let v = &self.vertices[vx as usize];
        let v_old_degree = self.get_internal_degree(vx, v.color);

        let mut max_degree_sum = 0.0;
        let mut max_degree_candidate = None;

        for &cx in candidates.iter() {
            let c = &self.vertices[cx as usize];
            if c.color != v.color {
                let c_old_degree = self.get_internal_degree(cx, c.color);

                let v_new_degree = self.get_internal_degree(vx, c.color);
                let c_new_degree = self.get_internal_degree(cx, v.color);

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
