# simple-partition
### A simple graph partitioning library for Rust
#### Features
- Direct k-way partitioning using simulated annealing based on the paper [Ja-be-Ja](https://www.diva-portal.org/smash/get/diva2:1043244/FULLTEXT01.pdf) - "A Distributed Algorithm for
  Balanced Graph Partitioning".
- Deserialization of graphs in the METIS format
- Graph utility functions and partitioning metrics

#### Note
If you're looking for higher quality partitioning using multilevel techniques, you might be interested in my cross-platform [METIS wrapper for Rust](https://github.com/vtroxi/metis-rs)