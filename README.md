# simple-partition
### A simple graph partitioning library for Rust
#### Features/Algorithms
- Direct k-way partitioning using simulated annealing based on the paper [Ja-be-Ja](https://www.diva-portal.org/smash/get/diva2:1043244/FULLTEXT01.pdf) - "A Distributed Algorithm for
  Balanced Graph Partitioning".
- [TODO] Recursive bisection with strict partition bounds based on an improved variation of the Kernighan-Lin algorithm: [Quick Cut](https://www.researchgate.net/publication/2505803_New_Faster_Kernighan-Lin-Type_Graph-Partitioning_Algorithms) - "New Faster Kernighan-Lin-Type Graph-Partitioning Algorithms"
- Deserialization of graphs in the METIS format