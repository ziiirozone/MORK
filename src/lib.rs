#![allow(non_snake_case)]
pub mod graph;
pub mod RK;
pub mod NDMORK;
pub mod GMORK;

use crate::graph::*;

const ERROR_FRACTION: f64 = 0.001;
const MAX_ITER: u32 = 500;
const MIN_ITER: u32 = 10;

pub trait Solver {
    fn approximate(
        &mut self,
        t0: f64,
        h: f64,
        f: &dyn Fn(f64, &Vec<Vec<f64>>) -> Vec<f64>,
        y0: &Vec<Vec<f64>>,
    ) -> Vec<Vec<f64>>;
}

#[derive(Debug, Clone)]
pub enum SCC {
    Implicit(Vec<usize>, Vec<usize>), // J and [|1,s|] \ J
    Explicit(usize),
}

pub fn create_computation_order(weight_graph: &Vec<Vec<bool>>) -> Vec<SCC> {
    let SCC = SCC(weight_graph);
    topological_sort(&contraction(weight_graph, &SCC))
        .into_iter()
        .map(|i| {
            if SCC[i].len() > 1 || weight_graph[SCC[i][0]][SCC[i][0]] {
                let J = SCC[i].clone();
                let comp_J = (0..(weight_graph.len() - 1))
                    .filter(|j| !J.contains(j))
                    .collect();
                SCC::Implicit(J, comp_J)
            } else {
                SCC::Explicit(SCC[i][0])
            }
        })
        .collect()
}
