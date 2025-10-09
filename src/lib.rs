//! # Multi-Order Runge-Kutta methods (MORK)
//!
//! `MORK` is an implementation of general multi-order Runge-Kutta methods ([GMORK]) and Runge-Kutta methods ([RK]), as described in the paper ["Multi-order Runge-Kutta methods or how to numerically solve initial value problems of any order"](https://doi.org/10.48550/arXiv.2509.23513). There also is an implementation of node-determined multi-order Runge Kutta methods ([NDMORK])
//! 
//! All class of methods implement the [Solver] trait, which consists of an [approximate][Solver::approximate] function. This function, given a differential equation function, initial instant, initial values, and a step size, returns the approximation of the method/struct which implements this trait. 
//! 
//! Some Runge-Kutta methods ([RK::list]) and node-determined multi-order Runge-Kutta methods ([NDMORK::list]) are already implemented.

#![allow(non_snake_case)]
pub mod GMORK;
pub mod NDMORK;
pub mod RK;
pub mod graph;

use crate::graph::*;

const ERROR_FRACTION: f64 = 0.001;
const MAX_ITER: u32 = 500;
const MIN_ITER: u32 = 50;

/// [Solver] is the trait that numerical schemes implement.
pub trait Solver {
    /// Given a differential equation function, initial instant, initial values, and a step size, [approximate][Solver::approximate] returns the approximation of the method/struct which implements this trait.
    fn approximate(
        &mut self,
        t0: f64,
        h: f64,
        f: &dyn Fn(f64, &Vec<Vec<f64>>) -> Vec<f64>,
        y0: &Vec<Vec<f64>>,
    ) -> Vec<Vec<f64>>;
}

/// [enum@SCC] allows to distinguish between implicit and explicit strongly connected components.
#[derive(Debug, Clone)]
pub enum SCC {
    Implicit(Vec<usize>, Vec<usize>), // J and [|1,s|] without J
    Explicit(usize),
}

/// [create_computation_order] takes the maximum weight digraph of a method and outputs an order of computation.
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