#![allow(non_snake_case)]

use std::{sync::Condvar, thread::JoinHandle};

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
pub enum Block {
    Implicit(Vec<usize>, Vec<usize>), // J and [|1,s|] \ J
    Explicit(usize),
}

pub fn create_queue(weight_graph: &Vec<Vec<bool>>) -> Vec<Block> {
    let SCC = scc(weight_graph);
    topological_sort(&contraction(weight_graph, &SCC))
        .into_iter()
        .map(|i| {
            if SCC[i].len() > 1 || weight_graph[SCC[i][0]][SCC[i][0]] {
                let J = SCC[i].clone();
                let comp_J = (0..(weight_graph.len() - 1))
                    .filter(|j| !J.contains(j))
                    .collect();
                Block::Implicit(J, comp_J)
            } else {
                Block::Explicit(SCC[i][0])
            }
        })
        .collect()
}

pub struct GMORK {
    pub nodes: Vec<f64>,
    pub main_weights: Vec<Vec<Vec<f64>>>,
    main_weights_f: Box<dyn Fn(u32) -> Vec<Vec<f64>>>,
    pub secondary_weights: Vec<Vec<Vec<f64>>>,
    secondary_weights_f: Box<dyn Fn(u32) -> Vec<Vec<f64>>>,
    pub queue: Vec<Block>,
    pub s: usize,
    pub stored_length: usize,
    pub cyclic_derivatives: Vec<Vec<bool>>, // [N-1][j]
    pub factorial: Vec<f64>,
    pub h_powers: Vec<f64>,
    pub h: f64,
    pub error_fraction: f64,
    pub min_iter: u32,
    pub max_iter: u32,
}

impl GMORK {
    pub fn new(
        main_weights_function: Box<dyn Fn(u32) -> Vec<Vec<f64>>>,
        secondary_weights_function: Box<dyn Fn(u32) -> Vec<Vec<f64>>>,
        nodes: Vec<f64>,
        maximum_weight_graph: Vec<Vec<bool>>,
    ) -> Self {
        let main_weights = vec![main_weights_function(1)];
        let secondary_weights = vec![secondary_weights_function(1)];
        let s = main_weights[0].len() - 1;
        let mut cyclic_derivatives: Vec<Vec<bool>> = vec![vec![false; s]];
        let queue = create_queue(&maximum_weight_graph);
        for task in queue.iter() {
            if let Block::Implicit(J, _) = task {
                for &j in J {
                    for &j1 in J {
                        if main_weights[0][j][j1] != 0. {
                            cyclic_derivatives[0][j] = true;
                            break;
                        }
                    }
                }
            }
        }
        GMORK {
            s,
            stored_length: 1,
            nodes,
            main_weights,
            main_weights_f: main_weights_function,
            secondary_weights,
            secondary_weights_f: secondary_weights_function,
            factorial: vec![1., 1.],
            h: 0.,
            h_powers: vec![1., 0.],
            queue,
            cyclic_derivatives,
            error_fraction: ERROR_FRACTION,
            min_iter: MIN_ITER,
            max_iter: MAX_ITER,
        }
    }

    pub fn set_minimum_length(
        n: usize,
        s: usize,
        stored_length: &mut usize,
        factorial: &mut Vec<f64>,
        main_weights: &mut Vec<Vec<Vec<f64>>>,
        main_weights_f: &dyn Fn(u32) -> Vec<Vec<f64>>,
        secondary_weights: &mut Vec<Vec<Vec<f64>>>,
        secondary_weights_f: &dyn Fn(u32) -> Vec<Vec<f64>>,
        h_powers: &mut Vec<f64>,
        h: f64,
        cyclic_derivatives: &mut Vec<Vec<bool>>,
        queue: &Vec<Block>,
    ) {
        if *stored_length >= n {
            return;
        }
        factorial.extend(vec![0.; n - *stored_length]);
        main_weights.extend((*stored_length..n).map(|N| (main_weights_f)(N as u32 + 1)));
        secondary_weights.extend((*stored_length..n).map(|N| (secondary_weights_f)(N as u32 + 1)));
        h_powers.extend((*stored_length..n).map(|N| h.powi(N as i32 + 1)));
        cyclic_derivatives.extend((*stored_length..n).map(|_| vec![false; s]));
        for N in *stored_length..n {
            factorial[N + 1] = factorial[N] * (N as f64 + 1.);
            for task in queue.iter() {
                if let Block::Implicit(J, _) = task {
                    for &j in J {
                        for &j1 in J {
                            if main_weights[N][j][j1] != 0. {
                                cyclic_derivatives[N][j] = true;
                                break;
                            }
                        }
                    }
                }
            }
        }
        *stored_length = n;
    }

    pub fn add_constant_part(
        J: &Vec<usize>,
        j: usize,
        y: &mut Vec<Vec<Vec<f64>>>,
        y0: &Vec<Vec<f64>>,
        F: &Vec<Vec<f64>>,
        main_weights: &Vec<Vec<Vec<f64>>>,
        secondary_weights: &Vec<Vec<Vec<f64>>>,
        h_powers: &Vec<f64>,
        factorial: &Vec<f64>,
    ) {
        let mut sum;
        for k in 0..y0.len() {
            for N in 0..y0[k].len() {
                for N1 in 1..=N {
                    y[j][k][N] += secondary_weights[N][N1][j] * h_powers[N1] * y0[k][N - N1]
                }
                sum = 0.;
                for &j1 in J {
                    sum += main_weights[N][j][j1] * F[j1][k];
                }
                y[j][k][N] += h_powers[N + 1] / factorial[N + 1] * sum;
            }
        }
    }

    pub fn picard_iterations(
        t: f64,
        h: f64,
        y0: &Vec<Vec<f64>>,
        y: &mut Vec<Vec<Vec<f64>>>,
        F: &mut Vec<Vec<f64>>,
        f: &dyn Fn(f64, &Vec<Vec<f64>>) -> Vec<f64>,
        J: &Vec<usize>,
        threshold: f64,
        min_iter: u32,
        max_iter: u32,
        nodes: &Vec<f64>,
        cyclic_derivatives: &Vec<Vec<bool>>,
        main_weights: &Vec<Vec<Vec<f64>>>,
        h_powers: &Vec<f64>,
        factorial: &Vec<f64>,
    ) {
        let mut sum;
        let constant = y.clone();
        let mut iter_count = 0;
        let mut d = threshold + 1.;
        let mut f_cache: Vec<f64>;
        while iter_count < min_iter || (d > threshold && iter_count < max_iter) {
            iter_count += 1;
            // update evaluations and calculate difference
            d = 0.;
            for &j in J {
                f_cache = f(t + nodes[j] * h, &y[j]);
                // calculate difference for threshold
                for k in 0..f_cache.len() {
                    d = d.max((f_cache[k] - F[j][k]).abs());
                }
                F[j] = f_cache;
            }
            // constant part
            *y = constant.clone();
            // add evaluations and calculate norm difference
            for &j in J {
                for k in 0..y0.len() {
                    for N in (0..y0[k].len()).filter(|&N| cyclic_derivatives[N][j]) {
                        sum = 0.;
                        for &j1 in J {
                            sum += main_weights[N][j][j1] * F[j1][k];
                        }
                        y[j][k][N] += h_powers[N + 1] / factorial[N + 1] * sum;
                    }
                }
            }
        }
    }
    pub fn approximate_GMORK(
        t: f64,
        h: f64,
        f: &dyn Fn(f64, &Vec<Vec<f64>>) -> Vec<f64>,
        y0: &Vec<Vec<f64>>,
        threshold: f64,
        s: usize,
        queue: &Vec<Block>,
        nodes: &Vec<f64>,
        main_weights: &Vec<Vec<Vec<f64>>>,
        secondary_weights: &Vec<Vec<Vec<f64>>>,
        h_powers: &Vec<f64>,
        factorial: &Vec<f64>,
        min_iter: u32,
        max_iter: u32,
        cyclic_derivatives: &Vec<Vec<bool>>,
    ) -> Vec<Vec<f64>> {
        let mut F: Vec<Vec<f64>> = (0..s).map(|_| vec![0.; y0.len()]).collect();
        let mut y: Vec<Vec<Vec<f64>>> = (0..=s).map(|_| y0.clone()).collect();
        for task in queue.iter() {
            match task {
                Block::Explicit(j) => {
                    let j = *j;
                    GMORK::add_constant_part(
                        &(0..s).collect(),
                        j,
                        &mut y,
                        y0,
                        &F,
                        main_weights,
                        secondary_weights,
                        h_powers,
                        factorial,
                    );
                    if j != s {
                        F[j] = f(t + nodes[j] * h, &y[j]);
                    }
                }
                Block::Implicit(J, comp_J) => {
                    // calculate constant terms
                    for &j in J {
                        GMORK::add_constant_part(
                            comp_J,
                            j,
                            &mut y,
                            y0,
                            &F,
                            main_weights,
                            secondary_weights,
                            h_powers,
                            factorial,
                        );
                    }
                    GMORK::picard_iterations(
                        t,
                        h,
                        y0,
                        &mut y,
                        &mut F,
                        f,
                        J,
                        threshold,
                        min_iter,
                        max_iter,
                        nodes,
                        cyclic_derivatives,
                        main_weights,
                        h_powers,
                        factorial,
                    );
                }
            }
        }
        return y[s].clone();
    }
}

impl Solver for GMORK {
    fn approximate(
        &mut self,
        t: f64,
        h: f64,
        f: &dyn Fn(f64, &Vec<Vec<f64>>) -> Vec<f64>,
        y0: &Vec<Vec<f64>>,
    ) -> Vec<Vec<f64>> {
        if h != self.h {
            self.h = h;
            for N in 1..=self.stored_length {
                self.h_powers[N] = h.powi(N as i32)
            }
        }
        // calculate difference threshold for picard iterations
        let mut threshold = y0[0][0].abs();
        for k in 0..y0.len() {
            // verify the length of the method is enough
            if y0[k].len() > self.stored_length {
                GMORK::set_minimum_length(
                    y0[k].len(),
                    self.s,
                    &mut self.stored_length,
                    &mut self.factorial,
                    &mut self.main_weights,
                    &self.main_weights_f,
                    &mut self.secondary_weights,
                    &self.secondary_weights_f,
                    &mut self.h_powers,
                    h,
                    &mut self.cyclic_derivatives,
                    &self.queue,
                );
            }
            for N in 0..y0[k].len() {
                if threshold < y0[k][N].abs() {
                    threshold = y0[k][N].abs()
                }
            }
        }
        threshold *= self.error_fraction;
        GMORK::approximate_GMORK(
            t,
            h,
            f,
            y0,
            threshold,
            self.s,
            &self.queue,
            &self.nodes,
            &self.main_weights,
            &self.secondary_weights,
            &self.h_powers,
            &self.factorial,
            self.min_iter,
            self.max_iter,
            &self.cyclic_derivatives,
        )
    }
}
pub struct NDMORK {
    pub s: usize,
    pub stored_length: usize,
    pub nodes: Vec<f64>,
    pub weights: Vec<Vec<Vec<f64>>>,
    pub weights_function: Box<dyn Fn(u32) -> Vec<Vec<f64>>>,
    pub factorial: Vec<f64>,
    pub coefficients: Vec<Vec<f64>>,
    pub h: f64,
    pub h_powers: Vec<f64>,
    pub queue: Vec<Block>,
    pub cyclic_derivatives: Vec<Vec<bool>>, // [N-1][j]
    pub error_fraction: f64,
    pub min_iter: u32,
    pub max_iter: u32,
}

impl NDMORK {
    pub fn new(
        weights_function: Box<dyn Fn(u32) -> Vec<Vec<f64>>>,
        nodes: Vec<f64>,
        maximum_weight_graph: Vec<Vec<bool>>,
    ) -> Self {
        let s = nodes.len() - 1;
        let weights = vec![weights_function(1)];
        let queue = create_queue(&maximum_weight_graph);
        let mut cyclic_derivatives: Vec<Vec<bool>> = vec![vec![false; s]];
        for task in queue.iter() {
            if let Block::Implicit(J, _) = task {
                for &j in J {
                    for &j1 in J {
                        if weights[0][j][j1] != 0. {
                            cyclic_derivatives[0][j] = true;
                            break;
                        }
                    }
                }
            }
        }
        NDMORK {
            s,
            stored_length: 1,
            nodes,
            weights,
            weights_function,
            factorial: vec![1., 1.],
            coefficients: vec![vec![1.; s + 1]],
            h: 0.,
            h_powers: vec![1., 0.],
            queue,
            cyclic_derivatives,
            error_fraction: ERROR_FRACTION,
            min_iter: MIN_ITER,
            max_iter: MAX_ITER,
        }
    }

    pub fn set_minimum_length(
        h_powers: &mut Vec<f64>,
        weights_function: &dyn Fn(u32) -> Vec<Vec<f64>>,
        weights: &mut Vec<Vec<Vec<f64>>>,
        h: f64,
        s: usize,
        factorial: &mut Vec<f64>,
        stored_length: &mut usize,
        n: usize,
        cyclic_derivatives: &mut Vec<Vec<bool>>,
        queue: &Vec<Block>,
        coefficients: &mut Vec<Vec<f64>>,
        nodes: &Vec<f64>,
    ) {
        if *stored_length >= n {
            return;
        }
        factorial.extend(vec![0.; n - *stored_length]);
        weights.extend((*stored_length..n).map(|N| (weights_function)(N as u32 + 1)));
        h_powers.extend((*stored_length..n).map(|N| h.powi(N as i32 + 1)));
        cyclic_derivatives.extend((*stored_length..n).map(|_| vec![false; s]));
        for N in *stored_length..n {
            factorial[N + 1] = factorial[N] * (N as f64 + 1.);
            for task in queue.iter() {
                if let Block::Implicit(J, _) = task {
                    for &j in J {
                        for &j1 in J {
                            if weights[N][j][j1] != 0. {
                                cyclic_derivatives[N][j] = true;
                                break;
                            }
                        }
                    }
                }
            }
        }
        coefficients.extend((*stored_length..n).map(|N| {
            (0..s + 1)
                .map(|j| nodes[j].powi(N as i32) / factorial[N])
                .collect()
        }));
        *stored_length = n;
    }

    fn add_constant_part(
        J: &Vec<usize>,
        j: usize,
        y: &mut Vec<Vec<Vec<f64>>>,
        y0: &Vec<Vec<f64>>,
        F: &Vec<Vec<f64>>,
        nodes: &Vec<f64>,
        weights: &Vec<Vec<Vec<f64>>>,
        coefficients: &Vec<Vec<f64>>,
        h_powers: &Vec<f64>,
        factorial: &Vec<f64>,
    ) {
        let mut sum;
        for k in 0..y0.len() {
            for N in 0..y0[k].len() {
                if nodes[j] != 0. {
                    for N1 in 1..=N {
                        y[j][k][N] += coefficients[N1][j] * h_powers[N1] * y0[k][N - N1]
                    }
                }
                sum = 0.;
                for &j1 in J {
                    sum += weights[N][j][j1] * F[j1][k];
                }
                y[j][k][N] += h_powers[N + 1] / factorial[N + 1] * sum;
            }
        }
    }

    fn picard_iterations(
        t: f64,
        h: f64,
        y0: &Vec<Vec<f64>>,
        y: &mut Vec<Vec<Vec<f64>>>,
        F: &mut Vec<Vec<f64>>,
        f: &dyn Fn(f64, &Vec<Vec<f64>>) -> Vec<f64>,
        J: &Vec<usize>,
        threshold: f64,
        min_iter: u32,
        max_iter: u32,
        nodes: &Vec<f64>,
        cyclic_derivatives: &Vec<Vec<bool>>,
        weights: &Vec<Vec<Vec<f64>>>,
        h_powers: &Vec<f64>,
        factorial: &Vec<f64>,
    ) {
        let constant = y.clone();
        // Picard iterations
        let mut iter_count = 0;
        let mut d = threshold + 1.;
        let mut f_cache: Vec<f64>;
        let mut sum;
        while iter_count < min_iter || (d > threshold && iter_count < max_iter) {
            iter_count += 1;
            // update evaluations and calculate difference
            d = 0.;
            for &j in J {
                f_cache = f(t + nodes[j] * h, &y[j]);
                // calculate difference for threshold
                for k in 0..f_cache.len() {
                    d = d.max((f_cache[k] - F[j][k]).abs());
                }
                F[j] = f_cache;
            }
            // add evaluations
            for &j in J {
                for k in 0..y0.len() {
                    for N in (0..y0[k].len()).filter(|&N| cyclic_derivatives[N][j]) {
                        sum = 0.;
                        for &j1 in J {
                            sum += weights[N][j][j1] * F[j1][k];
                        }
                        y[j][k][N] = constant[j][k][N] + h_powers[N + 1] / factorial[N + 1] * sum;
                    }
                }
            }
        }
    }

    pub fn approximate_ND(
        t: f64,
        h: f64,
        f: &dyn Fn(f64, &Vec<Vec<f64>>) -> Vec<f64>,
        y0: &Vec<Vec<f64>>,
        threshold: f64,
        min_iter: u32,
        max_iter: u32,
        s: usize,
        queue: &Vec<Block>,
        cyclic_derivatives: &Vec<Vec<bool>>,
        nodes: &Vec<f64>,
        weights: &Vec<Vec<Vec<f64>>>,
        coefficients: &Vec<Vec<f64>>,
        h_powers: &Vec<f64>,
        factorial: &Vec<f64>,
    ) -> Vec<Vec<f64>> {
        let mut F: Vec<Vec<f64>> = (0..s).map(|_| vec![0.; y0.len()]).collect();
        let mut y: Vec<Vec<Vec<f64>>> = (0..=s).map(|_| y0.clone()).collect();
        let J_explicit = (0..s).collect();
        for task in queue.iter() {
            match task {
                Block::Explicit(j) => {
                    let j = *j;
                    NDMORK::add_constant_part(
                        &J_explicit,
                        j,
                        &mut y,
                        y0,
                        &F,
                        nodes,
                        weights,
                        coefficients,
                        h_powers,
                        factorial,
                    );
                    if j != s {
                        F[j] = f(t + nodes[j] * h, &y[j]);
                    }
                }
                Block::Implicit(J, comp_J) => {
                    for &j in J {
                        NDMORK::add_constant_part(
                            comp_J,
                            j,
                            &mut y,
                            y0,
                            &F,
                            nodes,
                            weights,
                            coefficients,
                            h_powers,
                            factorial,
                        );
                    }
                    NDMORK::picard_iterations(
                        t,
                        h,
                        y0,
                        &mut y,
                        &mut F,
                        f,
                        J,
                        threshold,
                        min_iter,
                        max_iter,
                        nodes,
                        cyclic_derivatives,
                        weights,
                        h_powers,
                        factorial,
                    );
                }
            }
        }
        return y[s].clone();
    }
}

impl Solver for NDMORK {
    fn approximate(
        &mut self,
        t: f64,
        h: f64,
        f: &dyn Fn(f64, &Vec<Vec<f64>>) -> Vec<f64>,
        y0: &Vec<Vec<f64>>,
    ) -> Vec<Vec<f64>> {
        if h == 0. {
            return y0.clone();
        }
        if h != self.h {
            self.h = h;
            for N in 1..=self.stored_length {
                self.h_powers[N] = h.powi(N as i32)
            }
        }
        // calculate difference threshold for picard iterations
        let mut threshold = y0[0][0].abs();
        for k in 0..y0.len() {
            // verify the length of the method is enough
            // /*
            if y0[k].len() > self.stored_length {
                NDMORK::set_minimum_length(
                    &mut self.h_powers,
                    &self.weights_function,
                    &mut self.weights,
                    h,
                    self.s,
                    &mut self.factorial,
                    &mut self.stored_length,
                    y0[k].len(),
                    &mut self.cyclic_derivatives,
                    &self.queue,
                    &mut self.coefficients,
                    &self.nodes,
                );
            }
            // */
            for N in 0..y0[k].len() {
                if threshold < y0[k][N].abs() {
                    threshold = y0[k][N].abs()
                }
            }
        }
        threshold *= self.error_fraction;
        NDMORK::approximate_ND(
            t,
            h,
            f,
            y0,
            threshold,
            self.min_iter,
            self.max_iter,
            self.s,
            &self.queue,
            &self.cyclic_derivatives,
            &self.nodes,
            &self.weights,
            &self.coefficients,
            &self.h_powers,
            &self.factorial,
        )
    }
}

pub struct RK {
    pub s: usize,
    pub nodes: Vec<f64>,
    pub weights: Vec<Vec<f64>>,
    pub queue: Vec<Block>,
    pub error_fraction: f64,
    pub min_iter: u32,
    pub max_iter: u32,
}

impl RK {
    pub fn new(weights: Vec<Vec<f64>>, nodes: Vec<f64>) -> Self {
        let s = weights.len() - 1;
        let weight_graph = (0..=s)
            .map(|j| {
                (0..=s)
                    .map(|j1| {
                        if j1 == s || weights[j][j1] == 0. {
                            false
                        } else {
                            true
                        }
                    })
                    .collect()
            })
            .collect();
        let queue = create_queue(&weight_graph);
        let s = nodes.len() - 1;
        RK {
            s,
            nodes,
            weights,
            queue,
            error_fraction: ERROR_FRACTION,
            min_iter: MIN_ITER,
            max_iter: MAX_ITER,
        }
    }

    pub fn add_constant_part(
        J: &Vec<usize>,
        j: usize,
        y: &mut Vec<Vec<Vec<f64>>>,
        y0: &Vec<Vec<f64>>,
        F: &Vec<Vec<f64>>,
        h: f64,
        weights: &Vec<Vec<f64>>,
    ) {
        let mut sum;
        for k in 0..y0.len() {
            sum = 0.;
            for &j1 in J {
                sum += weights[j][j1] * F[j1][k];
            }
            y[j][k][0] += h * sum;
            for N in 1..y0[k].len() {
                sum = 0.;
                for &j1 in J {
                    sum += weights[j][j1] * y[j1][k][N - 1];
                }
                y[j][k][N] += h * sum;
            }
        }
    }

    pub fn picard_iterations(
        threshold: f64,
        min_iter: u32,
        max_iter: u32,
        t: f64,
        h: f64,
        nodes: &Vec<f64>,
        F: &mut Vec<Vec<f64>>,
        J: &Vec<usize>,
        f: &dyn Fn(f64, &Vec<Vec<f64>>) -> Vec<f64>,
        weights: &Vec<Vec<f64>>,
        y0: &Vec<Vec<f64>>,
        y: &mut Vec<Vec<Vec<f64>>>,
    ) {
        let constant = y.clone();
        let mut iter_count = 0;
        let mut d = threshold + 1.;
        let mut f_cache;
        let mut sum;
        while iter_count < min_iter || (d > threshold && iter_count < max_iter) {
            iter_count += 1;
            // update evaluations and calculate difference
            d = 0.;
            for &j in J {
                f_cache = f(t + nodes[j] * h, &y[j]);
                // calculate difference for threshold
                for k in 0..f_cache.len() {
                    d = d.max((f_cache[k] - F[j][k]).abs());
                }
                F[j] = f_cache;
            }
            // add evaluations
            for &j in J {
                for k in 0..y0.len() {
                    for N in (1..y0[k].len()).rev() {
                        sum = 0.;
                        for &j1 in J {
                            sum += weights[j][j1] * y[j1][k][N - 1];
                        }
                        y[j][k][N] = constant[j][k][N] + h * sum;
                    }
                    sum = 0.;
                    for &j1 in J {
                        sum += weights[j][j1] * F[j1][k];
                    }
                    y[j][k][0] = constant[j][k][0] + h * sum;
                }
            }
        }
    }

    pub fn approximate_RK(
        t: f64,
        h: f64,
        f: &dyn Fn(f64, &Vec<Vec<f64>>) -> Vec<f64>,
        y0: &Vec<Vec<f64>>,
        s: usize,
        threshold: f64,
        queue: &Vec<Block>,
        weights: &Vec<Vec<f64>>,
        nodes: &Vec<f64>,
        min_iter: u32,
        max_iter: u32,
    ) -> Vec<Vec<f64>> {
        let mut F: Vec<Vec<f64>> = (0..s).map(|_| vec![0.; y0.len()]).collect();
        let mut y: Vec<Vec<Vec<f64>>> = (0..=s).map(|_| y0.clone()).collect();
        for task in queue.iter() {
            match task {
                Block::Explicit(j) => {
                    let j = *j;
                    RK::add_constant_part(&(0..s).collect(), j, &mut y, y0, &F, h, weights);
                    if j != s {
                        F[j] = f(t + nodes[j] * h, &y[j]);
                    }
                }
                Block::Implicit(J, comp_J) => {
                    // calculate constant terms
                    for &j in J {
                        RK::add_constant_part(comp_J, j, &mut y, y0, &F, h, weights);
                    }
                    RK::picard_iterations(
                        threshold, min_iter, max_iter, t, h, nodes, &mut F, J, f, weights, y0,
                        &mut y,
                    );
                }
            }
        }
        return y[s].clone();
    }
}

impl Solver for RK {
    fn approximate(
        &mut self,
        t: f64,
        h: f64,
        f: &dyn Fn(f64, &Vec<Vec<f64>>) -> Vec<f64>,
        y0: &Vec<Vec<f64>>,
    ) -> Vec<Vec<f64>> {
        if h == 0. {
            return y0.clone();
        }
        // calculate difference threshold for picard iterations
        let mut threshold = y0[0][0].abs();
        for k in 0..y0.len() {
            for N in 0..y0[k].len() {
                if threshold < y0[k][N].abs() {
                    threshold = y0[k][N].abs()
                }
            }
        }
        threshold *= self.error_fraction;
        RK::approximate_RK(
            t,
            h,
            f,
            y0,
            self.s,
            threshold,
            &self.queue,
            &self.weights,
            &self.nodes,
            self.min_iter,
            self.max_iter,
        )
    }
}
