use crate::*;
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
    pub computation_order: Vec<SCC>,
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
        let computation_order = create_computation_order(&maximum_weight_graph);
        let mut cyclic_derivatives: Vec<Vec<bool>> = vec![vec![false; s]];
        for task in computation_order.iter() {
            if let SCC::Implicit(J, _) = task {
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
            computation_order,
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
        computation_order: &Vec<SCC>,
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
            for task in computation_order.iter() {
                if let SCC::Implicit(J, _) = task {
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
        computation_order: &Vec<SCC>,
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
        for task in computation_order.iter() {
            match task {
                SCC::Explicit(j) => {
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
                SCC::Implicit(J, comp_J) => {
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
        let mut threshold = y0[0][0].abs();
        for k in 0..y0.len() {
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
                    &self.computation_order,
                    &mut self.coefficients,
                    &self.nodes,
                );
            }
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
            &self.computation_order,
            &self.cyclic_derivatives,
            &self.nodes,
            &self.weights,
            &self.coefficients,
            &self.h_powers,
            &self.factorial,
        )
    }
}

pub mod list {
    use super::NDMORK;

    pub fn ENDMORK1_weight_function(_N: u32) -> Vec<Vec<f64>> {
        vec![vec![0.], vec![1.]]
    }

    pub fn ENDMORK1_nodes() -> Vec<f64> {
        vec![0., 1.]
    }

    pub fn ENDMORK1_weight_graph() -> Vec<Vec<bool>> {
        vec![vec![false, false], vec![true, false]]
    }

    pub fn ENDMORK1() -> NDMORK {
        NDMORK::new(
            Box::new(ENDMORK1_weight_function),
            ENDMORK1_nodes(),
            ENDMORK1_weight_graph(),
        )
    }

    pub fn INDMORK1_weight_function(_N: u32) -> Vec<Vec<f64>> {
        vec![vec![1.], vec![1.]]
    }

    pub fn INDMORK1_nodes() -> Vec<f64> {
        vec![1., 1.]
    }

    pub fn INDMORK1_weight_graph() -> Vec<Vec<bool>> {
        vec![vec![true, false], vec![true, false]]
    }

    pub fn INDMORK1() -> NDMORK {
        NDMORK::new(
            Box::new(INDMORK1_weight_function),
            INDMORK1_nodes(),
            INDMORK1_weight_graph(),
        )
    }

    pub fn ENDMORK2_weight_function(N: u32) -> Vec<Vec<f64>> {
        vec![
            vec![0., 0.],
            vec![2_f64.powi(-(N as i32)), 0.],
            vec![1. - 2. / (1. + N as f64), 2. / (1. + N as f64)],
        ]
    }

    pub fn ENDMORK2_nodes() -> Vec<f64> {
        vec![0., 0.5, 1.]
    }

    pub fn ENDMORK2_weight_graph() -> Vec<Vec<bool>> {
        vec![
            vec![false, false, false],
            vec![true, false, false],
            vec![true, true, false],
        ]
    }

    pub fn ENDMORK2() -> NDMORK {
        NDMORK::new(
            Box::new(ENDMORK2_weight_function),
            ENDMORK2_nodes(),
            ENDMORK2_weight_graph(),
        )
    }

    pub fn ENDMORK2bis_weight_function(N: u32) -> Vec<Vec<f64>> {
        vec![
            vec![0., 0.],
            vec![(2_f64 / 3_f64).powi(N as i32), 0.],
            vec![
                (2. * N as f64 - 1.) / (2. * (1. + N as f64)),
                3. / (2. * (1. + N as f64)),
            ],
        ]
    }

    pub fn ENDMORK2bis_nodes() -> Vec<f64> {
        vec![0., 2. / 3., 1.]
    }

    pub fn ENDMORK2bis_weight_graph() -> Vec<Vec<bool>> {
        vec![
            vec![false, false, false],
            vec![true, false, false],
            vec![true, true, false],
        ]
    }

    pub fn ENDMORK2bis() -> NDMORK {
        NDMORK::new(
            Box::new(ENDMORK2_weight_function),
            ENDMORK2_nodes(),
            ENDMORK2_weight_graph(),
        )
    }

    pub fn INDMORK2_weight_function(N: u32) -> Vec<Vec<f64>> {
        vec![vec![2_f64.powi(-(N as i32))], vec![1.]]
    }

    pub fn INDMORK2_nodes() -> Vec<f64> {
        vec![0.5, 1.]
    }

    pub fn INDMORK2_weight_graph() -> Vec<Vec<bool>> {
        vec![vec![true, false], vec![true, false]]
    }

    pub fn INDMORK2() -> NDMORK {
        NDMORK::new(
            Box::new(INDMORK2_weight_function),
            INDMORK2_nodes(),
            INDMORK2_weight_graph(),
        )
    }

    pub fn ENDMORK3_weight_function(N: u32) -> Vec<Vec<f64>> {
        let N = N as i32;
        let Nf = N as f64;
        vec![
            vec![0., 0., 0.],
            vec![3_f64.powi(-N), 0., 0.],
            vec![
                (2. / 3_f64).powi(N) * (Nf - 1.) / (1. + Nf),
                (2. / 3_f64).powi(N) * 2. / (1. + Nf),
                0.,
            ],
            vec![
                1. - 9. * Nf / (2. * (1. + Nf) * (2. + Nf)),
                6. * (Nf - 1.) / ((1. + Nf) * (2. + Nf)),
                3. * (4. - Nf) / (2. * (1. + Nf) * (2. + Nf)),
            ],
        ]
    }

    pub fn ENDMORK3_nodes() -> Vec<f64> {
        vec![0., 1. / 3., 2. / 3., 1.]
    }

    pub fn ENDMORK3_weight_graph() -> Vec<Vec<bool>> {
        vec![
            vec![false, false, false, false],
            vec![true, false, false, false],
            vec![true, true, false, false],
            vec![true, true, true, false],
        ]
    }

    pub fn ENDMORK3() -> NDMORK {
        NDMORK::new(
            Box::new(ENDMORK3_weight_function),
            ENDMORK3_nodes(),
            ENDMORK3_weight_graph(),
        )
    }

    pub fn INDMORK3_weight_function(_N: u32) -> Vec<Vec<f64>> {
        vec![vec![0., 0.], vec![0., 1.], vec![0.5, 0.5]]
    }

    pub fn INDMORK3_nodes() -> Vec<f64> {
        vec![0., 1., 1.]
    }

    pub fn INDMORK3_weight_graph() -> Vec<Vec<bool>> {
        vec![
            vec![false, false, false],
            vec![false, true, false],
            vec![true, true, false],
        ]
    }

    pub fn INDMORK3() -> NDMORK {
        NDMORK::new(
            Box::new(INDMORK3_weight_function),
            INDMORK3_nodes(),
            INDMORK3_weight_graph(),
        )
    }

    pub fn INDMORK3_1_weight_function(_N: u32) -> Vec<Vec<f64>> {
        vec![vec![0., 0.], vec![0.5, 0.5], vec![0.5, 0.5]]
    }

    pub fn INDMORK3_1_nodes() -> Vec<f64> {
        vec![0., 1., 1.]
    }

    pub fn INDMORK3_1_weight_graph() -> Vec<Vec<bool>> {
        vec![
            vec![false, false, false],
            vec![true, true, false],
            vec![true, true, false],
        ]
    }

    pub fn INDMORK3_1() -> NDMORK {
        NDMORK::new(
            Box::new(INDMORK3_1_weight_function),
            INDMORK3_1_nodes(),
            INDMORK3_1_weight_graph(),
        )
    }

    pub fn ENDMORK4_1_weight_function(N: u32) -> Vec<Vec<f64>> {
        let N = N as i32;
        let Nf = N as f64;
        vec![
            vec![0., 0., 0., 0.],
            vec![2_f64.powi(-N), 0., 0., 0.],
            vec![
                2_f64.powi(-N) * (Nf - 1.) / (1. + Nf),
                2_f64.powi(1 - N) / (1. + Nf),
                0.,
                0.,
            ],
            vec![(Nf - 1.) / (1. + Nf), (1. - Nf) / (1. + Nf), 1., 0.],
            vec![
                Nf.powi(2) / ((1. + Nf) * (2. + Nf)),
                2. * Nf / ((1. + Nf) * (2. + Nf)),
                2. * Nf / ((1. + Nf) * (2. + Nf)),
                (2. - Nf) / ((1. + Nf) * (2. + Nf)),
            ],
        ]
    }

    pub fn ENDMORK4_1_nodes() -> Vec<f64> {
        vec![0., 0.5, 0.5, 1., 1.]
    }

    pub fn ENDMORK4_1_weight_graph() -> Vec<Vec<bool>> {
        vec![
            vec![false, false, false, false, false],
            vec![true, false, false, false, false],
            vec![true, true, false, false, false],
            vec![true, true, true, false, false],
            vec![true, true, true, true, false],
        ]
    }

    pub fn ENDMORK4_1() -> NDMORK {
        NDMORK::new(
            Box::new(ENDMORK4_1_weight_function),
            ENDMORK4_1_nodes(),
            ENDMORK4_1_weight_graph(),
        )
    }

    pub fn ENDMORK4_2_weight_function(N: u32) -> Vec<Vec<f64>> {
        let N = N as i32;
        let Nf = N as f64;
        vec![
            vec![0., 0., 0., 0.],
            vec![2_f64.powi(-N), 0., 0., 0.],
            vec![
                2_f64.powi(-N) * Nf / (1. + Nf),
                2_f64.powi(-N) / (1. + Nf),
                0.,
                0.,
            ],
            vec![
                (Nf - 1.) / (1. + Nf),
                2. * (Nf - 2.) / (1. + Nf),
                2. * (3. - Nf) / (1. + Nf),
                0.,
            ],
            vec![
                Nf.powi(2) / ((1. + Nf) * (2. + Nf)),
                0.,
                4. * Nf / ((1. + Nf) * (2. + Nf)),
                (2. - Nf) / ((1. + Nf) * (2. + Nf)),
            ],
        ]
    }

    pub fn ENDMORK4_2_nodes() -> Vec<f64> {
        vec![0., 0.5, 0.5, 1., 1.]
    }

    pub fn ENDMORK4_2_weight_graph() -> Vec<Vec<bool>> {
        vec![
            vec![false, false, false, false, false],
            vec![true, false, false, false, false],
            vec![true, true, false, false, false],
            vec![true, true, true, false, false],
            vec![true, true, true, true, false],
        ]
    }

    pub fn ENDMORK4_2() -> NDMORK {
        NDMORK::new(
            Box::new(ENDMORK4_2_weight_function),
            ENDMORK4_2_nodes(),
            ENDMORK4_2_weight_graph(),
        )
    }

    pub fn INDMORK4_weight_function(N: u32) -> Vec<Vec<f64>> {
        let sqrt3 = 3_f64.sqrt();
        let no1 = 0.5 - 3_f64.sqrt() / 6.;
        let no2 = 0.5 + 3_f64.sqrt() / 6.;
        let Nf = N as f64;
        vec![
            vec![
                no1.powi(N as i32) / (1. + Nf) * (1. + Nf / 2. * (1. + sqrt3)),
                -sqrt3 * Nf / (1. + Nf) * no1.powi(N as i32 + 1),
            ],
            vec![
                sqrt3 * Nf / (1. + Nf) * no2.powi(N as i32 + 1),
                no2.powi(N as i32) / (1. + Nf) * (1. + Nf / 2. * (1. - sqrt3)),
            ],
            vec![
                0.5 + sqrt3 * (Nf - 1.) / (2. * (1. + Nf)),
                0.5 - sqrt3 * (Nf - 1.) / (2. * (1. + Nf)),
            ],
        ]
    }

    pub fn INDMORK4_nodes() -> Vec<f64> {
        vec![0.5 - 3_f64.sqrt() / 6., 0.5 + 3_f64.sqrt() / 6., 1.]
    }

    pub fn INDMORK4_weight_graph() -> Vec<Vec<bool>> {
        vec![
            vec![true, true, false],
            vec![true, true, false],
            vec![true, true, false],
        ]
    }

    pub fn INDMORK4() -> NDMORK {
        NDMORK::new(
            Box::new(INDMORK4_weight_function),
            INDMORK4_nodes(),
            INDMORK4_weight_graph(),
        )
    }
}