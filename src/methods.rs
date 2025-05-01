#![allow(non_snake_case)]

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
pub enum Task {
	Implicit(Vec<usize>),
	Explicit(usize),
}

pub struct GMORK {
	pub nodes: Vec<f64>,
	pub main_weights: Vec<Vec<Vec<f64>>>,
	main_weights_f: fn(u32) -> Vec<Vec<f64>>,
	pub secondary_weights: Vec<Vec<Vec<f64>>>,
	secondary_weights_f: fn(u32) -> Vec<Vec<f64>>,
	pub queue: Vec<Task>,
	pub s: usize,
	pub length: usize,
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
		main_weights_function: fn(u32) -> Vec<Vec<f64>>,
		secondary_weights_function: fn(u32) -> Vec<Vec<f64>>,
		nodes: Vec<f64>,
		maximum_weight_graph: Vec<Vec<bool>>,
	) -> Self {
		let SCC = scc(&maximum_weight_graph);
		let queue: Vec<Task> = topological_sort(&contraction(&maximum_weight_graph, &SCC))
			.into_iter()
			.map(|i| {
				if SCC[i].len() > 1 || maximum_weight_graph[SCC[i][0]][SCC[i][0]] {
					Task::Implicit(SCC[i].clone())
				} else {
					Task::Explicit(SCC[i][0])
				}
			})
			.collect();
		let main_weights = vec![main_weights_function(1)];
		let s = main_weights[0].len() - 1;
		let secondary_weights = vec![secondary_weights_function(1)];
		let mut cyclic_derivatives: Vec<Vec<bool>> = vec![vec![false; s]];
		for task in queue.iter() {
			if let Task::Implicit(J) = task {
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
			length: 1,
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

	pub fn set_minimum_length(&mut self, n: usize) {
		if self.length >= n {
			return;
		}
		self.factorial.extend(vec![0.; n - self.length]);
		self.main_weights
			.extend((self.length..n).map(|N| (self.main_weights_f)(N as u32 + 1)));
		self.secondary_weights
			.extend((self.length..n).map(|N| (self.secondary_weights_f)(N as u32 + 1)));
		self.h_powers
			.extend((self.length..n).map(|N| self.h.powi(N as i32 + 1)));
		self.cyclic_derivatives
			.extend((self.length..n).map(|_| vec![false; self.s]));
		for N in self.length..n {
			self.factorial[N + 1] = self.factorial[N] * (N as f64 + 1.);
			for task in self.queue.iter() {
				if let Task::Implicit(J) = task {
					for &j in J {
						for &j1 in J {
							if self.main_weights[N][j][j1] != 0. {
								self.cyclic_derivatives[N][j] = true;
								break;
							}
						}
					}
				}
			}
		}
		self.length = n;
	}

	pub fn compute_explicit(&self, J: &Vec<usize>,j: usize, y: &mut Vec<Vec<Vec<f64>>>, y0: &Vec<Vec<f64>>, F: &Vec<Vec<f64>>) {
		let mut sum;
		for k in 0..y0.len() {
			for N in 0..y0[k].len() {
				for N1 in 1..=N {
					y[j][k][N] += self.secondary_weights[N][N1][j]
						* self.h_powers[N1]
						* y0[k][N - N1]
				}
				sum = 0.;
				for &j1 in J {
					sum += self.main_weights[N][j][j1] * F[j1][k];
				}
				y[j][k][N] += self.h_powers[N + 1] / self.factorial[N + 1] * sum;
			}
		}
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
			for N in 1..=self.length {
				self.h_powers[N] = h.powi(N as i32)
			}
		}
		let mut F: Vec<Vec<f64>> = (0..self.s).map(|_| y0[0].clone()).collect();
		let mut y: Vec<Vec<Vec<f64>>> = (0..=self.s).map(|_| y0.clone()).collect();
		let mut sum;
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
		// verify the length of the method is enough
		for k in 0..y0.len() {
			if y0[k].len() > self.length {
				self.set_minimum_length(y0[k].len());
			}
		}
		for task in self.queue.iter() {
			match task {
				Task::Explicit(j) => {
					let j = *j;
					for k in 0..y0.len() {
						for N in 0..y0[k].len() {
							for N1 in 1..=N {
								y[j][k][N] += self.secondary_weights[N][N1][j]
									* self.h_powers[N1]
									* y0[k][N - N1]
							}
							sum = 0.;
							for j1 in 0..self.s {
								sum += self.main_weights[N][j][j1] * F[j1][k];
							}
							y[j][k][N] += self.h_powers[N + 1] / self.factorial[N + 1] * sum;
						}
					}
					if j != self.s {
						F[j] = f(t + self.nodes[j] * h, &y[j]);
					}
				}
				Task::Implicit(J) => {
					let comp_J: &Vec<usize> = &(0..self.s).filter(|j| !J.contains(j)).collect();
					// calculate constant terms
					for &j in J {
						for k in 0..y0.len() {
							for N in 0..y0[k].len() {
								for N1 in 1..=N {
									y[j][k][N] += self.secondary_weights[N][N1][j]
										* self.h_powers[N1]
										* y0[k][N - N1]
								}
								sum = 0.;
								for &j1 in comp_J {
									sum += self.main_weights[N][j][j1] * F[j1][k];
								}
								y[j][k][N] += self.h_powers[N + 1] / self.factorial[N + 1] * sum;
							}
						}
					}
					let constant = y.clone();
					// Picard iterations
					let mut iter_count = 0;
					let mut d = threshold + 1.;
					let mut f_cache: Vec<f64>;
					while iter_count < self.min_iter
						|| (d > threshold && iter_count < self.max_iter)
					{
						iter_count += 1;
						// update evaluations and calculate difference
						d = 0.;
						for &j in J {
							f_cache = f(t + self.nodes[j] * h, &y[j]);
							// calculate difference for threshold
							for k in 0..f_cache.len() {
								d = d.max((f_cache[k] - F[j][k]).abs());
							}
							F[j] = f_cache;
						}
						// constant part
						y = constant.clone();
						// add evaluations and calculate norm difference
						for &j in J {
							for k in 0..y0.len() {
								for N in (0..y0[k].len()).filter(|&N| self.cyclic_derivatives[N][j]) {
									sum = 0.;
									for &j1 in J {
										sum += self.main_weights[N][j][j1] * F[j1][k];
									}
									y[j][k][N] +=
										self.h_powers[N + 1] / self.factorial[N + 1] * sum;
								}
							}
						}
					}
				}
			}
		}
		return y[self.s].clone();
	}
}
pub struct NDMORK {
	pub s: usize,
	pub current_length: usize,
	pub nodes: Vec<f64>,
	pub weights: Vec<Vec<Vec<f64>>>,
	weights_function: fn(u32) -> Vec<Vec<f64>>,
	pub factorial: Vec<f64>,
	pub coefficients: Vec<Vec<f64>>,
	pub h: f64,
	pub h_powers: Vec<f64>,
	pub queue: Vec<Task>,
	pub cycle_derivative: Vec<Vec<bool>>, // [N-1][j]
	pub error_fraction: f64,
	pub min_iter: u32,
	pub max_iter: u32,
}

impl NDMORK {
	pub fn new(
		weights_function: fn(u32) -> Vec<Vec<f64>>,
		nodes: Vec<f64>,
		maximum_weight_graph: Vec<Vec<bool>>,
	) -> Self {
		let SCC = scc(&maximum_weight_graph);
		let queue: Vec<Task> = topological_sort(&contraction(&maximum_weight_graph, &SCC))
			.into_iter()
			.map(|i| {
				if SCC[i].len() > 1 || maximum_weight_graph[SCC[i][0]][SCC[i][0]] {
					Task::Implicit(SCC[i].clone())
				} else {
					Task::Explicit(SCC[i][0])
				}
			})
			.collect();
		let s = nodes.len() - 1;
		let weights = vec![weights_function(1)];
		let mut cycle_derivative: Vec<Vec<bool>> = vec![vec![false; s]];
		for task in queue.iter() {
			if let Task::Implicit(J) = task {
				for &j in J {
					for &j1 in J {
						if weights[0][j][j1] != 0. {
							cycle_derivative[0][j] = true;
							break;
						}
					}
				}
			}
		}
		NDMORK {
			s,
			current_length: 1,
			nodes,
			weights,
			weights_function,
			factorial: vec![1., 1.],
			coefficients: vec![vec![1.; s + 1]],
			h: 0.,
			h_powers: vec![1., 0.],
			queue,
			cycle_derivative,
			error_fraction: ERROR_FRACTION,
			min_iter: MIN_ITER,
			max_iter: MAX_ITER,
		}
	}

	pub fn set_minimum_length(&mut self, n: usize) {
		if self.current_length >= n {
			return;
		}
		self.factorial.extend(vec![0.; n - self.current_length]);
		self.weights
			.extend((self.current_length..n).map(|N| (self.weights_function)(N as u32 + 1)));
		self.h_powers
			.extend((self.current_length..n).map(|N| self.h.powi(N as i32 + 1)));
		self.cycle_derivative
			.extend((self.current_length..n).map(|_| vec![false; self.s]));
		for N in self.current_length..n {
			self.factorial[N + 1] = self.factorial[N] * (N as f64 + 1.);
			for task in self.queue.iter() {
				if let Task::Implicit(J) = task {
					for &j in J {
						for &j1 in J {
							if self.weights[N][j][j1] != 0. {
								self.cycle_derivative[N][j] = true;
								break;
							}
						}
					}
				}
			}
		}
		self.coefficients.extend((self.current_length..n).map(|N| {
			(0..self.s + 1)
				.map(|j| self.nodes[j].powi(N as i32) / self.factorial[N])
				.collect()
		}));
		self.current_length = n;
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
			for N in 1..=self.current_length {
				self.h_powers[N] = h.powi(N as i32)
			}
		}
		let mut F: Vec<Vec<f64>> = (0..self.s).map(|_| y0[0].clone()).collect();
		let mut y: Vec<Vec<Vec<f64>>> = (0..=self.s).map(|_| y0.clone()).collect();
		let mut sum;
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
		// verify the length of the method is enough
		for k in 0..y0.len() {
			if y0[k].len() > self.current_length {
				self.set_minimum_length(y0[k].len());
			}
		}
		for task in self.queue.iter() {
			match task {
				Task::Explicit(j) => {
					let j = *j;
					for k in 0..y0.len() {
						for N in 0..y0[k].len() {
							if self.nodes[j] != 0. {
								for N1 in 1..=N {
									y[j][k][N] +=
										self.coefficients[N1][j] * self.h_powers[N1] * y0[k][N - N1]
								}
							}
							sum = 0.;
							for j1 in 0..self.s {
								sum += self.weights[N][j][j1] * F[j1][k];
							}
							y[j][k][N] += self.h_powers[N + 1] / self.factorial[N + 1] * sum;
						}
					}
					if j != self.s {
						F[j] = f(t + self.nodes[j] * h, &y[j]);
					}
				}
				Task::Implicit(J) => {
					let comp_J: &Vec<usize> = &(0..self.s).filter(|j| !J.contains(j)).collect();
					// calculate constant terms
					for &j in J {
						for k in 0..y0.len() {
							for N in 0..y0[k].len() {
								if self.nodes[j] != 0. {
									for N1 in 1..=N {
										y[j][k][N] += self.coefficients[N1][j]
											* self.h_powers[N1]
											* y0[k][N - N1]
									}
								}
								sum = 0.;
								for &j1 in comp_J {
									sum += self.weights[N][j][j1] * F[j1][k];
								}
								y[j][k][N] += self.h_powers[N + 1] / self.factorial[N + 1] * sum;
							}
						}
					}
					let constant = y.clone();
					// Picard iterations
					let mut iter_count = 0;
					let mut d = threshold + 1.;
					let mut f_cache: Vec<f64>;
					while iter_count < self.min_iter
						|| (d > threshold && iter_count < self.max_iter)
					{
						iter_count += 1;
						// update evaluations and calculate difference
						d = 0.;
						for &j in J {
							f_cache = f(t + self.nodes[j] * h, &y[j]);
							// calculate difference for threshold
							for k in 0..f_cache.len() {
								d = d.max((f_cache[k] - F[j][k]).abs());
							}
							F[j] = f_cache;
						}
						// constant part
						y = constant.clone();
						// add evaluations and calculate norm difference
						for &j in J {
							for k in 0..y0.len() {
								for N in (0..y0[k].len()).filter(|&N| self.cycle_derivative[N][j]) {
									sum = 0.;
									for &j1 in J {
										sum += self.weights[N][j][j1] * F[j1][k];
									}
									y[j][k][N] +=
										self.h_powers[N + 1] / self.factorial[N + 1] * sum;
								}
							}
						}
					}
				}
			}
		}
		return y[self.s].clone();
	}
}

pub fn ENDMORK1() -> NDMORK {
	let weight_function: fn(_) -> _ = |_N: u32| vec![vec![0.], vec![1.]];
	let nodes = vec![0., 1.];
	let weight_graph = vec![vec![false, false], vec![true, false]];
	NDMORK::new(weight_function, nodes, weight_graph)
}

pub fn INDMORK1() -> NDMORK {
	let weight_function: fn(_) -> _ = |_N: u32| vec![vec![1.], vec![1.]];
	let nodes = vec![1., 1.];
	let weight_graph = vec![vec![true, false], vec![true, false]];
	NDMORK::new(weight_function, nodes, weight_graph)
}

pub fn ENDMORK2() -> NDMORK {
	let weight_function: fn(_) -> _ = |N: u32| {
		vec![
			vec![0., 0.],
			vec![2_f64.powi(-(N as i32)), 0.],
			vec![1. - 2. / (1. + N as f64), 2. / (1. + N as f64)],
		]
	};
	let nodes = vec![0., 0.5, 1.];
	let weight_graph = vec![
		vec![false, false, false],
		vec![true, false, false],
		vec![true, true, false],
	];
	NDMORK::new(weight_function, nodes, weight_graph)
}

pub fn INDMORK2() -> NDMORK {
	let weight_function: fn(_) -> _ = |N: u32| vec![vec![2_f64.powi(-(N as i32))], vec![1.]];
	let nodes = vec![0.5, 1.];
	let weight_graph = vec![vec![true, false], vec![true, false]];
	NDMORK::new(weight_function, nodes, weight_graph)
}

pub fn INDMORK3() -> NDMORK {
	let weight_function: fn(_) -> _ = |_N: u32| vec![vec![0., 0.], vec![0., 1.], vec![0.5, 0.5]];
	let nodes = vec![0., 1., 1.];
	let weight_graph = vec![
		vec![false, false, false],
		vec![false, true, false],
		vec![true, true, false],
	];
	NDMORK::new(weight_function, nodes, weight_graph)
}

pub fn INDMORK3_1() -> NDMORK {
	let weight_function: fn(_) -> _ = |_N: u32| vec![vec![0., 0.], vec![0.5, 0.5], vec![0.5, 0.5]];
	let nodes = vec![0., 1., 1.];
	let weight_graph = vec![
		vec![false, false, false],
		vec![true, true, false],
		vec![true, true, false],
	];
	NDMORK::new(weight_function, nodes, weight_graph)
}

pub fn INDMORK3_2() -> NDMORK {
	let weight_function: fn(_) -> _ = |N: u32| {
		if N == 1 {
			return vec![vec![0., 0.], vec![0.5, 0.5], vec![0.5, 0.5]];
		}
		vec![vec![0., 0.], vec![0., 0.], vec![0.5, 0.5]]
	};
	let nodes = vec![0., 1., 1.];
	let weight_graph = vec![
		vec![false, false, false],
		vec![true, true, false],
		vec![true, true, false],
	];
	NDMORK::new(weight_function, nodes, weight_graph)
}

pub fn ENDMORK3() -> NDMORK {
	let weight_function: fn(_) -> _ = |N: u32| {
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
	};
	let nodes = vec![0., 1. / 3., 2. / 3., 1.];
	let weight_graph = vec![
		vec![false, false, false, false],
		vec![true, false, false, false],
		vec![true, true, false, false],
		vec![true, true, true, false],
	];
	NDMORK::new(weight_function, nodes, weight_graph)
}

pub fn ENDMORK4_1() -> NDMORK {
	let weight_function: fn(_) -> _ = |N: u32| {
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
	};
	let nodes = vec![0., 0.5, 0.5, 1., 1.];
	let weight_graph = vec![
		vec![false, false, false, false, false],
		vec![true, false, false, false, false],
		vec![true, true, false, false, false],
		vec![true, true, true, false, false],
		vec![true, true, true, true, false],
	];
	NDMORK::new(weight_function, nodes, weight_graph)
}

pub fn ENDMORK4_2() -> NDMORK {
	let weight_function: fn(_) -> _ = |N: u32| {
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
	};
	let nodes = vec![0., 0.5, 0.5, 1., 1.];
	let weight_graph = vec![
		vec![false, false, false, false, false],
		vec![true, false, false, false, false],
		vec![true, true, false, false, false],
		vec![true, true, true, false, false],
		vec![true, true, true, true, false],
	];
	NDMORK::new(weight_function, nodes, weight_graph)
}

pub fn INDMORK4() -> NDMORK {
	let weight_function: fn(_) -> _ = |N: u32| {
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
	};
	let nodes = vec![0.5 - 3_f64.sqrt() / 6., 0.5 + 3_f64.sqrt() / 6., 1.];
	let weight_graph = vec![
		vec![true, true, false],
		vec![true, true, false],
		vec![true, true, false],
	];
	NDMORK::new(weight_function, nodes, weight_graph)
}

pub struct RK {
	pub s: usize,
	pub nodes: Vec<f64>,
	pub weights: Vec<Vec<f64>>,
	pub queue: Vec<Task>,
	pub error_fraction: f64,
	pub min_tries: u32,
	pub max_tries: u32,
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
		let SCC = scc(&weight_graph);
		let queue: Vec<Task> = topological_sort(&contraction(&weight_graph, &SCC))
			.into_iter()
			.map(|i| {
				if SCC[i].len() > 1 || weight_graph[SCC[i][0]][SCC[i][0]] {
					Task::Implicit(SCC[i].clone())
				} else {
					Task::Explicit(SCC[i][0])
				}
			})
			.collect();
		let s = nodes.len() - 1;
		RK {
			s,
			nodes,
			weights,
			queue,
			error_fraction: ERROR_FRACTION,
			min_tries: MIN_ITER,
			max_tries: MAX_ITER,
		}
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
		let mut F: Vec<Vec<f64>> = (0..self.s).map(|_| y0[0].clone()).collect();
		let mut y: Vec<Vec<Vec<f64>>> = (0..=self.s).map(|_| y0.clone()).collect();
		let mut sum;
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
		for task in self.queue.iter() {
			match task {
				Task::Explicit(j) => {
					let j = *j;
					for k in 0..y0.len() {
						sum = 0.;
						for j1 in 0..self.s {
							sum += self.weights[j][j1] * F[j1][k];
						}
						y[j][k][0] += h * sum;
						for N in 1..y0[k].len() {
							sum = 0.;
							for j1 in 0..self.s {
								sum += self.weights[j][j1] * y[j1][k][N - 1];
							}
							y[j][k][N] += h * sum;
						}
					}
					if j != self.s {
						F[j] = f(t + self.nodes[j] * h, &y[j]);
					}
				}
				Task::Implicit(J) => {
					let constant = y.clone();
					// Picard iterations
					let mut iter_count = 0;
					let mut d = threshold + 1.;
					let mut f_cache;
					while iter_count < self.min_tries
						|| (d > threshold && iter_count < self.max_tries)
					{
						iter_count += 1;
						// update evaluations and calculate difference
						d = 0.;
						for &j in J {
							f_cache = f(t + self.nodes[j] * h, &y[j]);
							// calculate difference for threshold
							for k in 0..f_cache.len() {
								d = d.max((f_cache[k] - F[j][k]).abs());
							}
							F[j] = f_cache;
						}
						// constant part
						y = constant.clone();
						// add evaluations and calculate norm difference
						for &j in J {
							for k in 0..y0.len() {
								for N in (1..y0[k].len()).rev() {
									sum = 0.;
									for j1 in 0..self.s {
										sum += self.weights[j][j1] * y[j1][k][N - 1];
									}
									y[j][k][N] += h * sum;
								}
								sum = 0.;
								for j1 in 0..self.s {
									sum += self.weights[j][j1] * F[j1][k];
								}
								y[j][k][0] += h * sum;
							}
						}
					}
				}
			}
		}
		return y[self.s].clone();
	}
}

pub fn ERK1() -> RK {
	let weights = vec![vec![0.], vec![1.]];
	let nodes = vec![0., 1.];
	RK::new(weights, nodes)
}

pub fn IRK1() -> RK {
	let weights = vec![vec![1.], vec![1.]];
	let nodes = vec![1., 1.];
	RK::new(weights, nodes)
}

pub fn ERK2() -> RK {
	let weights = vec![vec![0., 0.], vec![0.5, 0.], vec![0., 1.]];
	let nodes = vec![0., 0.5, 1.];
	RK::new(weights, nodes)
}

pub fn IRK2() -> RK {
	let weights = vec![vec![0.5], vec![1.]];
	let nodes = vec![0.5, 1.];
	RK::new(weights, nodes)
}

pub fn IRK3() -> RK {
	let weights = vec![vec![0., 0.], vec![0., 1.], vec![0.5, 0.5]];
	let nodes = vec![0., 1., 1.];
	RK::new(weights, nodes)
}

pub fn IRK3_1() -> RK {
	let weights = vec![vec![0., 0.], vec![0.5, 0.5], vec![0.5, 0.5]];
	let nodes = vec![0., 1., 1.];
	RK::new(weights, nodes)
}

pub fn IRK3_2() -> RK {
	let weights = vec![vec![0., 0.], vec![0.5, 0.5], vec![0.5, 0.5]];
	let nodes = vec![0., 1., 1.];
	RK::new(weights, nodes)
}

pub fn ERK3() -> RK {
	let weights = vec![
		vec![0., 0., 0.],
		vec![1. / 3., 0., 0.],
		vec![0., 2. / 3., 0.],
		vec![1. / 4., 0., 3. / 4.],
	];
	let nodes = vec![0., 1. / 3., 2. / 3., 1.];
	RK::new(weights, nodes)
}

pub fn ERK4_1() -> RK {
	let weights = {
		vec![
			vec![0., 0., 0., 0.],
			vec![0.5, 0., 0., 0.],
			vec![0., 0.5, 0., 0.],
			vec![0., 0., 1., 0.],
			vec![1. / 6., 1. / 3., 1. / 3., 1. / 6.],
		]
	};
	let nodes = vec![0., 0.5, 0.5, 1., 1.];
	RK::new(weights, nodes)
}

pub fn ERK4_2() -> RK {
	let weights = vec![
		vec![0., 0., 0., 0.],
		vec![0.5, 0., 0., 0.],
		vec![0.25, 0.25, 0., 0.],
		vec![0., -1., 2., 0.],
		vec![1. / 6., 0., 2. / 3., 1. / 6.],
	];
	let nodes = vec![0., 0.5, 0.5, 1., 1.];
	RK::new(weights, nodes)
}

pub fn IRK4() -> RK {
	let sqrt3 = 3_f64.sqrt();
	let weights = vec![
		vec![1. / 4., 1. / 4. - sqrt3 / 6.],
		vec![1. / 4. + sqrt3 / 6., 1. / 4.],
		vec![0.5, 0.5],
	];
	let nodes = vec![0.5 - sqrt3 / 6., 0.5 + sqrt3 / 6., 1.];
	RK::new(weights, nodes)
}
