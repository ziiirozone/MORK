use crate::ivp::IVPType;
use crate::methods::Solver;
use crate::plot::{log_plot, plot};
use eframe::egui::{self, Color32, RichText, Spinner, Widget};
use plotters::style::Color;
use plotters::style::{Palette, Palette99};

use serde_cbor::to_vec;
use std::fs::File;
use std::io::Write;

//const RENDER_PATH_WINDOWS: &str = "D:/programs/rust/MORK/src/render";
const RENDER_PATH_LINUX: &str = "/run/media/ziii/SSD_loris/programs/rust/MORK/src/render";
const RENDER_PATH: &str = RENDER_PATH_LINUX;

#[derive(Debug, PartialEq, Clone)]
enum Goal {
	Order,
	Plot,
}
pub struct ExperimentPicker {
	ivps: Vec<(String, IVPType, Vec<String>)>,
	solvers: Vec<(String, Box<dyn Solver>)>,
	extractors: Vec<(String, Box<dyn Fn(&Vec<Vec<f64>>) -> Vec<f64>>)>,
	distances: Vec<(
		String,
		Box<dyn Fn(&Vec<Vec<f64>>, &Vec<Vec<f64>>) -> Vec<f64>>,
	)>,
	ivp_index: usize,
	solvers_selected: Vec<bool>,
	extractor_index: usize,
	distance_index: usize,
	goal: Goal,
	plot_solution: bool,
	step_size: f64,
	left_h_interval: f64,
	right_h_interval: f64,
	order_samples: u32,
	delta_t: f64,
	t0: f64,
	order_plot_error: bool,
}

impl ExperimentPicker {
	pub fn new(
		experiments: Vec<(impl ToString, IVPType, Vec<impl ToString>)>,
		solvers: Vec<(impl ToString, Box<dyn Solver>)>,
		extractors: Vec<(impl ToString, Box<dyn Fn(&Vec<Vec<f64>>) -> Vec<f64>>)>,
		distances: Vec<(
			impl ToString,
			Box<dyn Fn(&Vec<Vec<f64>>, &Vec<Vec<f64>>) -> Vec<f64>>,
		)>,
	) -> Self {
		let len_solv = solvers.len();
		ExperimentPicker {
			ivps: experiments
				.into_iter()
				.map(|e| {
					(
						e.0.to_string(),
						e.1,
						e.2.into_iter().map(|s| s.to_string()).collect(),
					)
				})
				.collect(),
			solvers: solvers
				.into_iter()
				.map(|e| (e.0.to_string(), e.1))
				.collect(),
			extractors: extractors
				.into_iter()
				.map(|e| (e.0.to_string(), e.1))
				.collect(),
			distances: distances
				.into_iter()
				.map(|e| (e.0.to_string(), e.1))
				.collect(),
			ivp_index: 0,
			solvers_selected: vec![false; len_solv],
			extractor_index: 0,
			distance_index: 0,
			goal: Goal::Plot,
			plot_solution: true,
			step_size: 0.1,
			left_h_interval: 0.001,
			right_h_interval: 1.,
			order_samples: 20,
			delta_t: 10.,
			t0: 0.,
			order_plot_error: false,
		}
	}
}

impl eframe::App for ExperimentPicker {
	fn update(&mut self, ctx: &egui::Context, _frame: &mut eframe::Frame) {
		egui::CentralPanel::default().show(ctx, |ui| {
			// Select initial value problem
			egui::ComboBox::from_label("Experiments")
				.selected_text(&self.ivps[self.ivp_index].0)
				.show_ui(ui, |ui| {
					for i in 0..self.ivps.len() {
						ui.selectable_value(&mut self.ivp_index, i, &self.ivps[i].0);
					}
				});
			// list of solvers
			for i in 0..self.solvers.len() {
				ui.checkbox(&mut self.solvers_selected[i], &self.solvers[i].0);
			}
			// select goal
			ui.horizontal(|ui| {
				ui.selectable_value(&mut self.goal, Goal::Plot, "Plot");
				ui.selectable_value(&mut self.goal, Goal::Order, "Order");
			});
			// Options according to goal
			match &self.goal {
				Goal::Order => {
					// select distance function
					egui::ComboBox::from_label("Distances")
						.selected_text(&self.distances[self.distance_index].0)
						.show_ui(ui, |ui| {
							for i in 0..self.distances.len() {
								ui.selectable_value(
									&mut self.distance_index,
									i,
									&self.distances[i].0,
								);
							}
						});
					ui.label("Samples");
					ui.add(egui::DragValue::new(&mut self.order_samples).max_decimals(0));
					ui.label("Left interval");
					ui.add(egui::DragValue::new(&mut self.left_h_interval).max_decimals(99));
					ui.label("Right interval");
					ui.add(egui::DragValue::new(&mut self.right_h_interval).max_decimals(99));
				}
				Goal::Plot => {
					// select data extractor
					egui::ComboBox::from_label("Extractor")
						.selected_text(&self.extractors[self.extractor_index].0)
						.show_ui(ui, |ui| {
							for i in 0..self.extractors.len() {
								ui.selectable_value(
									&mut self.extractor_index,
									i,
									&self.extractors[i].0,
								);
							}
						});
					ui.checkbox(&mut self.plot_solution, "Plot solution if possible");
					ui.label("Step size");
					ui.add(egui::DragValue::new(&mut self.step_size).max_decimals(99));
					ui.label("Delta t");
					ui.add(egui::DragValue::new(&mut self.delta_t).max_decimals(99));
				}
			}
			// change paramaters of ivp
			match &mut self.ivps[self.ivp_index] {
				(_, IVPType::NotSolved(ivp), labels) => {
					for i in 0..ivp.parameters.len() {
						ui.label(&labels[i]);
						ui.add(egui::DragValue::new(&mut ivp.parameters[i]));
					}
				}
				(_, IVPType::Solved(ivp), labels) => {
					// If solved ivp select initial instant
					ui.label("t0");
					ui.add(egui::DragValue::new(&mut self.t0));
					for i in 0..ivp.parameters.len() {
						ui.label(&labels[i]);
						ui.add(egui::DragValue::new(&mut ivp.parameters[i]));
					}
				}
			}
			// runs the experiment
			if ui.button("Run").clicked() {
				let spinner = Spinner::default();
				spinner.ui(ui);
				let ivp_name = self.ivps[self.ivp_index].0.clone();
				match (&self.goal, &mut self.ivps[self.ivp_index].1) {
					(Goal::Order, IVPType::NotSolved(_)) => {
						self.order_plot_error = true;
					}
					(Goal::Order, IVPType::Solved(ivp)) => {
						self.order_plot_error = false;
						let plot_path = RENDER_PATH.to_string() + "/order.svg";
						let cbor_path = RENDER_PATH.to_string() + "/order.cbor";
						let title = format!("Order - {}", ivp_name);
						let mut names = Vec::new();
						let mut data: Vec<Vec<(f64, f64)>> = Vec::new();
						let ratio = self.right_h_interval / self.left_h_interval;
						let h_list: Vec<f64> = (0..self.order_samples)
							.map(|i| {
								self.left_h_interval
									* ratio.powf(i as f64 / (self.order_samples - 1) as f64)
							})
							.collect();
						let solution: Vec<Vec<Vec<f64>>> =
							h_list.iter().map(|t| ivp.solution(*t)).collect();
						let y0 = ivp.solution(self.t0);
						let f = |t: f64, y: &Vec<Vec<f64>>| ivp.f(t, y);
						let distance = &self.distances[self.distance_index].1;
						for solver_i in
							(0..self.solvers_selected.len()).filter(|i| self.solvers_selected[*i])
						{
							let temp: Vec<Vec<f64>> = h_list
								.iter()
								.enumerate()
								.map(|(i, h)| {
									let a =
										self.solvers[solver_i].1.approximate(self.t0, *h, &f, &y0);
									distance(&a, &solution[i])
								})
								.collect();
							for k in 0..temp[0].len() {
								data.push(
									(0..self.order_samples as usize)
										.map(|q| (h_list[q], temp[q][k]))
										.filter(|(x, y)| *x != 0. && *y != 0.)
										.collect(),
								);
								names.push(Some(self.solvers[solver_i].0.clone()));
							}
						}
						let cbor_data: Vec<(Vec<f64>, Vec<f64>)> = (0..data.len())
							.map(|q| {
								(
									(0..data[q].len()).map(|k| data[q][k].0).collect(),
									(0..data[q].len()).map(|k| data[q][k].1).collect(),
								)
							})
							.collect();
						let encoded = to_vec(&cbor_data).expect("failed serialization");
						let mut file = File::create(cbor_path).unwrap();
						file.write_all(&encoded).unwrap();
						log_plot(
							&data,
							&(0..data.len())
								.map(|i| Palette99::pick(i).stroke_width(7))
								.collect(),
							&names,
							plot_path.to_string(),
							title.to_string(),
							("h".to_string(), "e".to_string()),
						);
					}
					(Goal::Plot, IVPType::Solved(ivp)) => {
						self.order_plot_error = false;
						let plot_path = RENDER_PATH.to_string() + "/time.svg";
						let cbor_path = RENDER_PATH.to_string() + "/time.cbor";
						let title = format!("Plot - {}", ivp_name);
						let mut names = Vec::new();
						let mut data = Vec::new();
						let q = (self.delta_t / self.step_size) as usize;
						let t: Vec<f64> = (0..=q)
							.map(|q1| self.t0 + self.step_size * q1 as f64)
							.collect();
						let y0 = (ivp.solution)(self.t0, &ivp.parameters);
						let extractor = &self.extractors[self.extractor_index].1;
						if self.plot_solution {
							let y: Vec<Vec<f64>> = t
								.iter()
								.map(|t1| extractor(&(ivp.solution)(*t1, &ivp.parameters)))
								.collect();
							for k in 0..y[0].len() {
								let y1: Vec<(f64, f64)> =
									t.iter().enumerate().map(|(i, t1)| (*t1, y[i][k])).collect();
								data.push(y1);
								names.push(Some("Solution".to_string() + &k.to_string()));
							}
						}
						let f = |t: f64, y: &Vec<Vec<f64>>| ivp.f(t, y);
						for (_, (solver_name, solver)) in self
							.solvers
							.iter_mut()
							.enumerate()
							.filter(|x| self.solvers_selected[x.0])
						{
							let mut y: Vec<Vec<Vec<f64>>> = t.iter().map(|_| y0.clone()).collect();
							for k in 0..q {
								y[k + 1] = solver.approximate(t[k], self.step_size, &f, &y[k]);
							}
							let y1: Vec<Vec<f64>> = y.into_iter().map(|x| extractor(&x)).collect();
							for k in 0..y1[0].len() {
								let y2 = t
									.iter()
									.enumerate()
									.map(|(i, t1)| (*t1, y1[i][k]))
									.collect();
								data.push(y2);
								names.push(Some(solver_name.clone() + &k.to_string()));
							}
						}
						let cbor_data: Vec<(Vec<f64>, Vec<f64>)> = (0..data.len())
							.map(|q| {
								(
									(0..data[q].len()).map(|k| data[q][k].0).collect(),
									(0..data[q].len()).map(|k| data[q][k].1).collect(),
								)
							})
							.collect();
						let encoded = to_vec(&cbor_data).expect("failed serialization");
						let mut file = File::create(cbor_path).unwrap();
						file.write_all(&encoded).unwrap();
						plot(
							&data,
							&(0..data.len())
								.map(|i| Palette99::pick(i).stroke_width(7))
								.collect(),
							&names,
							plot_path,
							title,
							("t".to_string(), "y".to_string()),
						);
					}
					(Goal::Plot, IVPType::NotSolved(ivp)) => {
						self.order_plot_error = false;
						let plot_path = RENDER_PATH.to_string() + "/time.svg";
						let cbor_path = RENDER_PATH.to_string() + "/time.cbor";
						let title = format!("Plot - {}", ivp_name);
						let mut data = Vec::new();
						let mut names = Vec::new();
						let q = (self.delta_t / self.step_size) as usize;
						let t: Vec<f64> = (0..=q)
							.map(|q1| ivp.t0 + self.step_size * q1 as f64)
							.collect();
						let y0 = ivp.y0.clone();
						let extractor = &self.extractors[self.extractor_index].1;
						let f = |t: f64, y: &Vec<Vec<f64>>| {
							(ivp.differential_equation)(t, y, &ivp.parameters)
						};
						for (_, (solver_name, solver)) in self
							.solvers
							.iter_mut()
							.enumerate()
							.filter(|x| self.solvers_selected[x.0])
						{
							let mut y: Vec<Vec<Vec<f64>>> = t.iter().map(|_| y0.clone()).collect();
							for k in 0..q {
								y[k + 1] = solver.approximate(t[k], self.step_size, &f, &y[k]);
							}
							let y1: Vec<Vec<f64>> = y.into_iter().map(|x| extractor(&x)).collect();
							for k in 0..y1[0].len() {
								let y2: Vec<(f64, f64)> = t
									.iter()
									.enumerate()
									.map(|(i, t1)| (*t1, y1[i][k]))
									.collect();
								data.push(y2);
								names.push(Some(solver_name.clone() + &k.to_string()));
							}
						}
						let cbor_data: Vec<(Vec<f64>, Vec<f64>)> = (0..data.len())
							.map(|q| {
								(
									(0..data[q].len()).map(|k| data[q][k].0).collect(),
									(0..data[q].len()).map(|k| data[q][k].1).collect(),
								)
							})
							.collect();
						let encoded = to_vec(&cbor_data).expect("failed serialization");
						let mut file = File::create(cbor_path).unwrap();
						file.write_all(&encoded).unwrap();
						plot(
							&data,
							&(0..data.len())
								.map(|i| Palette99::pick(i).stroke_width(7))
								.collect(),
							&names,
							plot_path,
							title,
							("t".to_string(), "y".to_string()),
						);
					}
				}
			}
			if self.order_plot_error {
				ui.label(RichText::new("Can't plot unsolved ivp").color(Color32::RED));
			}
			if ui.button("quit").clicked() {
				ui.ctx().send_viewport_cmd(egui::ViewportCommand::Close);
			}
		});
	}
	fn on_exit(&mut self, _gl: Option<&eframe::glow::Context>) {}
}
