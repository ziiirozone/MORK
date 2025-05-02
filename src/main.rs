pub mod experiment_picker;
pub mod graph;
pub mod ivp;
pub mod methods;
pub mod plot;

use crate::experiment_picker::*;
use crate::{ivp::*, methods::*};
use eframe::egui;

#[allow(non_snake_case)]

fn main() {
    let solvers: Vec<(&str, Box<dyn Solver>)> = vec![
        ("ENDMORK1", Box::new(ENDMORK1())),
        ("ENDMORK2", Box::new(ENDMORK2())),
        ("ENDMORK3", Box::new(ENDMORK3())),
        ("ENDMORK4_1", Box::new(ENDMORK4_1())),
        ("ENDMORK4_2", Box::new(ENDMORK4_2())),
        ("INDMORK1", Box::new(INDMORK1())),
        ("INDMORK2", Box::new(INDMORK2())),
        ("INDMORK3", Box::new(INDMORK3())),
        ("INDMORK3_1", Box::new(INDMORK3_1())),
        ("INDMORK4", Box::new(INDMORK4())),
        ("ERK1", Box::new(ERK1())),
        ("ERK2", Box::new(ERK2())),
        ("ERK3", Box::new(ERK3())),
        ("ERK4_1", Box::new(ERK4_1())),
        ("ERK4_2", Box::new(ERK4_2())),
        ("IRK1", Box::new(IRK1())),
        ("IRK2", Box::new(IRK2())),
        ("IRK3", Box::new(IRK3())),
        ("IRK3_1", Box::new(IRK3_1())),
        ("IRK4", Box::new(IRK4())),
    ];
    let solver = IRK3_2();
    println!("{:?}",solver.queue);
    let experiments = vec![
        ("sin - cos", IVPType::Solved(ivp_sin_cos()), vec!["n"]),
        ("exponential", IVPType::Solved(ivp_exp()), vec!["n"]),
        ("exo fatma", IVPType::Solved(exo_fatma()), vec![]),
    ];
    let extractors: Vec<(&str, Box<dyn Fn(&Vec<Vec<f64>>) -> Vec<f64>>)> = vec![
        ("First coordinate", Box::new(first_coo_extractor)),
        (
            "First coordinate highest derivative",
            Box::new(highest_first_coo_extractor),
        ),
        (
            "First coordinate lowest derivative",
            Box::new(lowest_first_coo_extractor),
        ),
    ];
    let distances: Vec<(
        &str,
        Box<dyn Fn(&Vec<Vec<f64>>, &Vec<Vec<f64>>) -> Vec<f64>>,
    )> = vec![
        ("First coordinate", Box::new(first_coo_distance)),
        (
            "First coordinate highest derivative",
            Box::new(highest_first_coo_distance),
        ),
        (
            "First coordinate lowest derivative",
            Box::new(lowest_first_coo_distance),
        ),
    ];
    let options = eframe::NativeOptions {
        viewport: egui::ViewportBuilder::default().with_inner_size([320.0, 640.0]),
        vsync: false,
        ..Default::default()
    };
    eframe::run_native(
        "Experiment selector",
        options,
        Box::new(|_cc| {
            Ok(Box::new(ExperimentPicker::new(
                experiments,
                solvers,
                extractors,
                distances,
            )))
        }),
    )
    .unwrap();
}
