pub mod experiment_picker;
pub mod experiments;
pub mod plot;

use std::env;
use std::fs::create_dir;
use std::io::ErrorKind;

use crate::experiment_picker::*;
use crate::experiments::*;
use MORK::solvers::NDMORK_methods::*;
use MORK::solvers::RK_methods::*;
use MORK::solvers::Solver;
use eframe::egui;

#[allow(non_snake_case)]

fn main() {
    let solvers: Vec<(&str, Box<dyn Solver>)> = vec![
        ("ENDMORK1", Box::new(ENDMORK1())),
        ("ENDMORK2", Box::new(ENDMORK2())),
        ("ENDMORK2 bis", Box::new(ENDMORK2bis())),
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
    let experiments: Vec<Box<dyn Experiment>> = vec![
        Box::new(experiment_sin_cos()),
        Box::new(experiment_exp()),
        Box::new(experiment_linear_single_root()),
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
    let metrics: Vec<(
        &str,
        Box<dyn Fn(&Vec<Vec<f64>>, &Vec<Vec<f64>>) -> Vec<f64>>,
    )> = vec![
        ("First coordinate", Box::new(first_coo_metric)),
        (
            "First coordinate highest derivative",
            Box::new(highest_first_coo_metric),
        ),
        (
            "First coordinate lowest derivative",
            Box::new(lowest_first_coo_metric),
        ),
    ];
    let options = eframe::NativeOptions {
        viewport: egui::ViewportBuilder::default().with_inner_size([600.0, 800.0]),
        vsync: false,
        ..Default::default()
    };
    let mut path = env::current_dir().expect("couldn't get current directory");
    path.push("output");
    match create_dir(&path) {
        Ok(_) => {
            println!("Created folder \"output\" at {path:?}");
        }
        Err(error) => {
            if error.kind() == ErrorKind::AlreadyExists {
                println!("Folder \"output\" already existed at {path:?}");
            } else {
                panic!("Couldn't create folder \"output\" at {path:?}, error : \n {error:?}")
            }
        }
    };
    eframe::run_native(
        "Experiment selector",
        options,
        Box::new(|_cc| {
            Ok(Box::new(ExperimentPicker::new(
                path,
                experiments,
                solvers,
                extractors,
                metrics,
            )))
        }),
    )
    .unwrap();
}
