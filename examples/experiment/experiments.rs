#![allow(non_snake_case)]

pub trait Experiment {
    fn differential_equation(&self,t:f64, y: &Vec<Vec<f64>>) -> Vec<f64>;
    fn is_solved(&self) -> bool;
    fn solution(&self,t : f64) -> Option<Vec<Vec<f64>>>; // Should return None if not solution
    fn initial_values(&self) -> (f64,Vec<Vec<f64>>); // t0 and y0
    fn get_parameters(&self) -> (Vec<f64>,Vec<String>);
    fn change_parameters(&mut self, parameters : Vec<f64>);
    fn apply_parameters(&mut self); // If changing parameters involves heavy computations this function will be called before the experiment is used as a trigger for those computations
    fn name(&self) -> String;
}

pub struct SolvedIVP {
    pub differential_equation: fn(f64, &Vec<Vec<f64>>, &Vec<f64>,f64,&Vec<f64>) -> Vec<f64>, // t,y,parameters,t0,cache
    pub solution: fn(f64, &Vec<f64>,f64, &Vec<f64>) -> Vec<Vec<f64>>, // t, parameter,t0, cache
    pub compute_cache: fn(&Vec<f64>,f64) -> Vec<f64>, // parameters,t0
    pub t0 : f64,
    pub parameters: Vec<f64>,
    pub parameters_names : Vec<String>,
    pub ivp_name : String,
    pub cache : Vec<f64>,
}

impl Experiment for SolvedIVP {
    fn differential_equation(&self,t:f64, y: &Vec<Vec<f64>>) -> Vec<f64> {
        (self.differential_equation)(t,y,&self.parameters,self.t0,&self.cache)
    }
    fn is_solved(&self) -> bool {
        true
    }
    fn apply_parameters(&mut self) {
        self.cache = (self.compute_cache)(&self.parameters,self.t0)
    }
    fn solution(&self,t : f64) -> Option<Vec<Vec<f64>>> {
        Some((self.solution)(t,&self.parameters,self.t0,&self.cache))
    }
    fn change_parameters(&mut self, mut parameters : Vec<f64>) {
        let t0 = parameters.remove(0);
        self.parameters = parameters;
        self.t0 = t0;
    }
    fn get_parameters(&self) -> (Vec<f64>,Vec<String>) {
        let mut parameters = self.parameters.clone();
        parameters.insert(0,self.t0);
        let mut parameters_name = self.parameters_names.clone();
        parameters_name.insert(0,"t0".to_string());
        (parameters,parameters_name)
    }
    fn initial_values(&self) -> (f64,Vec<Vec<f64>>) {
        (self.t0,(self.solution)(self.t0,&self.parameters,self.t0,&self.cache))
    }
    fn name(&self) -> String {
        self.ivp_name.clone()
    }
}
pub struct IVP {
    pub differential_equation: fn(f64, &Vec<Vec<f64>>, &Vec<f64>,f64,&Vec<f64>) -> Vec<f64>, // t,y,parameters,t0,cache
    pub compute_cache: fn(&Vec<f64>,f64) -> Vec<f64>, // parameters,t0
    pub t0 : f64,
    pub y0 : Vec<Vec<f64>>,
    pub parameters: Vec<f64>,
    pub parameters_names : Vec<String>,
    pub ivp_name : String,
    pub cache : Vec<f64>,
}

impl Experiment for IVP {
    fn differential_equation(&self,t:f64, y: &Vec<Vec<f64>>) -> Vec<f64> {
        (self.differential_equation)(t,y,&self.parameters,self.t0,&self.cache)
    }
    fn is_solved(&self) -> bool {
        false
    }
    fn apply_parameters(&mut self) {
        self.cache = (self.compute_cache)(&self.parameters,self.t0)
    }
    fn solution(&self,_t : f64) -> Option<Vec<Vec<f64>>> {
        None
    }
    fn change_parameters(&mut self, parameters : Vec<f64>) {
        self.parameters = parameters;
    }
    fn get_parameters(&self) -> (Vec<f64>,Vec<String>) {
        (self.parameters.clone(),self.parameters_names.clone())
    }
    fn initial_values(&self) -> (f64,Vec<Vec<f64>>) {
        (self.t0,self.y0.clone())
    }
    fn name(&self) -> String {
        self.ivp_name.clone()
    }
}

pub fn experiment_sin_cos() -> SolvedIVP {
    let differential_equation=
        |_, y: &Vec<Vec<f64>>, parameters: &Vec<f64>,_t0: f64,_cache: &Vec<f64>| {
            let n = parameters[0] as usize;
            let q = n / 2;
            let r = n % 2;
            let mut s = 0.;
            for k in 0..q {
                s += (-1_f64).powi((q - k) as i32) * y[0][n - 2 * k - r - 1];
            }
            vec![s / (q as f64)]
        };
    let solution = |t: f64, parameters : &Vec<f64>,_to: f64,_cache: &Vec<f64>| -> Vec<Vec<f64>> {
        let n = parameters[0] as usize;
        let v = vec![t.sin() + t.cos(), t.cos() - t.sin()];
        vec![
            (1..=n)
                .map(|N| (-1_f64).powi((n - N) as i32 / 2) * v[(n - N) % 2])
                .collect(),
        ]
    };
    let compute_cache = |_parameters: &Vec<f64>,_t0: f64| Vec::new();
    SolvedIVP {
        differential_equation,
        compute_cache,
        t0:0.,
        solution,
        parameters: vec![2.],
        parameters_names: vec!["n".to_string()],
        ivp_name: "sin + cos".to_string(),
        cache : Vec::new()
    }
}


pub fn experiment_exp() -> SolvedIVP {
    let differential_equation=
        |_, y: &Vec<Vec<f64>>, parameters: &Vec<f64>,_t0: f64,_cache: &Vec<f64>| {
            let n = parameters[0] as usize;
            let s: f64 = y[0].iter().sum();
            vec![s / n as f64]
        };
    let solution = |t: f64, parameters : &Vec<f64>,_to: f64,_cache: &Vec<f64>| -> Vec<Vec<f64>> {
        let n = parameters[0] as usize;
        vec![vec![t.exp(); n]]
    };
    let compute_cache = |_parameters: &Vec<f64>,_t0: f64| Vec::new();
    SolvedIVP {
        differential_equation,
        compute_cache,
        t0:0.,
        solution,
        parameters: vec![2.],
        parameters_names: vec!["n".to_string()],
        ivp_name: "exp".to_string(),
        cache : Vec::new()
    }
}

pub fn first_coo_extractor(y: &Vec<Vec<f64>>) -> Vec<f64> {
    y[0].clone()
}

pub fn highest_first_coo_extractor(y: &Vec<Vec<f64>>) -> Vec<f64> {
    vec![y[0][0]]
}

pub fn lowest_first_coo_extractor(y: &Vec<Vec<f64>>) -> Vec<f64> {
    vec![y[0][y[0].len() - 1]]
}

pub fn first_coo_metric(x: &Vec<Vec<f64>>, y: &Vec<Vec<f64>>) -> Vec<f64> {
    y[0].iter()
        .zip(x[0].iter())
        .map(|(x0, y0)| (x0 - y0).abs())
        .collect()
}

pub fn highest_first_coo_metric(x: &Vec<Vec<f64>>, y: &Vec<Vec<f64>>) -> Vec<f64> {
    vec![(y[0][0] - x[0][0]).abs()]
}

pub fn lowest_first_coo_metric(x: &Vec<Vec<f64>>, y: &Vec<Vec<f64>>) -> Vec<f64> {
    vec![(y[0][y[0].len() - 1] - x[0][y[0].len() - 1]).abs()]
}
