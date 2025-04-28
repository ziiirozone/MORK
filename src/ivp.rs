#![allow(non_snake_case)]

pub struct SolvedIVP {
    pub differential_equation: fn(f64, &Vec<Vec<f64>>, &Vec<f64>) -> Vec<f64>,
    pub solution: fn(f64, &Vec<f64>) -> Vec<Vec<f64>>,
    pub parameters: Vec<f64>,
}

impl SolvedIVP {
    pub fn f(&self, t: f64, y: &Vec<Vec<f64>>) -> Vec<f64> {
        (self.differential_equation)(t, y, &self.parameters)
    }
    pub fn solution(&self, t: f64) -> Vec<Vec<f64>> {
        (self.solution)(t, &self.parameters)
    }
}

pub enum IVPType {
    Solved(SolvedIVP),
    NotSolved(IVP),
}

pub struct IVP {
    pub differential_equation: fn(f64, &Vec<Vec<f64>>, &Vec<f64>) -> Vec<f64>,
    pub y0: Vec<Vec<f64>>,
    pub t0: f64,
    pub parameters: Vec<f64>,
}

impl IVP {
    pub fn f(&self, t: f64, y: &Vec<Vec<f64>>) -> Vec<f64> {
        (self.differential_equation)(t, y, &self.parameters)
    }
}

pub fn ivp_sin_cos() -> SolvedIVP {
    let differential_equation: for<'a, 'b> fn(f64, &'a Vec<Vec<f64>>, &'b Vec<f64>) -> _ =
        |_t, y, parameters| {
            let n = parameters[0] as usize;
            let q = n / 2;
            let r = n % 2;
            let mut s = 0.;
            for k in 0..q {
                s += (-1_f64).powi((q - k) as i32) * y[0][n - 2 * k - r - 1];
            }
            vec![s / (q as f64)]
        };
    let solution: for<'a, 'b> fn(f64, &'a Vec<f64>) -> _ = |t, parameters| -> Vec<Vec<f64>> {
        let n = parameters[0] as usize;
        let v = vec![t.sin() + t.cos(), t.cos() - t.sin()];
        vec![
            (1..=n)
                .map(|N| (-1_f64).powi((n - N) as i32 / 2) * v[(n - N) % 2])
                .collect(),
        ]
    };
    SolvedIVP {
        differential_equation,
        solution,
        parameters: vec![2.],
    }
}

pub fn ivp_exp() -> SolvedIVP {
    let differential_equation: for<'a, 'b> fn(f64, &'a Vec<Vec<f64>>, &'b Vec<f64>) -> _ =
        |_t, y, parameters| {
            let n = parameters[0] as usize;
            let s: f64 = y[0].iter().sum();
            vec![s / n as f64]
        };
    let solution: for<'a, 'b> fn(f64, &'a Vec<f64>) -> _ = |t, parameters| -> Vec<Vec<f64>> {
        let n = parameters[0] as usize;
        vec![vec![t.exp(); n]]
    };
    SolvedIVP {
        differential_equation,
        solution,
        parameters: vec![2.],
    }
}

pub fn exo_fatma() -> SolvedIVP {
    let differential_equation: for<'a, 'b> fn(f64, &'a Vec<Vec<f64>>, &'b Vec<f64>) -> _ =
        |t, y, _parameters| vec![9. * t - 3. * y[0][1] - 4. * y[0][0]];
    let solution: for<'a, 'b> fn(f64, &'a Vec<f64>) -> _ = |t, _parameters| -> Vec<Vec<f64>> {
        vec![vec![
            -2. * ((-t).exp() + 3. * (-3. * t).exp()) + 3.,
            2. * ((-t).exp() + (-3. * t).exp()) + 3. * t - 4.,
        ]]
    };
    SolvedIVP {
        differential_equation,
        solution,
        parameters: vec![],
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

pub fn first_coo_distance(x: &Vec<Vec<f64>>, y: &Vec<Vec<f64>>) -> Vec<f64> {
    y[0].iter()
        .zip(x[0].iter())
        .map(|(x0, y0)| (x0 - y0).abs())
        .collect()
}

pub fn highest_first_coo_distance(x: &Vec<Vec<f64>>, y: &Vec<Vec<f64>>) -> Vec<f64> {
    vec![(y[0][0] - x[0][0]).abs()]
}

pub fn lowest_first_coo_distance(x: &Vec<Vec<f64>>, y: &Vec<Vec<f64>>) -> Vec<f64> {
    vec![(y[0][y[0].len() - 1] - x[0][y[0].len() - 1]).abs()]
}
