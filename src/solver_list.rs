
pub mod NDMORK_methods {
    use crate::solvers::NDMORK;

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
            vec![(2_f64/3_f64).powi(N as i32), 0.],
            vec![(2. * N as f64 - 1.) / (2.*(1. + N as f64)), 3. / (2.*(1. + N as f64))],
        ]
    }

    pub fn ENDMORK2bis_nodes() -> Vec<f64> {
        vec![0., 2./3., 1.]
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

pub mod RK_methods {
    use crate::solvers::RK;

    pub fn ERK1_weights() -> Vec<Vec<f64>> {
        vec![vec![0.], vec![1.]]
    }

    pub fn ERK1_nodes() -> Vec<f64> {
        vec![0., 1.]
    }

    pub fn ERK1() -> RK {
        RK::new(ERK1_weights(), ERK1_nodes())
    }

    pub fn IRK1_weights() -> Vec<Vec<f64>> {
        vec![vec![1.], vec![1.]]
    }

    pub fn IRK1_nodes() -> Vec<f64> {
        vec![1., 1.]
    }

    pub fn IRK1() -> RK {
        RK::new(IRK1_weights(), IRK1_nodes())
    }

    pub fn ERK2_weights() -> Vec<Vec<f64>> {
        vec![vec![0., 0.], vec![0.5, 0.], vec![0., 1.]]
    }

    pub fn ERK2_nodes() -> Vec<f64> {
        vec![0., 0.5, 1.]
    }

    pub fn ERK2() -> RK {
        RK::new(ERK2_weights(), ERK2_nodes())
    }

    pub fn IRK2_weights() -> Vec<Vec<f64>> {
        vec![vec![0.5], vec![1.]]
    }

    pub fn IRK2_nodes() -> Vec<f64> {
        vec![0.5, 1.]
    }

    pub fn IRK2() -> RK {
        RK::new(IRK2_weights(), IRK2_nodes())
    }

    pub fn IRK3_weights() -> Vec<Vec<f64>> {
        vec![vec![0., 0.], vec![0., 1.], vec![0.5, 0.5]]
    }

    pub fn IRK3_nodes() -> Vec<f64> {
        vec![0., 1., 1.]
    }

    pub fn IRK3() -> RK {
        RK::new(IRK3_weights(), IRK3_nodes())
    }

    pub fn IRK3_1_weights() -> Vec<Vec<f64>> {
        vec![vec![0., 0.], vec![0.5, 0.5], vec![0.5, 0.5]]
    }

    pub fn IRK3_1_nodes() -> Vec<f64> {
        vec![0., 1., 1.]
    }

    pub fn IRK3_1() -> RK {
        RK::new(IRK3_1_weights(), IRK3_1_nodes())
    }

    pub fn ERK3_weights() -> Vec<Vec<f64>> {
        vec![
            vec![0., 0., 0.],
            vec![1. / 3., 0., 0.],
            vec![0., 2. / 3., 0.],
            vec![1. / 4., 0., 3. / 4.],
        ]
    }

    pub fn ERK3_nodes() -> Vec<f64> {
        vec![0., 1. / 3., 2. / 3., 1.]
    }

    pub fn ERK3() -> RK {
        RK::new(ERK3_weights(), ERK3_nodes())
    }

    pub fn ERK4_1_weights() -> Vec<Vec<f64>> {
        {
            vec![
                vec![0., 0., 0., 0.],
                vec![0.5, 0., 0., 0.],
                vec![0., 0.5, 0., 0.],
                vec![0., 0., 1., 0.],
                vec![1. / 6., 1. / 3., 1. / 3., 1. / 6.],
            ]
        }
    }

    pub fn ERK4_1_nodes() -> Vec<f64> {
        vec![0., 0.5, 0.5, 1., 1.]
    }

    pub fn ERK4_1() -> RK {
        RK::new(ERK4_1_weights(), ERK4_1_nodes())
    }

    pub fn ERK4_2_weights() -> Vec<Vec<f64>> {
        vec![
            vec![0., 0., 0., 0.],
            vec![0.5, 0., 0., 0.],
            vec![0.25, 0.25, 0., 0.],
            vec![0., -1., 2., 0.],
            vec![1. / 6., 0., 2. / 3., 1. / 6.],
        ]
    }

    pub fn ERK4_2_nodes() -> Vec<f64> {
        vec![0., 0.5, 0.5, 1., 1.]
    }

    pub fn ERK4_2() -> RK {
        RK::new(ERK4_2_weights(), ERK4_2_nodes())
    }

    pub fn IRK4_weights() -> Vec<Vec<f64>> {
        let sqrt3 = 3_f64.sqrt();
        vec![
            vec![1. / 4., 1. / 4. - sqrt3 / 6.],
            vec![1. / 4. + sqrt3 / 6., 1. / 4.],
            vec![0.5, 0.5],
        ]
    }

    pub fn IRK4_nodes() -> Vec<f64> {
        let sqrt3 = 3_f64.sqrt();
        vec![0.5 - sqrt3 / 6., 0.5 + sqrt3 / 6., 1.]
    }

    pub fn IRK4() -> RK {
        RK::new(IRK4_weights(), IRK4_nodes())
    }
}
