use std::path::Path;

use plotters::{prelude::*, style::full_palette::GREY};

const RES4K: (u32, u32) = (3840, 2160);
const POINT_SIZE: u32 = 10;
const STROKE_WIDTH: u32 = 4;
const STROKE_LEGEND_WIDTH: u32 = 10;

fn plot_min(data: &Vec<Vec<(f64, f64)>>) -> ((f64, f64), (f64, f64)) {
    let mut y_min = f64::MAX;
    let mut y_max = f64::MIN;
    let mut x_min = f64::MAX;
    let mut x_max = f64::MIN;
    for j in 0..data.len() {
        for i in 0..data[j].len() {
            if data[j][i].1 < y_min {
                y_min = data[j][i].1;
            }
            if data[j][i].1 > y_max {
                y_max = data[j][i].1;
            }
            if data[j][i].0 < x_min {
                x_min = data[j][i].0;
            }
            if data[j][i].0 > x_max {
                x_max = data[j][i].0;
            }
        }
    }
    ((x_min, x_max), (y_min, y_max))
}

fn log_plot_key_points((min, max): (f64, f64)) -> Vec<f64> {
    ((min.log10()).floor() as i32..(max.log10()).ceil() as i32)
        .into_iter()
        .map(|power| 10_f64.powi(power))
        .collect()
}

pub fn plot<T>(
    data: &Vec<Vec<(f64, f64)>>,
    styles: &Vec<ShapeStyle>,
    names: &Vec<Option<String>>,
    path: &T,
    title: String,
    desc: (String, String),
) where
    T: AsRef<Path> + ?Sized,
{
    let ((x_min, x_max), (y_min, y_max)) = plot_min(data);
    let root = SVGBackend::new(&path, RES4K).into_drawing_area();
    root.fill(&WHITE).unwrap();
    let mut ctx = ChartBuilder::on(&root)
        .set_label_area_size(LabelAreaPosition::Left, 10.percent())
        .set_label_area_size(LabelAreaPosition::Bottom, 10.percent())
        .caption(title.clone(), ("sans-serif", 5.percent()))
        .margin_right(5.percent())
        .build_cartesian_2d(x_min..x_max, y_min..y_max)
        .unwrap();

    let line_style = ShapeStyle {
        color: GREY.into(),
        filled: true,
        stroke_width: 2,
    };

    ctx.configure_mesh()
        .bold_line_style(line_style)
        .x_label_style(("sans-serif", 4.percent()))
        .y_label_style(("sans-serif", 4.percent()))
        .x_desc(desc.0)
        .y_desc(desc.1)
        .draw()
        .unwrap();

    for ((d, s), legend) in data.iter().zip(styles).zip(names) {
        match legend {
            None => {
                ctx.draw_series(
                    LineSeries::new(
                        d.into_iter().map(|i| -> (f64, f64) { (i.0, i.1) }),
                        s.stroke_width(STROKE_WIDTH),
                    )
                    .point_size(POINT_SIZE),
                )
                .unwrap();
            }
            Some(n) => {
                ctx.draw_series(
                    LineSeries::new(
                        d.into_iter().map(|i| -> (f64, f64) { (i.0, i.1) }),
                        s.stroke_width(STROKE_WIDTH),
                    )
                    .point_size(POINT_SIZE),
                )
                .unwrap()
                .label(n)
                .legend(|(x, y)| {
                    PathElement::new(
                        vec![(x, y), (x + 30, y)],
                        s.stroke_width(STROKE_LEGEND_WIDTH),
                    )
                });
            }
        }
    }
    ctx.configure_series_labels()
        .position(SeriesLabelPosition::UpperRight)
        .border_style(&BLACK)
        .background_style(&WHITE.mix(0.8))
        .label_font(("sans-serif", 4.percent()))
        .draw()
        .unwrap();
    root.present().unwrap();
}

pub fn log_plot<T>(
    data: &Vec<Vec<(f64, f64)>>,
    styles: &Vec<ShapeStyle>,
    names: &Vec<Option<String>>,
    path: &T,
    title: String,
    desc: (String, String),
) where
    T: AsRef<Path> + ?Sized,
{
    let ((x_min, x_max), (y_min, y_max)) = plot_min(data);

    let root = SVGBackend::new(&path, RES4K).into_drawing_area();
    root.fill(&WHITE).unwrap();
    let mut ctx = ChartBuilder::on(&root)
        .set_label_area_size(LabelAreaPosition::Left, 10.percent())
        .set_label_area_size(LabelAreaPosition::Bottom, 10.percent())
        .margin_right(5.percent())
        .caption(title.clone(), ("sans-serif", 5.percent()))
        .build_cartesian_2d(
            (x_min..x_max)
                .log_scale()
                .with_key_points(log_plot_key_points((x_min, x_max))),
            (y_min..y_max)
                .log_scale()
                .with_key_points(log_plot_key_points((y_min, y_max))),
        )
        .unwrap();
    let line_style = ShapeStyle {
        color: GREY.into(),
        filled: true,
        stroke_width: 2,
    };
    ctx.configure_mesh()
        .bold_line_style(line_style)
        .x_label_style(("sans-serif", 4.percent()))
        .y_label_style(("sans-serif", 4.percent()))
        .x_desc(desc.0)
        .y_desc(desc.1)
        .draw()
        .unwrap();

    for ((d, s), legend) in data.iter().zip(styles).zip(names) {
        match legend {
            None => {
                ctx.draw_series(
                    LineSeries::new(
                        d.into_iter().map(|i| -> (f64, f64) { (i.0, i.1) }),
                        s.stroke_width(STROKE_WIDTH),
                    )
                    .point_size(POINT_SIZE),
                )
                .unwrap();
            }
            Some(n) => {
                ctx.draw_series(
                    LineSeries::new(
                        d.into_iter().map(|i| -> (f64, f64) { (i.0, i.1) }),
                        s.stroke_width(STROKE_WIDTH),
                    )
                    .point_size(POINT_SIZE),
                )
                .unwrap()
                .label(n)
                .legend(|(x, y)| {
                    PathElement::new(
                        vec![(x, y), (x + 20, y)],
                        s.stroke_width(STROKE_LEGEND_WIDTH),
                    )
                });
            }
        }
    }
    ctx.configure_series_labels()
        .position(SeriesLabelPosition::UpperLeft)
        .border_style(&BLACK)
        .label_font(("sans-serif", 4.percent()))
        .background_style(&WHITE.mix(0.8))
        .draw()
        .unwrap();
    root.present().unwrap();
}
