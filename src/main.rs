#![allow(unused_imports)]
#![allow(dead_code)]

mod functions;
use functions::{g};

use ndarray::{array,Array2};
use itertools_num::linspace;
use plotly::common::{ColorScale, ColorScalePalette, Title};
use plotly::contour::Contours;
use plotly::{Contour, HeatMap, Layout, Plot};
use easy_ml::differentiation;

fn heat_map(z: Vec<Vec<f64>>) {
    let trace = HeatMap::new_z(z);
    let mut plot = Plot::new();
    plot.add_trace(trace);    
    plot.show();
}

fn det2x2(m: Array2<f64>)->f64{
        m[[0,0]]*m[[1,1]]-m[[0,1]]*m[[1,0]]
}


fn contour_plot_detg(n: usize, xc: f64, yc: f64, dist: f64) {
    let x: Vec<f64> = linspace(xc-dist,xc+dist,n).collect();
    let y: Vec<f64> = linspace(yc-dist,yc+dist,n).collect();
    let mut z: Vec<Vec<f64>> = Vec::new();

    for yi in 0..n {
        let mut row = Vec::<f64>::new();
        for xi in 0..n {
            let zv = det2x2(g(x[xi],y[yi])).atan();
            row.push(zv);
        }
        z.push(row);
    }

    // let layout = Layout::new().title(Title::new("Metric tensor determinant"));
    let trace = HeatMap::new(x,y,z).color_scale(ColorScale::Palette(ColorScalePalette::Jet)).auto_color_scale(false);

    let mut plot = Plot::new();

    plot.add_trace(trace);
    plot.show();
}


fn main(){
    // heat_map(z);
    
    contour_plot_detg(10, -0.5, f64::sqrt(3./5.), 1.);


}