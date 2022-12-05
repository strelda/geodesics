#![allow(unused_imports)]
#![allow(dead_code)]

mod functions;
use std::array;

use functions::{g};

use nalgebra::Matrix2;
use ndarray::{array,Array2};
use itertools_num::linspace;
use plotly::common::{ColorScale, ColorScalePalette, Title};
use plotly::{HeatMap, Layout, Plot};
use rgsl::numerical_differentiation::{self, deriv_central};
use nalgebra::base::SquareMatrix;

fn det2x2(m: &Array2<f64>)->f64{
        m[[0,0]]*m[[1,1]]-m[[0,1]]*m[[1,0]]
}
fn inverse2x2(m: &Array2<f64>) -> Array2<f64>{
    let det: f64 = det2x2(m);
    array![
        [m[[1,1]]/det,-m[[0,1]]/det],
        [-m[[1,0]]/det,m[[0,0]]/det]
    ]
}

fn heat_plot_detg(n: usize, xc: f64, yc: f64, dist: f64) {
    let x: Vec<f64> = linspace(xc-dist,xc+dist,n).collect();
    let y: Vec<f64> = linspace(yc-dist,yc+dist,n).collect();
    let mut z: Vec<Vec<f64>> = Vec::new();

    for yi in 0..n {
        let mut row = Vec::<f64>::new();
        for xi in 0..n {
            let zv = det2x2(&g(x[xi],y[yi])).atan();
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


fn gfixc_11(c: f64) -> impl Fn(f64)->f64{
    move |l: f64| g(l,c)[[0,0]]
}
fn gfixc_12(c: f64) -> impl Fn(f64)->f64{
    move |l: f64| g(l,c)[[0,1]]
}
fn gfixc_22(c: f64) -> impl Fn(f64)->f64{
    move |l: f64| g(l,c)[[1,1]]
}
fn gfixl_11(l: f64) -> impl Fn(f64)->f64{
    move |c: f64| g(l,c)[[0,0]]
}
fn gfixl_12(l: f64) -> impl Fn(f64)->f64{
    move |c: f64| g(l,c)[[0,1]]
}
fn gfixl_22(l: f64) -> impl Fn(f64)->f64{
    move |c: f64| g(l,c)[[1,1]]
}

fn g_dl(l: f64, c: f64) -> Array2<f64>{
    let deriv_step: f64 = 1e-5;
    let g11_dl: f64 = deriv_central(gfixc_11(c), l, deriv_step).1;
    let g12_dl: f64 = deriv_central(gfixc_12(c), l, deriv_step).1;
    let g22_dl: f64 = deriv_central(gfixc_22(c), l, deriv_step).1;
    array![
        [g11_dl,g12_dl],
        [g12_dl,g22_dl]
        ]
}
fn g_dc(l: f64, c: f64) -> Array2<f64>{
    let deriv_step: f64 = 1e-5;
    let g11_dc: f64 = deriv_central(gfixl_11(l), c, deriv_step).1;
    let g12_dc: f64 = deriv_central(gfixl_12(l), c, deriv_step).1;
    let g22_dc: f64 = deriv_central(gfixl_22(l), c, deriv_step).1;
    array![
        [g11_dc,g12_dc],
        [g12_dc,g22_dc]
        ]
}
fn gamma(l: f64, c: f64) -> Array2<f64>{
    let g_dl: Array2<f64> = g_dl(l,c); 
    let g_dc: Array2<f64> = g_dc(l,c);
    let g_inv_2: Array2<f64> = 0.5*inverse2x2(&g(l,c)); 
    //make it matrix multiplication
    let g11 = g_inv_2[[0,0]] * (2.*g_dl[[0,0]] - g_dl[[0,0]]) 
                +  g_inv_2[[0,1]] * (2.*g_dl[[0,1]] - g_dc[[0,0]]);
    let g12 = g_inv_2[[0,0]] * (g_dc[[0,0]]) 
                +  g_inv_2[[0,1]] * (g_dl[[1,1]]);
    let g13 = g_inv_2[[0,0]] * (2.*g_dc[[1,0]]- g_dl[[1,1]]) 
                +  g_inv_2[[0,1]] * (g_dc[[1,1]]);
    let g21 = g_inv_2[[1,0]] * (2.*g_dl[[0,0]] - g_dl[[0,0]]) 
                +  g_inv_2[[1,1]] * (2.*g_dl[[0,1]] - g_dc[[0,0]]);
    let g22 = g_inv_2[[1,0]] * (g_dc[[0,0]]) 
                +  g_inv_2[[1,1]] * (g_dl[[1,1]]);
    let g23 = g_inv_2[[1,0]] * (2.*g_dc[[1,0]]- g_dl[[1,1]]) 
                +  g_inv_2[[1,1]] * (g_dc[[1,1]]);
    
    array![
        [g11,g12,g13],
        [g21,g22,g23]]
}

fn main(){
    // heat_plot_detg(10, -0.5, f64::sqrt(3./5.), 1.);

    // dbg!(deriv_central(gfixl(1.)[[0,0]], 1., 0.00001));
    dbg!(gamma(2.,3.));
    // deriv_central(g(1.), 1., 0.00001);
}