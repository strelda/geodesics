#![allow(unused_imports)]
#![allow(dead_code)]
use ndarray::{array, Array1, Array2, Axis};
use ndarray_linalg::{Eigh, UPLO};
use std::f64;

fn h(l:f64, c:f64)-> Array2<f64> {
    let dvat: f64 = 2.*f64::sqrt(3.);
    array![
        [-(l+6.)/4. , -c/dvat    ,-l/dvat    , 0.],
        [-c/dvat    , (-6.-7.*l-4.*c.powi(2))/12.  ,-c   ,-l/dvat],
        [-l/dvat    , -c         , (6.-7.*l-16.*c.powi(2))/12.  ,-5.*c/dvat],
        [0.,-l/dvat , -5.*c/dvat          , 3./2. - l/4.-3.*c.powi(2)]
    ]
}

fn dh1(_l:f64, _c:f64)-> Array2<f64> {
    let fr: f64 = -1./(2.*f64::sqrt(3.));
    array![
        [-1./4. , 0.     ,fr      ,0.    ],
        [0.     ,-7./12. ,0.      ,fr    ],
        [fr     ,0.      ,-7./12. ,0.    ],
        [0.     ,fr      ,0.      ,-1./4.]
    ]
}
fn dh2(_l:f64, c:f64)-> Array2<f64> {
    let fr: f64 = -1./(2.*f64::sqrt(3.));
    array![
        [0.     , fr      ,0.       ,0.    ],
        [fr     ,-2.*c/3. ,-1.      ,0.    ],
        [0.     ,-1.      ,-8.*c/3. ,5.*fr ],
        [0.     ,0.       ,5.*fr    ,-6.*c ]
    ]
}

// fn Mk(eigv,eigs,k){

// }

fn g(l:f64, c:f64) -> Array2<f64> {
    let (eigs, eigv) = h(l,c).eigh(UPLO::Lower).unwrap();
    let dim = eigv.select(Axis(1),&[0]).len();

    let mut g11: f64 = 0.0;
    let mut g12: f64 = 0.0;
    let mut g22: f64 = 0.0;
    for k in 1..dim{
        let m1 = &eigv.select(Axis(1),&[0]).reversed_axes()
            .dot(&dh1(l,c))
            .dot(&eigv.select(Axis(1),&[k])
                )[[0,0]];
        let m1t = &eigv.select(Axis(1),&[k]).reversed_axes()
            .dot(&dh1(l,c))
            .dot(&eigv.select(Axis(1),&[0])
                )[[0,0]];
        let m2 = &eigv.select(Axis(1),&[0]).reversed_axes()
            .dot(&dh2(l,c))
            .dot(&eigv.select(Axis(1),&[k])
                )[[0,0]];
        let m2t = &eigv.select(Axis(1),&[k]).reversed_axes()
            .dot(&dh2(l,c))
            .dot(&eigv.select(Axis(1),&[0])
                )[[0,0]];
        g11 += m1*m1t/(&eigs[0]-&eigs[k]).powi(2);
        g12 += m1*m2t/(&eigs[0]-&eigs[k]).powi(2);
        g22 += m2*m2t/(&eigs[0]-&eigs[k]).powi(2);
    }
    array![
        [g11,g12],
        [g12,g22]
    ]
}
fn main(){
    dbg!(g(2.0,3.0));
}