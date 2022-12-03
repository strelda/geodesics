fn kronecker_delta(a: i32, b: i32) -> f64 {
    if a==b {
        1 as f64
    } else {
        0 as f64
    }
}

fn je(j: i32, m: i32, s:i32) -> f64 {
    f64::sqrt((j*(j+1)-m*(m+s)) as f64 )
}


fn hamiltonian_value(l: f64, c: f64, jp: i32, mp: i32, j: i32, m: i32) -> f64  {
    let v1: f64 =  0.25*(je(1, j, m + 1)*je(1, j, m)*kronecker_delta(mp, m + 2) + je(-1, j, m - 1) *je(-1, j, m)*kronecker_delta(mp, m - 2) + kronecker_delta(mp, m)*(je(-1, j, m)*je(1, j, m - 1) + je(1, j, m)*je(-1, j, m + 1) ));
    
    let v2: f64 =  0.5*(m + j) as f64*(je(1, j, m) *kronecker_delta(mp, m + 1) + je(-1, j, m)*kronecker_delta(mp, m - 1)) + 0.5*(m + 1 + j) as f64 *je(1, j, m) *kronecker_delta(mp, m + 1) + 0.5*(m - 1 + j) as f64 * je(-1, j, m) *kronecker_delta(mp, m - 1);

    let v3: f64 = kronecker_delta(mp, m)*i32::pow(j + m, 2) as f64;


    m as f64 * kronecker_delta(jp, j) * kronecker_delta(mp, m) - kronecker_delta(jp, j)/(2.0 * j as f64)*(l * v1 + c * v2 + c*c* v3)
}


fn hamiltonian(nn: u32, l: f64, c: f64) -> f64 {
    let jj = nn as f64 /2.0; 
    let dim = nn + 1;

    // hamiltonian_value(jj, mp, jj, m) for mp in (-jj,jj) m in (-jj,jj)
}