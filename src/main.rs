use std::f64::consts::{PI};
use std::fs::File;
use std::io::prelude::*;
use gauss_quad::GaussLegendre;
use std::str::FromStr;
use std::env;
mod alp;

fn main() -> std::io::Result<()>{
    let c = 1.0;

    let amp = 0.5;
    let radius = 1.0;

    let n_mu: usize= 64;
    let n_theta: usize= 128;

    let quad = GaussLegendre::new(n_mu).expect("Error calculating Gauss-Legendre quadrature.");

    let mu = quad.nodes();

    let mut theta: Vec<f64> = Vec::new();
    let mut phi: Vec<f64> = Vec::new();

    let dtheta: f64 = 2.0*PI/((n_theta - 1) as f64);
    let dphi: f64 = 1.0*PI/((n_mu - 1) as f64);

    let mut params = Vec::new();
    for arg in env::args().skip(1) {
        params.push(usize::from_str(&arg)
        .expect("Error parsing arguments."));
    }

    if params.len() == 0{
        eprintln!("Usage: laplacian_fun NMAX");
        std::process::exit(1);
    }

    let nmax = params[0];

    for i in 0..n_mu{
        phi.push(-PI/2.0 + (i as f64) * dphi);
    }

    for i in 0..n_theta{
        theta.push(0.0 + (i as f64) * dtheta);
    }

    let n_size = (nmax + 1) * (nmax + 2) / 2;

    let mut a = vec![0.0; n_size];
    let mut b = vec![0.0; n_size];
    let mut p = vec![0.0; n_size];
    
    alp::precompute_ab(nmax, &mut a, &mut b);

    for n in 0..nmax{
        for m in 0..n+1{
            println!("n,m = {:?},{:?}",n,m);
            let omega_nm = c * ((n as f64) * ((n + 1) as f64)).sqrt();
            let filename = format!("data/{:02}_{:02}.txt", n,m);
            let mut file = File::create(filename)?;
            write!(file, "{:?},{:?},{:?},{:?},{:?}\n", omega_nm,c,radius,n_mu,n_theta)?;
            for mu_val in mu.clone(){
                for j in 0..n_theta{
                    let theta_val = theta[j];
                    alp::compute_p(nmax, &a, &b, &mut p, *mu_val);
                    let pmn_val = p[alp::pt(n,m)];
                    let cos_theta = (m as f64 * theta_val).cos();
                    let u_val = radius + amp * pmn_val * cos_theta;
                    write!(file, "{:?},{:?},{:?}\n",theta_val, *mu_val, u_val)?;
                }
            }
        }
    }
    Ok(())
}

