use scirs2_special::{legendre_assoc};
use std::f64::consts::{PI};
use std::fs::File;
use std::io::prelude::*;
use gauss_quad::GaussLegendre;
use std::str::FromStr;
use std::env;

fn main() -> std::io::Result<()>{
    let c = 1.0;
    let t_period = 2.0*PI/c;

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
        eprintln!("Usage: laplacian_fun NMAX NT");
        std::process::exit(1);
    }

    let nmax = params[0];
    let nt = params[1];

    for i in 0..n_mu{
        phi.push(-PI/2.0 + (i as f64) * dphi);
    }

    for i in 0..n_theta{
        theta.push(0.0 + (i as f64) * dtheta);
    }

    for n in 0..nmax{
        for m in 0..n+1{
            println!("n,m = {:?},{:?}",n,m);
            let omega_nm = c * ((n as f64) * ((n + 1) as f64)).sqrt();

            for i_t in 0..nt{ 
                let filename = format!("scripts/data/{:02}_{:02}_{:04}.txt", n,m,i_t);
                let mut file = File::create(filename)?;
                let t = 0.0 + (i_t as f64) * t_period/((nt - 1 )as f64);
                for mu_val in mu.clone(){
                    for j in 0..n_theta{
                        let theta_val = theta[j];
                        let pmn_val = legendre_assoc(n,m as i32,*mu_val);
                        let cos_theta = (m as f64 * theta_val).cos();
                        let cos_omega = (omega_nm * t).cos();
                        let u_val = radius + amp * pmn_val * cos_theta * cos_omega;
                        write!(file, "{:?},{:?},{:?}\n",theta_val, *mu_val, u_val)?;
                    }
                }
            }
        }
    }
    Ok(())
}

