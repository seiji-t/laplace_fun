use std::f64::consts::PI;

fn pt(l: usize, m: usize) -> usize {
    m + (l * (l + 1)) / 2
}

fn precompute_ab(ll: usize, a: &mut Vec<f64>, b: &mut Vec<f64>) {
    for l in 2..=ll {
        let ls = (l * l) as f64;
        let lm1s = ((l - 1) * (l - 1)) as f64;
        for m in 0..(l - 1) {
            let ms = (m * m) as f64;
            a[pt(l, m)] = ((4.0 * ls - 1.0) / (ls - ms)).sqrt();
            b[pt(l, m)] = -((lm1s - ms) / (4.0 * lm1s - 1.0)).sqrt();
        }
    }
}

fn compute_p(
    l_max: usize,
    a: &Vec<f64>,
    b: &Vec<f64>,
    p: &mut Vec<f64>,
    x: f64,
) {
    let sintheta = (1.0 - x * x).sqrt();
    let mut temp = (0.5 / PI).sqrt();
    p[pt(0, 0)] = temp;

    if l_max > 0 {
        let sqrt3 = 1.7320508075688772935;
        p[pt(1, 0)] = x * sqrt3 * temp;

        let sqrt3div2 = -1.2247448713915890491;
        temp = sqrt3div2 * sintheta * temp;
        p[pt(1, 1)] = temp;

        for l in 2..=l_max {
            for m in 0..(l - 1) {
                p[pt(l, m)] = a[pt(l, m)] * (x * p[pt(l - 1, m)] + b[pt(l, m)] * p[pt(l - 2, m)]);
            }
            p[pt(l, l - 1)] = x * ((2 * (l - 1) + 3) as f64).sqrt() * temp;
            temp = -((1.0 + 0.5 / (l as f64)).sqrt()) * sintheta * temp;
            p[pt(l, l)] = temp;
        }
    }
}


#[cfg(test)]
mod tests{
    use super::*;
    use gauss_quad::GaussLegendre;
    #[test]
    fn test_calc_alp_01(){
        let l_max = 5;
        let n = (l_max + 1) * (l_max + 2) / 2;

        let mut a = vec![0.0; n];
        let mut b = vec![0.0; n];
        let mut p = vec![0.0; n];

        let x = 0.1;

        precompute_ab(l_max, &mut a, &mut b);
        compute_p(l_max, &a, &b, &mut p, x);

        assert!( (p[pt(0,0)] - 0.3989422804014327).abs() <= 1e-10);
        assert!( (p[pt(1,0)] - 0.690988298942671 * x).abs() <= 1e-10);
        assert!( (p[pt(1,1)] - 0.4886025119029199 * (-1.0*(1.0 - x*x).sqrt())).abs() <= 1e-10);
    }

    #[test]
    fn test_normalization(){
    let number_of_nodes = 64;
    let quad = GaussLegendre::new(number_of_nodes).expect("Error calculating Gauss-Legendre quadrature.");

    let l_max = 5;
    let n = (l_max + 1) * (l_max + 2) / 2;

    let mut a = vec![0.0; n];
    let mut b = vec![0.0; n];
    let mut p = vec![0.0; n];

    let norm = 1.0/PI;
    
    precompute_ab(l_max, &mut a, &mut b);
    let mut temp = 0.0;
    for (mu_val, weight) in quad.iter(){
        println!("{:?},{:?}", mu_val, weight);
        temp += weight;
    } 
    assert!((temp-2.0).abs() <= 1e-8);
    for n in 0..=l_max{
        for m in 0..=n{
            let mut inner_prod = 0.0;
            for (mu_val, weight) in quad.iter(){
                compute_p(l_max, &a, &b, &mut p, *mu_val);
                inner_prod += weight * p[pt(n,m)] * p[pt(n,m)];
            }
            println!("{:?}", inner_prod);
            assert!((inner_prod-norm).abs() <= 1e-8);
        }
    }
    }

    #[test]
    fn test_orthogonality(){
    let number_of_nodes = 64;
    let quad = GaussLegendre::new(number_of_nodes).expect("Error calculating Gauss-Legendre quadrature.");

    let l_max = 5;
    let n = (l_max + 1) * (l_max + 2) / 2;

    let mut a = vec![0.0; n];
    let mut b = vec![0.0; n];
    let mut p = vec![0.0; n];
    
    precompute_ab(l_max, &mut a, &mut b);

    for n1 in 0..=l_max{
        for m1 in 0..=n1{
            for n2 in 0..=l_max{
                for m2 in 0..=n2{
                    if n1 != n2 && m1 == m2{
                        let mut inner_prod = 0.0;
                        for (mu_val, weight) in quad.iter(){
                            compute_p(l_max, &a, &b, &mut p, *mu_val);
                            inner_prod += weight * p[pt(n1,m1)] * p[pt(n2,m2)];
                        }
                        println!("{:?},{:?};{:?},{:?} = {:?}",n1,m1, n2,m2, inner_prod);
                        assert!((inner_prod-0.0).abs() <= 1e-8);
                        }
                    }
                }
            }
        }
    }
}