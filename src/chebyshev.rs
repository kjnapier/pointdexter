use nalgebra::{DMatrix, DVector, Matrix3, Matrix3xX, Unit, Vector3, matrix};


#[inline(never)]
pub fn fit_chebyshev_direct(xs: &[f64], ys: &[f64], degree: usize) -> Vec<f64> {
    
    let n = xs.len();
    let mut vander = DMatrix::zeros(n, degree + 1);
    let mut y_vec = DVector::from_element(n, 0.0);

    for (i, &x) in xs.iter().enumerate() {
        let mut t_prev = 1.0;
        let mut t_curr = x;
        vander[(i, 0)] = t_prev;

        if degree >= 1 {
            vander[(i, 1)] = t_curr;
        }
        for k in 2..=degree {
            let t_next = 2.0 * x * t_curr - t_prev;
            vander[(i, k)] = t_next;
            t_prev = t_curr;
            t_curr = t_next;
        }
        y_vec[i] = ys[i];
    }

    let vt_v = &vander.transpose() * &vander;
    let vt_y = &vander.transpose() * &y_vec;

    let coeffs = vt_v.lu().solve(&vt_y).expect("Singular matrix");

    coeffs.iter().copied().collect()
}


// #[inline(never)]
pub fn chebyshev_eval(c: &[f64], x: f64) -> f64 {
    if c.is_empty() {
        return 0.0;
    }
    let mut sum = c[0];
    if c.len() >= 2 {
        sum += c[1] * x;
    }
    let mut t_prev = 1.0;
    let mut t_curr = x;
    for k in 2..c.len() {
        let t_next = 2.0 * x * t_curr - t_prev;
        sum += c[k] * t_next;
        t_prev = t_curr;
        t_curr = t_next;
    }
    sum
}

// #[inline(always)]
// pub fn chebyshev_eval(c: &[f64], x: f64) -> f64 {
//     let n = c.len();
//     if n == 0 {
//         return 0.0;
//     }
//     if n == 1 {
//         return c[0];
//     }

//     let two_x = 2.0 * x;

//     // Clenshaw recurrence:
//     // b_{k} = 2x b_{k+1} - b_{k+2} + c_k
//     // p(x) = b_0 - x b_1
//     let mut b_k1 = 0.0; // b_{k+1}
//     let mut b_k2 = 0.0; // b_{k+2}

//     for &ck in c[1..].iter().rev() {
//         let b_k = two_x.mul_add(b_k1, ck) - b_k2;
//         b_k2 = b_k1;
//         b_k1 = b_k;
//     }

//     // incorporate c0 and finalize
//     let b0 = two_x.mul_add(b_k1, c[0]) - b_k2;
//     b0 - x * b_k1
// }
