#[allow(dead_code)]
pub mod integrals {
    use std::f64::consts::PI;

    use ndarray_linalg::Scalar;
    use statrs::function::erf::erf;
    pub fn overlap_integral(a: f64, b: f64, r: f64) -> f64 {
        (PI / (a + b)).powf(1.5) * ((-a * b / (a + b)) * r.powf(2.)).exp()
    }

    pub fn kinetic_integral(a: f64, b: f64, r_ab: f64) -> f64 {
        (a * b / (a + b)) * (3. - 2. * a * b / (a + b) * r_ab.powf(2.)) * (PI / (a + b)).powf(1.5) * (-a * b / (a + b) * r_ab.pow(2.)).exp()
    }

    pub fn nuclear_attraction_integral(a: f64, b: f64, r_ab: f64, r_pc: f64, z_c: u32) -> f64 {
        // Computes the nuclear attraction integral for <A|-Z/r|B> two basis functions A and B at separation r
        -2. * PI / (a + b)
            * z_c as f64
            * (-a * b / (a + b) * r_ab.powf(2.)).exp()
            * f_0((a + b) * r_pc.powf(2.))
    }

    pub fn f_0(t: f64) -> f64 {
        if t < 1e-6 {
            1. - t / 3.
        } else {
            0.5 * (PI / t).powf(0.5) * erf(t.powf(0.5))
        }
    }

    pub fn two_electron_integral(
        a: f64,
        b: f64,
        c: f64,
        d: f64,
        r_ij: f64,
        r_kl: f64,
        r_pq: f64,
    ) -> f64 {
        2. * PI.powf(2.5) / ((a + b) * (c * d) * (a + b + c + d).sqrt())
            * (-a * b / (a + b) * r_ij.powf(2.) - c * d / (c + d) * r_kl.powf(2.)).exp()
            * f_0((a + b) * (c + d) / (a + b + c + d) * r_pq)
    }
}
