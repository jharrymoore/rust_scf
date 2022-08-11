use clap::Parser;
use ndarray::Array;
use rust_scf::integrals::{nuclear_attraction_integral, overlap_integral};
use std::f64::consts::PI;

#[derive(Parser, Debug)]
struct Args {
    #[clap(long, value_parser, help = "number of electrons", default_value = "3")]
    n: u32,
    #[clap(long, value_parser, help = "atomic separation", default_value="1.463200")]
    r: f64,
    #[clap(long, value_parser, help = "atom A nuclear charge", default_value="2")]
    za: u32,
    #[clap(long, value_parser, help = "atom B nuclear charge", default_value="1")]
    zb: u32,
    #[clap(long, value_parser, help = "basis function 1 Slater exponent", default_value="2.09250")]
    zeta1: f64,
    #[clap(long, value_parser, help = "basis function 2 Slater exponent", default_value="1.2400")]
    zeta2: f64,
}

fn compute_integrals(n: u32, r: f64, zeta_1: f64, zeta_2: f64, z_a: u32, z_b: u32) {
    // compute and store one- and two-electron integrals for the full basis set - assume no memory bottlebeck
    // let coeffs = Array::from_shape_vec(
    //     (3, 3),
    //     vec![
    //         1.0, 0.0, 0.0, 0.678914, 0.430129, 0.0, 0.444635, 0.535328, 0.154329,
    //     ],
    // )
    // .unwrap();
    // let alpha = Array::from_shape_vec(
    //     (3, 3),
    //     vec![
    //         0.270950, 0.0, 0.0, 0.1151623, 0.851819, 0.0, 0.109818, 0.405771, 2.22766,
    //     ],
    // )
    // .unwrap();

	// implement only sto-3g linear fits for now
	// these are the unnormalized coeffs
    let expons = Array::from_vec(vec![0.109818, 0.405771, 2.22766]);

    let coeffs = Array::from_vec(vec![0.444635, 0.535328, 0.154329]);

	// basis functions on atom 1
	let mut a1 = expons.clone();
	a1 = a1.iter_mut().map(|x| *x * zeta_1.powf(2.)).collect();
	let mut d1 = coeffs.clone();
	d1 = d1.indexed_iter_mut().map(|(idx, x)| *x * &a1[idx] * (2. / PI).powf(0.75)).collect();

	// basis functions on atom 2
	let mut a2 = expons.clone();
	a2 = a2.iter_mut().map(|x| *x * zeta_1.powf(2.)).collect();
	let mut d2 = coeffs.clone();
	d2 = d2.indexed_iter_mut().map(|(idx, x)| *x * &a2[idx] * (2. / PI).powf(0.75)).collect();
	println!("{:?}", d2);

	// let all_coeffs = Array::from_shape_vec((3,1), )

	// now we have the normalized basis functions.  Compute the standard integrals


	



	// normalize the coeffs and the 

}

fn main() {
    let args = Args::parse();

    // collect the required integrals
    compute_integrals(args.n, args.r, args.zeta1, args.zeta2, args.za, args.zb)

    // Run SCF iterations until convergence
}
