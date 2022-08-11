use clap::Parser;
use ndarray::Array;
use rust_scf::integrals::{nuclear_attraction_integral, overlap_integral};

#[derive(Parser, Debug)]
struct Args {
    #[clap(long, value_parser, help = "number of electrons")]
    n: u32,
    #[clap(long, value_parser, help = "atomic separation")]
    r: f64,
    #[clap(long, value_parser, help = "atom A nuclear charge")]
    za: u32,
    #[clap(long, value_parser, help = "atom B nuclear charge")]
    zb: u32,
    #[clap(long, value_parser, help = "basis function 1 Slater exponent")]
    zeta1: f64,
    #[clap(long, value_parser, help = "basis function 2 Slater exponent")]
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

    let coeffs = Array::from_vec(vec![0.444635, 0.535328, 0.154329]);

    let alpha = Array::from_vec(vec![0.109818, 0.405771, 2.22766]);

	// normalize the coeffs
}

fn main() {
    let args = Args::parse();

    // collect the required integrals
    compute_integrals(args.n, args.r, args.zeta1, args.zeta2, args.za, args.zb)

    // Run SCF iterations until convergence
}
