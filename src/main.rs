use clap::Parser;
use ndarray::{Array, Array2, Array4};
use ndarray_linalg::{UPLO, Eigh};
use rust_scf::integrals::{kinetic_integral, nuclear_attraction_integral, overlap_integral, two_electron_integral};
use std::f64::consts::PI;

#[derive(Parser, Debug)]
struct Args {
    #[clap(long, value_parser, help = "number of electrons", default_value = "3")]
    n: u32,
    #[clap(
        long,
        value_parser,
        help = "atomic separation",
        default_value = "1.463200"
    )]
    r: f64,
    #[clap(
        long,
        value_parser,
        help = "atom A nuclear charge",
        default_value = "2"
    )]
    za: u32,
    #[clap(
        long,
        value_parser,
        help = "atom B nuclear charge",
        default_value = "1"
    )]
    zb: u32,
    #[clap(
        long,
        value_parser,
        help = "basis function 1 Slater exponent",
        default_value = "2.09250"
    )]
    zeta1: f64,
    #[clap(
        long,
        value_parser,
        help = "basis function 2 Slater exponent",
        default_value = "1.2400"
    )]
    zeta2: f64,
}

fn hf_calc(n: u32, r: f64, zeta_1: f64, zeta_2: f64, z_a: u32, z_b: u32) {
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
    d1 = d1
        .indexed_iter_mut()
        .map(|(idx, x)| *x * &a1[idx] * (2. / PI).powf(0.75))
        .collect();

    // basis functions on atom 2
    let mut a2 = expons.clone();
    a2 = a2.iter_mut().map(|x| *x * zeta_1.powf(2.)).collect();
    let mut d2 = coeffs.clone();
    d2 = d2
        .indexed_iter_mut()
        .map(|(idx, x)| *x * &a2[idx] * (2. / PI).powf(0.75))
        .collect();
    println!("{:?}", d2);

    // let all_coeffs = Array::from_shape_vec((3,1), )

    // now we have the normalized basis functions.  Compute the standard integrals
    /* Integrals are the following
    OVERLAP
    S12
    KINETIC ENERGY
    T11, T12, T22
    ELECTRON-NUCLEAR ATTRACTION
    V11A, V11B, V22A, V22B, V12A, V12B
    TWO ELECTRON INTEGRALS
    V1111, V2111, V2121, V2211, V2221, V2222
    */
    // let (s12, t11, t12, t22, v11a, v11b, v22a, v22b, v12a, v12b, v1111, v2111, v2121, v2211, v2221, v2222) = ()
    let mut s12 = 0.;
    let mut t11 = 0.;
    let mut t12 = 0.;
    let mut t22 = 0.;
    let mut v11a = 0.;
    let mut v11b = 0.;
    let mut v22a = 0.;
    let mut v22b = 0.;
    let mut v12a = 0.;
    let mut v12b = 0.;

    // loop over the number of electrons
    for i  in 0..n {
        for j in 0..n {
			let i: usize = i.try_into().unwrap();
			let j: usize = j.try_into().unwrap();
			let r_ap = a2[j] * r / (a1[i] + a2[j]);
			let r_bp = r - r_ap;
            s12 += overlap_integral(a1[i], a2[j], r) * d1[i] * d2[j];
            t11 += kinetic_integral(a1[i], a1[j], 0.) * d1[i] * d1[j];
            t12 += kinetic_integral(a1[i], a2[j], r) * d1[i] * d2[j];
            t22 += kinetic_integral(a1[i], a2[j], 0.) * d2[i] * d2[j];
            v11a += nuclear_attraction_integral(a1[i], a1[j], 0., 0., z_a) * d1[i] * d1[j];
            v11b += nuclear_attraction_integral(a1[i], a1[j], 0., 0., z_b)* d1[i] * d1[j];
            v22a += nuclear_attraction_integral(a2[i], a2[j], 0., 0., z_a)* d2[i] * d2[j];
            v22b += nuclear_attraction_integral(a2[i], a2[j], 0., 0., z_b)* d2[i] * d2[j];
            v12a += nuclear_attraction_integral(a1[i], a2[j], r, r_ap, z_a)* d1[i] * d2[j];
            v12b += nuclear_attraction_integral(a1[i], a2[j], r, r_bp, z_b)* d1[i] * d2[j];

           
        }
    }
	let mut v1111 = 0.;
	let mut v2111 = 0.;
	let mut v2121 = 0.;
	let mut v2211 = 0.;
	let mut v2221 = 0.;
	let mut v2222 = 0.;

	for i in 0..n {
		for j in 0..n {
			for k in 0..n {
				for l in 0..n {
					let i: usize = i.try_into().unwrap();
					let j: usize = j.try_into().unwrap();
					let k: usize = k.try_into().unwrap();
					let l: usize = l.try_into().unwrap();
					let r_ap = a2[j] * r / (a2[i] + a1[j]);
					let r_bp = r - r_ap;
					let r_aq = a2[k] * r / (a2[k] + a1[l]);
					let r_bq = r - r_aq;
					let r_pq  = r_ap - r_aq;
					v1111 += two_electron_integral(a1[i], a1[j], a1[k], a1[l], 0., 0., 0.) * d1[i] * d1[j] * d1[k] * d1[l];
					v2111 += two_electron_integral(a2[i], a1[j], a1[k], a1[l], r, 0., r_ap) * d2[i] * d1[j] * d1[k] * d1[l];
					v2121 += two_electron_integral(a2[i], a1[j], a2[k], a1[l], r, r, r_pq) * d2[i] * d1[j] * d2[k] * d1[l];
					v2211 += two_electron_integral(a2[i], a2[j], a1[k], a1[l], 0., 0., r) * d2[i] * d2[j] * d1[k] * d1[l];
					v2221 += two_electron_integral(a2[i], a2[j], a2[k], a1[l], 0., r, r_bq) * d2[i] * d2[j] * d2[k] * d1[l];
					v2222 += two_electron_integral(a2[i], a2[j], a2[k], a2[l], 0., 0., 0.) * d2[i] * d2[j] * d2[k] * d2[l];


				}
			}
		}
	}
	
	let (h, s, x, xt, tt) = collect(s12, t11, t12, t22, v11a, v11b, v22a, v22b, v12a, v12b, v1111, v2111, v2121, v2211, v2221, v2222);


	scf(r, z_a, z_b, h, s, x, xt, tt);
    
}

fn collect(s12: f64, t11: f64, t12: f64, t22: f64, v11a: f64, v11b: f64, v22a: f64, v22b: f64, v12a: f64, v12b: f64, v1111: f64, v2111: f64, v2121: f64, v2211: f64, v2221: f64, v2222: f64) -> (Array2<f64>, Array2<f64>, Array2<f64>, Array2<f64>, Array4<f64>){
	// put the values into the necessary arrays
	let mut h = Array::zeros((2,2));
	let mut s = Array::eye(2);
	let mut x = Array::zeros((2,2));
	// store two-electron integrals
	let mut t = Array::zeros((2,2,2,2));
	println!("{},{},{}", t11, v11a, v11b);
	h[[0,0]] = t11 + v11a + v11b;
	h[[0,1]] = t12 + v12a + v12b;
	h[[1,0]] = h[[0,1]];
	h[[1,1]] = t22 + v22a + v22b;

	s[[0,1]] = s12;
	s[[1,0]] = s12;

	t[[0,0,0,0]] = v1111;

	t[[1,0,0,0]] = v2111;
	t[[0,1,0,0]] = v2111;
	t[[0,0,1,0]] = v2111;
	t[[0,0,0,1]] = v2111;
	
	t[[1,0,1,0]] = v2121;
	t[[0,1,1,0]] = v2121;
	t[[1,0,0,1]] = v2121;
	t[[0,1,0,1]] = v2121;

	t[[1,1,0,0]] = v2211;
	t[[0,0,1,1]] = v2211;

	t[[1,1,1,0]] = v2221;
	t[[1,1,0,1]] = v2221;
	t[[1,0,1,1]] = v2221;
	t[[0,1,1,1]] = v2221;

	t[[1,1,1,1]] = v2222;


	// canonical normalization
	x[[0,0]] = 1. / (2. * (1. + s12)).sqrt();
	x[[1,0]] = x[[0,0]];
	x[[0,1]] = 1. / (2. * (1. - s12)).sqrt();
	x[[1,1]] = -x[[0,1]];

	let xt = x.clone().reversed_axes();

	println!("h: {:?}\ns: {:?}\nt: {:?}\n\n", h, s, t);

	return (h, s, x, xt, t)

}

fn compute_g_matrix(p: &Array2<f64>, tt: &Array4<f64>) -> Array2<f64> {
	let mut g = Array::zeros((2,2));
	for i in 0..2 {
		for j in 0..2 {
			for k in 0..2 {
				for l in 0..2 {
					// coulomb - exchange * P matrix
					g[[i,j]] = p[[k,l]] * (tt[[i,j,k,l]] - 0.5 * tt[[i,j,k,l]])
				}
			}
		}
	}
	g
}

fn scf(r: f64, z_a: u32, z_b: u32, h: Array2<f64>, s: Array2<f64>, x: Array2<f64>, xt: Array2<f64>, tt: Array4<f64>) {

	let criteria = 1e-3;
	let maxiter = 50;
	let mut iter = 0;


	// initial guess at the density matrix, i.e. use core hamiltonian term only
	let mut p = Array::zeros((2,2));
	let mut f = h.clone();
	let mut e_elec = 0.;
	while iter < maxiter {
		iter += 1;
		// generate the G matrix from P and the two electron integrals
		let mut g = compute_g_matrix(&p, &tt);

		// generate fock matrix 
		
		f.indexed_iter_mut().for_each(|((i,j), x)| *x += g[[i,j]]);
		

		// compute electronic energy = $E_0 = 0.5 * \Sum_{u,v} P_{\mu,\nu} (H_{\mu, \nu} + F_{\mu,\nu}$ 
		e_elec = 0.;
		for i in 0..2 {
			for j in 0..2 {
				e_elec += 0.5 * p[[i,j]] * (h[[i,j]] + f[[i,j]]);
			}
		}
		println!("Iteration: {}, \nEnergy: {}\n\n", &iter, &e_elec);

		// Transform the fock matrix 
		g = f.clone().dot(&x);
		let f_prime = xt.clone().dot(&g);

		// eigenvalues e (orbital energies), eigenvectors c_prime (transformed coefficient matrix)
		let (e, c_prime) = f_prime.eigh(UPLO::Lower).unwrap();
		// transform eigenvectors 
		let c = x.clone().dot(&c_prime);

		// compute new density matrix P
		let p_prev = p.clone();
		// eqn 3.145
		p = Array::zeros((2,2));
		for i in 0..2 {
			for j in 0..2 {
				// use zero index since only one occupied orbital here
				p[[i,j]] += 2. * c[[i,0]] * c[[j,0]]
			}
		}

		// compute delta in P matrix
		let mut  delta = 0.;
		for i in 0..2 {
			for j in 0..2 {
				delta += (p[[i,j]] - p_prev[[i,j]]).powf(2.)
			}
		}
		delta = (delta / 4.).sqrt();

		println!("Iter: {}, \nDelta: {}", iter, delta);
		if delta < criteria || iter > maxiter {
			println!("SCF finished!");
			break
		}
	}

	// compute total energy & expectation values from the C P and F matrices
	let e_total = e_elec + z_a as f64 * z_b as f64 / r;

	println!("Converged energies:\n");
	println!("e_total: {}", e_total);
	println!("e_elec: {}", e_elec);


	// print mulliken populations
	println!("Mulliken populations: {:?}", p.dot(&s))
}


fn main() {
    let args = Args::parse();

    // collect the required integrals
    hf_calc(args.n, args.r, args.zeta1, args.zeta2, args.za, args.zb)

    // Run SCF iterations until convergence
}
