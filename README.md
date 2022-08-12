# rust_scf
Implementation of the SCF algorithm for two-electron system using sto-3g basis set

Follows the original fortran implementation from Szabo & Ostlund, with some idiomatic changes

TODO:
* Expand to > 2 electron systems
* Extend the basis function library to use Gaussian Lobe fits from Whitten (https://doi.org/10.1063/1.1726470)
* Refactor two electron computation and storage if this becomes an issue for larger systems
