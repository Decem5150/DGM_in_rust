use std::collections::HashMap;
use ndarray::{Array, Ix1};
use ndarray::{Ix2, Ix3};
pub use super::gauss_point::GaussPoints;
#[derive(Debug)]
pub struct DubinerBasis {
    pub dof: usize,
    pub phis_cell_gps: Array<f64, Ix2>,
    pub phis_edge_gps: Array<f64, Ix3>,
    pub derivatives: Array<HashMap<(usize, usize), f64>, Ix2>,
    pub mass_mat_diag: Array<f64, Ix1>,
}
impl DubinerBasis {
    pub fn new (dof: usize, gauss_points: &GaussPoints) -> DubinerBasis {
        let gp_number = gauss_points.cell_gp_number;
        let edge_gp_number = gauss_points.edge_gp_number;
        let phis_cell_gps = Array::zeros((gp_number, dof));
        let phis_edge_gps = Array::zeros((3, edge_gp_number, dof));
        let derivatives: Array<HashMap<(usize, usize), f64>, Ix2> = Array::from_elem((gp_number, dof), HashMap::new());
        let mass_mat_diag = Array::zeros(dof);
        let mut basis = DubinerBasis {
            dof,
            phis_cell_gps,
            phis_edge_gps,
            derivatives,
            mass_mat_diag,
        };
        basis.compute_phi_cell(gauss_points);
        basis.compute_phi_edge(gauss_points);
        basis.compute_derivatives(gauss_points);
        basis.compute_mass_mat(gauss_points);
        //basis.debug_basis(gauss_points);
        basis
    }
    fn debug_basis(&self, gauss_points: &GaussPoints) {
        let mut mass_mat: Array<f64, Ix2> = Array::zeros((self.dof, self.dof));
        let nbasis = self.dof;
        let ngp = gauss_points.cell_gp_number;
        for ibasis in 0..nbasis {
            for jbasis in 0..nbasis{
                mass_mat[[ibasis, jbasis]] = {
                    let mut sum = 0.0;
                    for igp in 0..ngp {
                        sum += gauss_points.cell_weights[igp] * self.phis_cell_gps[[igp, ibasis]] * self.phis_cell_gps[[igp, jbasis]];
                    }
                    0.5 * sum
                };
            }
        }
        dbg!(mass_mat);
    }
    fn jacobi_poly(n: i32, alpha: f64, beta: f64, x: f64) -> f64 {
        match n {
            0 => 1.0,
            1 => (alpha + 1.0) + (alpha + beta + 2.0) * (x - 1.0) / 2.0,
            _ => {
                let n = n as f64;
                let a1 = 2.0 * n * (n + alpha + beta) * (2.0 * n + alpha + beta - 2.0);
                let a2 = 2.0 * n + alpha + beta - 1.0;
                let a3 = (2.0 * n + alpha + beta) * (2.0 * n + alpha + beta - 2.0) * x + alpha.powi(2) - beta.powi(2);
                let a4 = 2.0 * (n + alpha - 1.0) * (n + beta - 1.0) * (2.0 * n + alpha + beta);
    
                let n = n as i32;
                let pn_1 = Self::jacobi_poly(n - 1, alpha, beta, x);
                let pn_2 = Self::jacobi_poly(n - 2, alpha, beta, x);
    
                (a2 * a3 * pn_1 - a4 * pn_2) / a1
            }
        }
    }
    fn dubiner_basis(l: i32, k: i32, xi: f64, eta: f64) -> f64 {
        let legendre = Self::jacobi_poly(k, 0.0, 0.0, 2.0 * xi / (1.0 - eta) - 1.0);
        let jacobi = Self::jacobi_poly(l, 2.0 * k as f64 + 1.0, 0.0, 2.0 * eta - 1.0);
        2.0_f64.powi(k) * legendre * (1.0 - eta).powi(k) * jacobi
    }
    fn compute_phi_cell(&mut self, gauss_points: &GaussPoints) {
        let gp_number = gauss_points.cell_gp_number;
        let cell_points = &gauss_points.cell_points;
        for i in 0..gp_number {
            let xi = cell_points[[i, 0]];
            let eta = cell_points[[i, 1]];
            self.phis_cell_gps[[i, 0]] = Self::dubiner_basis(0, 0, xi, eta);
            self.phis_cell_gps[[i, 1]] = Self::dubiner_basis(0, 1, xi, eta);
            self.phis_cell_gps[[i, 2]] = Self::dubiner_basis(1, 0, xi, eta);
            /* 
            self.phis_cell_gps[[i, 3]] = Self::dubiner_basis(0, 2, xi, eta);
            self.phis_cell_gps[[i, 4]] = Self::dubiner_basis(1, 1, xi, eta);
            self.phis_cell_gps[[i, 5]] = Self::dubiner_basis(2, 0, xi, eta);
            */
        }
    }
    fn compute_phi_edge(&mut self, gauss_points: &GaussPoints) {
        let gp_number = gauss_points.edge_gp_number;
        let edge_points = &gauss_points.edge_points;
        for iedge in 0..3 {
            for igp in 0..gp_number {
                let xi = edge_points[[iedge, igp, 0]];
                let eta = edge_points[[iedge, igp, 1]];
                self.phis_edge_gps[[iedge, igp, 0]] = Self::dubiner_basis(0, 0, xi, eta);
                self.phis_edge_gps[[iedge, igp, 1]] = Self::dubiner_basis(0, 1, xi, eta);
                self.phis_edge_gps[[iedge, igp, 2]] = Self::dubiner_basis(1, 0, xi, eta);
                /* 
                self.phis_edge_gps[[iedge, igp, 3]] = Self::dubiner_basis(0, 2, xi, eta);
                self.phis_edge_gps[[iedge, igp, 4]] = Self::dubiner_basis(1, 1, xi, eta);
                self.phis_edge_gps[[iedge, igp, 5]] = Self::dubiner_basis(2, 0, xi, eta);
                */
            }
        }
    }
    fn compute_derivatives(&mut self, gauss_points: &GaussPoints) {
        let ngp = gauss_points.cell_gp_number;
        let gauss_points = &gauss_points.cell_points;
        for i in 0..ngp {
            let xi = gauss_points[[i, 0]];
            let eta = gauss_points[[i, 1]];
            self.derivatives[[i, 0]].insert((1, 0), 0.0);
            self.derivatives[[i, 0]].insert((0, 1), 0.0);
            self.derivatives[[i, 0]].insert((2, 0), 0.0);
            self.derivatives[[i, 0]].insert((1, 1), 0.0);
            self.derivatives[[i, 0]].insert((0, 2), 0.0);
            self.derivatives[[i, 1]].insert((1, 0), 4.0);
            self.derivatives[[i, 1]].insert((0, 1), 2.0);
            self.derivatives[[i, 1]].insert((2, 0), 0.0);
            self.derivatives[[i, 1]].insert((1, 1), 0.0);
            self.derivatives[[i, 1]].insert((0, 2), 0.0);
            self.derivatives[[i, 2]].insert((1, 0), 0.0);
            self.derivatives[[i, 2]].insert((0, 1), 3.0);
            self.derivatives[[i, 2]].insert((2, 0), 0.0);
            self.derivatives[[i, 2]].insert((1, 1), 0.0);
            self.derivatives[[i, 2]].insert((0, 2), 0.0);
            /* 
            self.derivatives[[i, 3]].insert((1, 0), 48.0 * xi + 24.0 * eta - 24.0);
            self.derivatives[[i, 3]].insert((0, 1), 24.0 * xi + 8.0 * eta - 8.0);
            self.derivatives[[i, 3]].insert((2, 0), 48.0);
            self.derivatives[[i, 3]].insert((1, 1), 24.0);
            self.derivatives[[i, 3]].insert((0, 2), 8.0);
            self.derivatives[[i, 4]].insert((1, 0), 20.0 * eta - 4.0);
            self.derivatives[[i, 4]].insert((0, 1), 20.0 * xi + 20.0 * eta - 12.0);
            self.derivatives[[i, 4]].insert((2, 0), 0.0);
            self.derivatives[[i, 4]].insert((1, 1), 20.0);
            self.derivatives[[i, 4]].insert((0, 2), 20.0);
            self.derivatives[[i, 5]].insert((1, 0), 0.0);
            self.derivatives[[i, 5]].insert((0, 1), 20.0 * eta - 8.0);
            self.derivatives[[i, 5]].insert((2, 0), 0.0);
            self.derivatives[[i, 5]].insert((1, 1), 0.0);
            self.derivatives[[i, 5]].insert((0, 2), 20.0);
            */
        }
    }
    fn compute_mass_mat(&mut self, gauss_points: &GaussPoints) {
        //self.mass_mat_diag.iter_mut().for_each (|mass| *mass = 0.0);
        let nbasis = self.dof;
        let ngp = gauss_points.cell_gp_number;
        for ibasis in 0..nbasis {
            for igp in 0..ngp {
                self.mass_mat_diag[ibasis] += gauss_points.cell_weights[igp] * self.phis_cell_gps[[igp, ibasis]].powf(2.0);
            }
        }
    }
}