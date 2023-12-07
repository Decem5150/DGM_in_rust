use std::rc::Rc;
use ndarray::Array;
use ndarray::{Ix2, Ix3};
use sprs::{TriMat, CsMat};
pub use super::gauss_point::GaussPoints;
pub struct DubinerBasis {
    pub dof: usize,
    pub phis_cell_gps: Array<f64, Ix2>,
    pub phis_edge_gps: Array<f64, Ix3>,
    pub derivatives: Array<CsMat<f64>, Ix2>,
    pub gauss_points: Rc<GaussPoints>,
}
impl DubinerBasis {
    fn new (dof: usize, gauss_points: Rc<GaussPoints>) -> DubinerBasis {
        let gp_number = gauss_points.cell_gp_number;
        let edge_gp_number = gauss_points.edge_gp_number;
        let phis_cell_gps = Array::zeros((gp_number, dof));
        let phis_edge_gps = Array::zeros((3, edge_gp_number, dof));
        let derivatives: Array<CsMat<f64>, _> = Array::default((gp_number, dof));
        let mut basis = DubinerBasis {
            dof,
            phis_cell_gps,
            phis_edge_gps,
            derivatives,
            gauss_points,
        };
        basis.compute_phi_cell();
        basis.compute_phi_edge();
        basis.compute_derivatives();
        basis
    }
    fn compute_phi_cell(&mut self) {
        let gp_number = self.gauss_points.cell_gp_number;
        let cell_points = &self.gauss_points.cell_points;
        for i in 0..gp_number {
            let xi = cell_points[[i, 0]];
            let eta = cell_points[[i, 1]];
            self.phis_cell_gps[[i, 0]] = 1.0;
            self.phis_cell_gps[[i, 1]] = xi;
            self.phis_cell_gps[[i, 2]] = eta;
            self.phis_cell_gps[[i, 3]] = xi * (2.0 * xi - 1.0);
            self.phis_cell_gps[[i, 4]] = 4.0 * xi * eta;
            self.phis_cell_gps[[i, 5]] = eta * (2.0 * eta - 1.0);
        }
    }
    fn compute_phi_edge(&mut self) {
        let gp_number = self.gauss_points.edge_gp_number;
        let edge_points = &self.gauss_points.edge_points;
        for iedge in 0..3 {
            for igp in 0..gp_number {
                let xi = edge_points[[iedge, igp, 0]];
                let eta = edge_points[[iedge, igp, 1]];
                self.phis_edge_gps[[iedge, igp, 0]] = 1.0;
                self.phis_edge_gps[[iedge, igp, 1]] = xi;
                self.phis_edge_gps[[iedge, igp, 2]] = eta;
                self.phis_edge_gps[[iedge, igp, 3]] = xi * (2.0 * xi - 1.0);
                self.phis_edge_gps[[iedge, igp, 4]] = 4.0 * xi * eta;
                self.phis_edge_gps[[iedge, igp, 5]] = eta * (2.0 * eta - 1.0);
            }
        }
    }
    fn compute_derivatives(&mut self) {
        let ngp = self.gauss_points.cell_gp_number;
        let nbasis = self.dof;
        let gauss_points = &self.gauss_points.cell_points;
        for i in 0..ngp {
            let xi = gauss_points[[i, 0]];
            let eta = gauss_points[[i, 1]];
            // 1st basis
            let mut derivatives = TriMat::new((3, 3));
            self.derivatives[[i, 0]] = derivatives.to_csr();
            // 2nd basis
            let mut derivatives = TriMat::new((3, 3));
            derivatives.add_triplet(1, 0, 1.0);
            self.derivatives[[i, 1]] = derivatives.to_csr();
            // 3rd basis
            let mut derivatives = TriMat::new((3, 3));
            derivatives.add_triplet(0, 1, 1.0);
            self.derivatives[[i, 2]] = derivatives.to_csr();
            // 4th basis 
            let mut derivatives = TriMat::new((3, 3));
            derivatives.add_triplet(1, 0, 4.0 * xi - 1.0);
            derivatives.add_triplet(2, 0, 4.0);
            self.derivatives[[i, 3]] = derivatives.to_csr();
            // 5th basis
            let mut derivatives = TriMat::new((3, 3));
            derivatives.add_triplet(1, 0, 4.0 * eta);
            derivatives.add_triplet(0, 1, 4.0 * xi);
            derivatives.add_triplet(1, 1, 4.0);
            self.derivatives[[i, 4]] = derivatives.to_csr();
            // 6th basis
            let mut derivatives = TriMat::new((3, 3));
            derivatives.add_triplet(0, 1, 4.0 * eta - 1.0);
            derivatives.add_triplet(0, 2, 4.0);
            self.derivatives[[i, 5]] = derivatives.to_csr();
        }
    }
}