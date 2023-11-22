use ndarray::Array;
use ndarray::{Ix2, Ix3};
pub use super::gauss_point::GaussPoints;
pub struct DubinerBasis<'a> {
    pub dof: usize,
    pub phis_cell_gps: Array<f64, Ix2>,
    pub phis_edge_gps: Array<f64, Ix3>,
    pub dphis_dxi: Array<f64, Ix2>,
    pub dphis_deta: Array<f64, Ix2>,
    pub gauss_points: &'a GaussPoints,
}
impl DubinerBasis<'_> {
    fn new (dof: usize, gauss_points: &GaussPoints) -> DubinerBasis {
        let gp_number = gauss_points.cell_gp_number;
        let edge_gp_number = gauss_points.edge_gp_number;
        let phis_cell_gps = Array::zeros((gp_number, dof));
        let phis_edge_gps = Array::zeros((3, edge_gp_number, dof));
        let dphis_dxi = Array::zeros((gp_number, dof));
        let dphis_deta = Array::zeros((gp_number, dof));
        let mut basis = DubinerBasis {
            dof,
            phis_cell_gps,
            phis_edge_gps,
            dphis_dxi,
            dphis_deta,
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
            self.phis_cell_gps[[i, 3]] = xi * xi - 1.0 / 3.0;
            self.phis_cell_gps[[i, 4]] = xi * eta;
            self.phis_cell_gps[[i, 5]] = eta * eta - 1.0 / 3.0;
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
        let gp_number = self.gauss_points.cell_gp_number;
        let gauss_points = &self.gauss_points.cell_points;
        for i in 0..gp_number {
            let xi = gauss_points[[i, 0]];
            let eta = gauss_points[[i, 1]];
            self.dphis_dxi[[i, 0]] = 0.0;
            self.dphis_dxi[[i, 1]] = 1.0;
            self.dphis_dxi[[i, 2]] = 0.0;
            self.dphis_dxi[[i, 3]] = 2.0 * xi;
            self.dphis_dxi[[i, 4]] = eta;
            self.dphis_dxi[[i, 5]] = 0.0;
            self.dphis_deta[[i, 0]] = 0.0;
            self.dphis_deta[[i, 1]] = 0.0;
            self.dphis_deta[[i, 2]] = 1.0;
            self.dphis_deta[[i, 3]] = 0.0;
            self.dphis_deta[[i, 4]] = xi;
            self.dphis_deta[[i, 5]] = 2.0 * eta;
        }
    }
}