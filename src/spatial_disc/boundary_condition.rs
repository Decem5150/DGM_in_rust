use ndarray::{Ix1, Ix2, Ix3};
use ndarray::{Array, Axis};
use ndarray::s;
use crate::basis_function::DubinerBasis;
use crate::mesh::{BoundaryEdge, BoundaryType, Patch, Mesh};
use super::local_characteristics;
pub struct BoundaryCondition {
    pub ngp: usize,
    pub neq: usize,
    pub nbasis: usize,
    pub hcr: f64,
}
impl BoundaryCondition {
    pub fn apply_bc(&self, residuals: &mut Array<f64, Ix3>, solutions: &Array<f64, Ix3>, mesh: &Mesh, basis: &DubinerBasis) {
        for patch in mesh.patches.iter() {
            match patch.boundary_type {
                BoundaryType::Wall => self.wall(residuals, solutions, mesh, patch, basis),
                BoundaryType::FarField => self.far_field(residuals, solutions, mesh, patch, basis),
                _ => panic!("Boundary type not implemented"),
            }
        }
    }
    fn wall(&self, residuals: &mut Array<f64, Ix3>, solutions: &Array<f64, Ix3>, mesh: &Mesh, patch: &Patch, basis: &DubinerBasis) {
        for &iedge in patch.iedges.iter() {
            let ielem = mesh.boundary_edges[iedge].ielement;
            let nx = mesh.boundary_edges[iedge].normal[0];
            let ny = mesh.boundary_edges[iedge].normal[1];
            let left_values_gps: Array<f64, Ix2> = self.compute_boundary_edges_values(&mesh.boundary_edges[iedge], solutions, basis);
            for igp in 0..self.ngp {
                for ibasis in 0..self.nbasis {
                    let pressure = (self.hcr - 1.0) * (left_values_gps[[igp, 3]] - 0.5 * (left_values_gps[[igp, 1]] * left_values_gps[[igp, 1]] + left_values_gps[[igp, 2]] * left_values_gps[[igp, 2]]) / left_values_gps[[igp, 0]]);
                    residuals[[ielem, 1, ibasis]] -= pressure * basis.phis_edge_gps[[ielem, igp, ibasis]] * nx * 0.5 * mesh.boundary_edges[iedge].jacob_det;
                    residuals[[ielem, 2, ibasis]] -= pressure * basis.phis_edge_gps[[ielem, igp, ibasis]] * ny * 0.5 * mesh.boundary_edges[iedge].jacob_det;
                }
            }
        }
    }
    fn far_field(&self, residuals: &mut Array<f64, Ix3>, solutions: &Array<f64, Ix3>, mesh: &Mesh, patch: &Patch, basis: &DubinerBasis) {
        let boundary_quantity = patch.boundary_quantity.as_ref().unwrap();
        let boundary_consvar = Array::from_vec(vec![
            boundary_quantity.rho, 
            boundary_quantity.rho * boundary_quantity.u, 
            boundary_quantity.rho * boundary_quantity.v, 
            boundary_quantity.p / (self.hcr - 1.0) + 0.5 * boundary_quantity.rho * (boundary_quantity.u.powi(2) + boundary_quantity.v.powi(2))
            ]);
        for &iedge in patch.iedges.iter() {
            let ielem = mesh.boundary_edges[iedge].ielement;
            let nx = mesh.boundary_edges[iedge].normal[0];
            let ny = mesh.boundary_edges[iedge].normal[1];
            let left_values_gps: Array<f64, Ix2> = self.compute_boundary_edges_values(&mesh.boundary_edges[iedge], solutions, basis);
            for igp in 0..self.ngp {
                for ibasis in 0..self.nbasis {
                    let u = left_values_gps[[igp, 1]] / left_values_gps[[igp, 0]];
                    let v = left_values_gps[[igp, 2]] / left_values_gps[[igp, 0]];
                    let p = (self.hcr - 1.0) * (left_values_gps[[igp, 3]] - 0.5 * (left_values_gps[[igp, 1]] * left_values_gps[[igp, 1]] + left_values_gps[[igp, 2]] * left_values_gps[[igp, 2]]) / left_values_gps[[igp, 0]]);
                    let c = (self.hcr * p / left_values_gps[[igp, 0]]).sqrt();
                    let vn = u * nx + v * ny;
                    let (left_elem_lmatrix, left_elem_rmatrix) = local_characteristics::compute_eigenmatrix(left_values_gps.slice(s![igp, ..]), nx, ny, self.hcr);
                    let left_elem_eigenvalues = [vn - c, vn, vn, vn + c];
                    let alpha: Array<f64, Ix2> = left_elem_lmatrix.dot(&left_values_gps.slice(s![igp, ..]).insert_axis(Axis(1)));
                    let beta: Array<f64, Ix2> = left_elem_lmatrix.dot(&boundary_consvar.clone().insert_axis(Axis(1)));
                    let mut q_new: Array<f64, Ix1> = Array::zeros(self.neq);
                    for i in 0..self.neq {
                        q_new[i] = if left_elem_eigenvalues[i] >= 0.0 {alpha[[i, 0]]} else {beta[[i, 0]]};
                    }
                    q_new = left_elem_rmatrix.dot(&q_new);
                    let u_new = q_new[1] / q_new[0];
                    let v_new = q_new[2] / q_new[0];
                    let vn_new = u_new * nx + v_new * ny;
                    let h_new = (q_new[3] + p) / q_new[0];
                    let flux = [
                        q_new[0] * vn_new,
                        q_new[1] * vn_new + p * nx,
                        q_new[2] * vn_new + p * ny,
                        q_new[0] * h_new * vn_new
                    ];
                    for ivar in 0..self.neq {
                        residuals[[ielem, ivar, ibasis]] -= flux[ivar] * basis.phis_edge_gps[[ielem, igp, ibasis]] * 0.5 * mesh.boundary_edges[iedge].jacob_det;
                    }
                }
            }
        }
    }
    fn compute_boundary_edges_values(&self, edge: &BoundaryEdge, solutions: &Array<f64, Ix3>, basis: &DubinerBasis) -> Array<f64, Ix2> {
        let mut edges_values = Array::zeros((self.ngp, self.neq));
        let internal_element = edge.ielement;
        let ielem = edge.in_cell_index;
        let left_sol = solutions.slice(s![internal_element, .., ..]);
        for igp in 0..self.ngp {
            for ivar in 0..self.neq {
                for ibasis in 0..self.nbasis {
                    edges_values[[igp, ivar]] += left_sol[[ivar, ibasis]] * basis.phis_edge_gps[[ielem, igp, ibasis]];
                }
            }
        }
        edges_values
    }
}
