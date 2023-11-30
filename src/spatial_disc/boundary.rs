use ndarray::Array;
use ndarray::{Ix1, Ix2, Ix3};
use ndarray::{ArrayView, ArrayViewMut, Axis};
use ndarray::s;
use crate::mesh::{Element, Edge, BoundaryEdge, BoundaryType, Patch};
use crate::solver::{SolverParameters, FlowParameters};
use crate::spatial_disc::basis_function;
use crate::spatial_disc::local_characteristics;
pub struct BoundaryCondition<'a> {
    pub basis: basis_function::DubinerBasis<'a>,
    pub solver_param: &'a SolverParameters,
    pub flow_param: &'a FlowParameters,
    pub patches: &'a Array<Patch<'a>, Ix1>,
}
impl<'a> BoundaryCondition<'a> {
    pub fn apply_bc(&self, solutions: ArrayView<f64, Ix3>, residuals: ArrayViewMut<f64, Ix3>) {
        for patch in self.patches.iter() {
            match patch.boundary_type {
                BoundaryType::NoSlipWall => self.no_slip_wall(solutions, residuals, patch),
                BoundaryType::FreeSlipWall => self.free_slip_wall(solutions, residuals, patch),
                BoundaryType::FarField => self.far_field(solutions, residuals, patch),
                _ => panic!("Boundary type not implemented"),
            }
        }
    }
    pub fn no_slip_wall(&self, solutions: ArrayView<f64, Ix3>, mut residuals: ArrayViewMut<f64, Ix3>, patch: &Patch) {
        let nbasis = self.solver_param.number_of_basis_functions;
        let ngp = self.solver_param.number_of_edge_gp;
        let neq = self.solver_param.number_of_equations;
        let hcr = self.flow_param.hcr;
        for edge in patch.edges.iter() {
            let ielem = edge.ind_in_internal_elem;
            let nx = edge.normal[0];
            let ny = edge.normal[1];
            let mut left_values_gps: Array<f64, Ix2> = Array::zeros((ngp, neq));
            self.compute_boundary_edges_values(edge, solutions, left_values_gps.view_mut()); 
            for igp in 0..ngp {
                for ibasis in 0..nbasis {
                    let pressure = (hcr - 1.0) * (left_values_gps[[igp, 3]] - 0.5 * (left_values_gps[[igp, 1]] * left_values_gps[[igp, 1]] + left_values_gps[[igp, 2]] * left_values_gps[[igp, 2]]) / left_values_gps[[igp, 0]]);
                    residuals[[ielem, 1, ibasis]] -= pressure * self.basis.phis_edge_gps[[ielem, igp, ibasis]] * nx * 0.5 * edge.jacob_det;
                    residuals[[ielem, 2, ibasis]] -= pressure * self.basis.phis_edge_gps[[ielem, igp, ibasis]] * ny * 0.5 * edge.jacob_det;
                }
            }
        }
    }
    pub fn free_slip_wall(&self, solutions: ArrayView<f64, Ix3>, mut residuals: ArrayViewMut<f64, Ix3>, patch: &Patch) {
        let nbasis = self.solver_param.number_of_basis_functions;
        let ngp = self.solver_param.number_of_edge_gp;
        let neq = self.solver_param.number_of_equations;
        let hcr = self.flow_param.hcr;
        for edge in patch.edges.iter() {
            let ielem = edge.ind_in_internal_elem;
            let nx = edge.normal[0];
            let ny = edge.normal[1];
            let mut left_values_gps: Array<f64, Ix2> = Array::zeros((ngp, neq));
            self.compute_boundary_edges_values(edge, solutions, left_values_gps.view_mut()); 
            for igp in 0..ngp {
                for ibasis in 0..nbasis {
                    let pressure = (hcr - 1.0) * (left_values_gps[[igp, 3]] - 0.5 * (left_values_gps[[igp, 1]] * left_values_gps[[igp, 1]] + left_values_gps[[igp, 2]] * left_values_gps[[igp, 2]]) / left_values_gps[[igp, 0]]);
                    residuals[[ielem, 1, ibasis]] -= pressure * self.basis.phis_edge_gps[[ielem, igp, ibasis]] * nx * 0.5 * edge.jacob_det;
                    residuals[[ielem, 2, ibasis]] -= pressure * self.basis.phis_edge_gps[[ielem, igp, ibasis]] * ny * 0.5 * edge.jacob_det;
                }
            }
        }
    }
    pub fn far_field(&self, solutions: ArrayView<f64, Ix3>, residuals: ArrayViewMut<f64, Ix3>, patch: &Patch) {
        let nbasis = self.solver_param.number_of_basis_functions;
        let ngp = self.solver_param.number_of_edge_gp;
        let neq = self.solver_param.number_of_equations;
        let hcr = self.flow_param.hcr;
        let boundary_quantity = patch.boundary_quantity.unwrap();
        let boundary_consvar = Array::from_vec(vec![boundary_quantity.rho, boundary_quantity.rho * boundary_quantity.u, boundary_quantity.rho * boundary_quantity.v, boundary_quantity.p / (hcr - 1.0) + 0.5 * boundary_quantity.rho * (boundary_quantity.u.powi(2) + boundary_quantity.v.powi(2))]);
        for edge in patch.edges.iter() {
            let ielem = edge.ind_in_internal_elem;
            let nx = edge.normal[0];
            let ny = edge.normal[1];
            let mut left_values_gps: Array<f64, Ix2> = Array::zeros((ngp, neq));
            self.compute_boundary_edges_values(edge, solutions, left_values_gps.view_mut()); 
            for igp in 0..ngp {
                for ibasis in 0..nbasis {
                    let u = left_values_gps[[igp, 1]] / left_values_gps[[igp, 0]];
                    let v = left_values_gps[[igp, 2]] / left_values_gps[[igp, 0]];
                    let p = (hcr - 1.0) * (left_values_gps[[igp, 3]] - 0.5 * (left_values_gps[[igp, 1]] * left_values_gps[[igp, 1]] + left_values_gps[[igp, 2]] * left_values_gps[[igp, 2]]) / left_values_gps[[igp, 0]]);
                    let c = (hcr * p / left_values_gps[[igp, 0]]).sqrt();
                    let vn = u * nx + v * ny;
                    let (left_elem_lmatrix, left_elem_rmatrix) = local_characteristics::compute_eigenmatrix(left_values_gps.slice(s![igp, ..]), nx, ny, hcr);
                    let left_elem_eigenvalues = [vn - c, vn, vn, vn + c];
                    let alpha = left_elem_lmatrix.dot(&left_values_gps.slice(s![igp, ..]).insert_axis(Axis(1)));
                    let beta = left_elem_lmatrix.dot(&boundary_consvar.insert_axis(Axis(1)));
                    let mut q_new: Array<f64, Ix1> = Array::zeros(neq);
                    for i in 0..neq {
                        q_new[i] = if left_elem_eigenvalues[i] >= 0.0 {
                                alpha[[i, 0]]
                            } else {
                                beta[[i, 0]]
                            };
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
                    for ivar in 0..neq {
                        residuals[[ielem, ivar, ibasis]] -= flux[ivar] * self.basis.phis_edge_gps[[ielem, igp, ibasis]] * 0.5 * edge.jacob_det;
                    }
                }
            }
        }
    }
    pub fn compute_boundary_edges_values(&self, edge: &BoundaryEdge<'a>, solutions: ArrayView<f64, Ix3>, mut left_values_gps: ArrayViewMut<f64, Ix2>) {
        let nbasis = self.solver_param.number_of_basis_functions;
        let ngp = self.solver_param.number_of_edge_gp;
        let neq = self.solver_param.number_of_equations;
        let internal_element = edge.internal_element;
        let ielem = edge.ind_in_internal_elem;
        let left_sol = solutions.slice(s![ielem, .., ..]);
        for igp in 0..ngp {
            for ivar in 0..neq {
                for ibasis in 0..nbasis {
                    left_values_gps[[igp, ivar]] += left_sol[[ivar, ibasis]] * self.basis.phis_edge_gps[[ielem, igp, ibasis]];  
                }
            }
        }
    }
}
