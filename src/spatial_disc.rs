pub mod flux;
pub mod gauss_point;
pub mod basis_function;
pub mod boundary;
use ndarray::Array;
use ndarray::{ArrayView, ArrayViewMut};
use ndarray::{Ix2, Ix3};
use ndarray::{array, s};
pub use crate::mesh::{Element, Edge, Mesh, BoundaryEdge, Patch};
pub use crate::solver::{SolverParameters, FlowParameters};
pub struct SpatialDisc<'a> {
    pub gauss_point: gauss_point::GaussPoints,
    pub inviscid_flux: Box<dyn flux::InvisFluxScheme<'a>>,
    // pub viscous_flux: Box<dyn flux::VisFluxScheme>,
    // pub solutions: &'a Array<f64, Ix3>,
    pub basis: basis_function::DubinerBasis<'a>,
    pub mesh: &'a Mesh<'a>,
    pub solver_param: &'a SolverParameters,
    pub flow_param: &'a FlowParameters,
}
impl<'a> SpatialDisc<'_> {
    pub fn compute_residuals(&mut self, mut residuals: ArrayViewMut<f64, Ix3>, solutions: ArrayView<f64, Ix3>) {
        self.integrate_over_cell(residuals, solutions);
        self.integrate_over_edges(residuals, solutions);
        self.divide_residual_by_mass_mat_diag(residuals);
    }
    pub fn integrate_over_cell(&self, mut residuals: ArrayViewMut<f64, Ix3>, solutions: ArrayView<f64, Ix3>) {
        let nbasis = self.solver_param.number_of_basis_functions;
        let ngp = self.solver_param.number_of_cell_gp;
        let neq = self.solver_param.number_of_equations;
        let nelem = self.solver_param.number_of_elements;
        for ielem in 0..nelem {
            let mut sol_gp:Array<f64, Ix2> = Array::zeros((ngp, neq));
            for igp in 0..ngp {
                for ivar in 0..neq {
                    for ibasis in 0..nbasis {
                        sol_gp[[igp, ivar]] += solutions[[ielem, ivar, ibasis]] * self.basis.phis_cell_gps[[igp, ibasis]];
                    }
                }
            }
            for igp in 0..ngp {
                let q = sol_gp.slice(s![igp, ..]);
                let (f, g) = flux::flux(
                    q, 
                    self.flow_param.hcr
                );
                for ivar in 0..neq {
                    for ibasis in 0..nbasis {
                        residuals[[ielem, ivar, ibasis]] += (f[ivar] * self.basis.dphis_dxi[[igp, ibasis]] + g[ivar] * self.basis.dphis_deta[[igp, ibasis]]) * self.gauss_point.cell_weights[igp] * 0.5 * self.mesh.elements[ielem].jacob_det;
                    }
                }
            }
        }
    }
    pub fn integrate_over_edges(&self, mut residuals: ArrayViewMut<f64, Ix3>, solutions: ArrayView<f64, Ix3>) {
        let nbasis = self.solver_param.number_of_basis_functions;
        let ngp = self.solver_param.number_of_edge_gp;
        let neq = self.solver_param.number_of_equations;
        for edge in self.mesh.edges.iter() {
            let ilelem = edge.ind_in_left_elem;
            let irelem = edge.ind_in_right_elem;
            let mut left_values_gps = Array::zeros([ngp, neq]);
            let mut right_values_gps = Array::zeros([ngp, neq]);
            self.compute_edges_values(edge, solutions, left_values_gps.view_mut(), right_values_gps.view_mut());
            for igp in 0..ngp {
                let num_flux = self.inviscid_flux.compute(
                    left_values_gps.slice(s![igp, ..]), 
                    right_values_gps.slice(s![igp, ..]), 
                    &edge.normal, 
                    &self.flow_param.hcr
                );
                for ivar in 0..neq {
                    for ibasis in 0..nbasis {
                        residuals[[ilelem, ivar, ibasis]] -= num_flux[ivar] * self.basis.phis_edge_gps[[ilelem, igp, ibasis]] * 0.5 * edge.jacob_det;
                        residuals[[irelem, ivar, ibasis]] += num_flux[ivar] * self.basis.phis_edge_gps[[irelem, igp, ibasis]] * 0.5 * edge.jacob_det;
                    }
                }
            }
        }
    }
    pub fn integrate_over_boundary_edges(&self, mut residuals: ArrayViewMut<f64, Ix3>, solutions: ArrayView<f64, Ix3>) {
        let nbasis = self.solver_param.number_of_basis_functions;
        let ngp = self.solver_param.number_of_edge_gp;
        let neq = self.solver_param.number_of_equations;
        for patch in self.mesh.patches.iter() {
            for edge in patch.boundary_edges.iter() {
                let ielem = edge.ind_in_internal_elem;
                let mut left_values_gps = Array::zeros([ngp, neq]);
                let mut right_values_gps = Array::zeros([ngp, neq]);
                self.compute_boundary_edges_values(edge, solutions, left_values_gps.view_mut());
                // patch.
                for igp in 0..ngp {
                    let num_flux = self.inviscid_flux.compute(
                        left_values_gps.slice(s![igp, ..]), 
                        right_values_gps.slice(s![igp, ..]), 
                        &edge.normal, 
                        &self.flow_param.hcr
                    );
                    for ivar in 0..neq {
                        for ibasis in 0..nbasis {
                            residuals[[ielem, ivar, ibasis]] -= num_flux[ivar] * self.basis.phis_edge_gps[[ielem, igp, ibasis]] * 0.5 * edge.jacob_det;
                        }
                    }
                }
            }
        }
    }
    pub fn compute_edges_values(&self, edge: &Edge<'a>, solutions: ArrayView<f64, Ix3>, mut left_value_gp: ArrayViewMut<f64, Ix2>, mut right_value_gp: ArrayViewMut<f64, Ix2>) {
        let nbasis = self.solver_param.number_of_basis_functions;
        let ngp = self.solver_param.number_of_edge_gp;
        let neq = self.solver_param.number_of_equations;
        let left_element = edge.elements[0];
        let right_element = edge.elements[1];
        let ilelem = edge.ind_in_left_elem;
        let irelem = edge.ind_in_right_elem;
        let left_sol = solutions.slice(s![ilelem, .., ..]);
        let right_sol = solutions.slice(s![irelem, .., ..]);
        for igp in 0..ngp {
            for ivar in 0..neq {
                for ibasis in 0..nbasis {
                    left_value_gp[[igp, ivar]] += left_sol[[ivar, ibasis]] * self.basis.phis_edge_gps[[ilelem, igp, ibasis]];  
                }
            }
        }
        for igp in (0..ngp).rev() {
            for ivar in 0..neq {
                for ibasis in 0..nbasis {
                    right_value_gp[[igp, ivar]] += right_sol[[ivar, ibasis]] * self.basis.phis_edge_gps[[irelem, igp, ibasis]];
                }
            }
        }
    }
    pub fn compute_boundary_edges_values(&self, edge: &BoundaryEdge<'a>, solutions: ArrayView<f64, Ix3>, mut left_value_gp: ArrayViewMut<f64, Ix2>) {
        let nbasis = self.solver_param.number_of_basis_functions;
        let ngp = self.solver_param.number_of_edge_gp;
        let neq = self.solver_param.number_of_equations;
        let internal_element = edge.internal_element;
        let ielem = edge.ind_in_internal_elem;
        let left_sol = solutions.slice(s![ielem, .., ..]);
        for igp in 0..ngp {
            for ivar in 0..neq {
                for ibasis in 0..nbasis {
                    left_value_gp[[igp, ivar]] += left_sol[[ivar, ibasis]] * self.basis.phis_edge_gps[[ielem, igp, ibasis]];  
                }
            }
        }
    }
    pub fn divide_residual_by_mass_mat_diag(&mut self, mut residuals: ArrayViewMut<f64, Ix3>) {
        let nbasis = self.solver_param.number_of_basis_functions;
        let neq = self.solver_param.number_of_equations;
        let nelem = self.solver_param.number_of_elements;
        for ielem in 0..nelem {
            for ivar in 0..neq {
                for ibasis in 0..nbasis {
                    residuals[[ielem, ivar, ibasis]] /= self.mesh.elements[ielem].mass_mat_diag[ibasis];
                }
            }
        }
    }
}