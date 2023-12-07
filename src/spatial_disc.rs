pub mod flux;
pub mod boundary;
pub mod limiter;
pub mod local_characteristics;
use std::cell::RefCell;
use std::rc::{Weak, Rc};
use ndarray::Array;
use ndarray::{ArrayView, ArrayViewMut};
use ndarray::{Ix2, Ix3};
use ndarray::s;
use super::mesh::{Edge, Mesh};
use super::solver::{SolverParameters, FlowParameters};
use crate::basis_function::DubinerBasis;
use crate::gauss_point::GaussPoints;
pub struct SpatialDisc<'a> {
    pub inviscid_flux: Box<dyn flux::InvisFluxScheme<'a>>,
    // pub viscous_flux: Box<dyn flux::VisFluxScheme>,
    pub boundary_condition: boundary::BoundaryCondition<'a>,
    pub basis: Rc<DubinerBasis>,
    pub gauss_point: Rc<GaussPoints>,
    pub mesh: Rc<Mesh>,
    pub solver_param: Rc<SolverParameters>,
    pub flow_param: Rc<FlowParameters>,
}
impl<'a> SpatialDisc<'_> {
    pub fn compute_residuals(&mut self, mut residuals: ArrayViewMut<f64, Ix3>, solutions: ArrayView<f64, Ix3>) {
        self.integrate_over_cell(residuals, solutions);
        self.integrate_over_edges(residuals, solutions);
        self.boundary_condition.apply_bc(residuals, solutions);
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
                        residuals[[ielem, ivar, ibasis]] += (f[ivar] * self.basis.dphis_dxi[[igp, ibasis]] 
                            + g[ivar] * self.basis.dphis_deta[[igp, ibasis]]) 
                            * self.gauss_point.cell_weights[igp] * 0.5 * self.mesh.elements[ielem].jacob_det;
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
            let ilelem = edge.in_cell_index[0];
            let irelem = edge.in_cell_index[1];
            let mut left_values_gps: Array<f64, Ix2> = Array::zeros((ngp, neq));
            let mut right_values_gps: Array<f64, Ix2> = Array::zeros((ngp, neq));
            self.compute_edge_values(edge, solutions, left_values_gps.view_mut(), right_values_gps.view_mut());
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
    pub fn compute_edge_values(&self, edge: &Edge, solutions: ArrayView<f64, Ix3>, mut left_values_gps: ArrayViewMut<f64, Ix2>, mut right_values_gps: ArrayViewMut<f64, Ix2>) {
        let nbasis = self.solver_param.number_of_basis_functions;
        let ngp = self.solver_param.number_of_edge_gp;
        let neq = self.solver_param.number_of_equations;
        let left_element = edge.ielements[0];
        let right_element = edge.ielements[1];
        let ilelem = edge.in_cell_index[0];
        let irelem = edge.in_cell_index[1];
        let left_sol = solutions.slice(s![left_element, .., ..]);
        let right_sol = solutions.slice(s![right_element, .., ..]);
        for igp in 0..ngp {
            for ivar in 0..neq {
                for ibasis in 0..nbasis {
                    left_values_gps[[igp, ivar]] += left_sol[[ivar, ibasis]] * self.basis.phis_edge_gps[[ilelem, igp, ibasis]];  
                }
            }
        }
        for igp in 0..ngp {
            for ivar in 0..neq {
                for ibasis in 0..nbasis {
                    right_values_gps[[igp, ivar]] += right_sol[[ivar, ibasis]] * self.basis.phis_edge_gps[[irelem, ngp - 1 - igp, ibasis]];
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