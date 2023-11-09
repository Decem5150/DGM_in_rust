use ndarray::Array;
use ndarray::{Ix2, Ix3};
use ndarray::{array, s};
pub mod flux;
pub mod quadrature;
pub mod gauss_point;
pub mod basis_function;
pub mod boundary;
pub use crate::mesh::{Element, Edge, Mesh};
pub use crate::solver::{SolverParameters, FlowParameters, };
use itertools::izip;
use std::iter::zip;
pub struct SpatialDisc<'a> {
    pub gauss_point: gauss_point::GaussPoints,
    pub inviscid_flux: Box<dyn flux::InvisFluxScheme<'a>>,
    // pub viscous_flux: Box<dyn flux::VisFluxScheme>,
    pub solutions: &'a Array<f64, Ix3>,
    pub basis: basis_function::DubinerBasis<'a>,
    pub mesh: &'a Mesh<'a>,
    pub solver_param: &'a SolverParameters,
    pub flow_param: &'a FlowParameters,
}
impl<'a> SpatialDisc<'_> {
    pub fn compute_residuals(&mut self, residuals: &mut Array<f64, Ix3>) {
        self.integrate_over_cell(residuals);
        self.integrate_over_edges(residuals);
        for element in self.mesh.elements.iter_mut() {
            element.divide_residual_by_mass_mat_diag(residuals);
        }
    }
    pub fn integrate_over_cell(&self, residuals: &mut Array<f64, Ix3>) {
        let nbasis = self.solver_param.number_of_basis_functions;
        let ngp = self.solver_param.number_of_cell_gp;
        let neq = self.solver_param.number_of_equations;
        let nelem = self.solver_param.number_of_elements;
        for ielem in 0..nelem {
            let mut sol_gp:Array<f64, Ix2> = Array::zeros((ngp, neq));
            for igp in 0..ngp {
                for ivar in 0..neq {
                    for ibasis in 0..nbasis {
                        sol_gp[[igp, ivar]] += self.solutions[[ielem, ivar, ibasis]] * self.basis.phis_cell_gps[[igp, ibasis]];
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
    pub fn integrate_over_edges(&self, residuals: &mut Array<f64, Ix3>) {
        for edge in self.mesh.edges.iter() {
            let left_element = edge.elements[0];
            let right_element = edge.elements[1];
            let indl = edge.ind_in_left_elem;
            let indr = edge.ind_in_right_elem;
            let mut left_values_gps = vec![ConsVar::default(); self.solver_param.number_of_edge_gp];
            let mut right_values_gps = vec![ConsVar::default(); self.solver_param.number_of_edge_gp];
            let left_phis_gps_iter = self.basis.phis_edge_gps[indl].iter();
            let right_phis_gps_iter = self.basis.phis_edge_gps[indr].iter();
            self.compute_edges_values(edge, &mut left_values_gps, &mut right_values_gps);
            for ((left_value_gp, right_value_gp), (left_phi_gp, right_phi_gp)) in 
            left_values_gps.iter().zip(right_values_gps.iter()).zip(left_phis_gps_iter.zip(right_phis_gps_iter.rev())) {
                let num_flux = self.inviscid_flux.compute(left_value_gp, right_value_gp, &edge.normal, &self.flow_param.hcr);
                for (&flux, left_res_vars, right_res_vars) in izip!(&num_flux, &mut left_element.residual, &mut right_element.residual) {
                    for (left_res, right_res, &left_phi, &right_phi) in izip!(left_res_vars.iter_mut(), right_res_vars.iter_mut(), left_phi_gp, right_phi_gp) {
                        *left_res -= flux * left_phi * 0.5 * edge.jacob_det;
                        *right_res += flux * right_phi * 0.5 * edge.jacob_det;
                    }
                }
            }
        }
    }
    pub fn compute_edges_values(&self, edge: &Edge<'a>, left_value_gp: &mut Vec<ConsVar>, right_value_gp: &mut Vec<ConsVar>) {
        let left_element = edge.elements[0];
        let right_element = edge.elements[1];
        let indl = edge.ind_in_left_elem;
        let indr = edge.ind_in_right_elem;
        let left_sol = left_element.solution;
        let right_sol = right_element.solution;
        for (value_gp, phi_gp) in left_value_gp.iter_mut().zip(self.basis.phis_edge_gps[indl].iter()) {
            for (&phi, &density, &x_momentum, &y_momentum, &energy) in izip!(
                phi_gp, 
                &left_sol.density, 
                &left_sol.x_momentum, 
                &left_sol.y_momentum, 
                &left_sol.energy
            ) {
                value_gp.density += density * phi;
                value_gp.x_momentum += x_momentum * phi;
                value_gp.y_momentum += y_momentum * phi;
                value_gp.energy += energy * phi;
            }
        }
        for (value_gp, phi_gp) in right_value_gp.iter_mut().zip(self.basis.phis_edge_gps[indr].iter().rev()) {
            for (&phi, &density, &x_momentum, &y_momentum, &energy) in izip!(
                phi_gp, 
                &right_sol.density, 
                &right_sol.x_momentum, 
                &right_sol.y_momentum, 
                &right_sol.energy
            ) {
                value_gp.density += density * phi;
                value_gp.x_momentum += x_momentum * phi;
                value_gp.y_momentum += y_momentum * phi;
                value_gp.energy += energy * phi;
            }
        }
    }
}