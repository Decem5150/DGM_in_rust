pub mod flux;
pub mod quadrature;
pub mod gauss_point;
pub mod basis_function;
pub mod boundary;
pub use crate::mesh::{Element, Edge, Mesh};
pub use crate::solver::{SolverParameters, ConsVar, FlowParameters};
use itertools::izip;
pub struct SpatialDisc<'a> {
    pub gauss_point: gauss_point::GaussPoints,
    pub inviscid_flux: Box<dyn flux::InvisFluxScheme<'a>>,
    // pub viscous_flux: Box<dyn flux::VisFluxScheme>,
    pub basis: basis_function::DubinerBasis<'a>,
    pub mesh: &'a Mesh<'a>,
    pub solver_param: &'a SolverParameters,
    pub flow_param: &'a FlowParameters,
}
impl<'a> SpatialDisc<'_> {
    pub fn integrate_over_cell(&self) {
        let nbasis = self.solver_param.number_of_basis_functions;
        let ngp = self.solver_param.number_of_cell_gp;
        for element in self.mesh.elements.iter() {
            let mut sol_gp = vec![ConsVar::default(); ngp];
            for (value_gp, phi_gp) in sol_gp.iter_mut().zip(self.basis.phi_cell_gp.iter()) {
                for (&phi, &density, &x_momentum, &y_momentum, &energy) in izip!(
                    phi_gp, 
                    &element.solution.density, 
                    &element.solution.x_momentum, 
                    &element.solution.y_momentum, 
                    &element.solution.energy
                ) {
                    value_gp.density += density * phi;
                    value_gp.x_momentum += x_momentum * phi;
                    value_gp.y_momentum += y_momentum * phi;
                    value_gp.energy += energy * phi;
                }
            }
            for (value_gp, weight, dphis_dx_gps, dphis_dy_gps)  in izip!(&mut sol_gp, &self.gauss_point.cell_weights, &element.dphis_dx, &element.dphis_dy) {
                let (f1, f2, f3, f4, g1, g2, g3, g4) = flux::flux(
                    value_gp.density, 
                    value_gp.x_momentum, 
                    value_gp.y_momentum, 
                    value_gp.energy, 
                    self.flow_param.hcr
                );
                for (residual, &dphi_dx, &dphi_dy) in izip!(&mut element.residual, dphis_dx_gps, dphis_dy_gps) {
                    residual.density += (f1 * dphi_dx + g1 * dphi_dy) * weight * 0.5 * element.jacob_det;
                    residual.x_momentum += (f2 * dphi_dx + g2 * dphi_dy) * weight * 0.5 * element.jacob_det;
                    residual.y_momentum += (f3 * dphi_dx + g3 * dphi_dy) * weight * 0.5 * element.jacob_det;
                    residual.energy += (f4 * dphi_dx + g4 * dphi_dy) * weight * 0.5 * element.jacob_det;
                }
            }
        }
    }
    pub fn integrate_over_edges(&self) {
        for edge in self.mesh.edges.iter() {
            let left_element = edge.elements[0];
            let right_element = edge.elements[1];
            let indl = edge.ind_in_left_elem;
            let indr = edge.ind_in_right_elem;
            let mut left_values_gps = vec![ConsVar::default(); self.solver_param.number_of_edge_gp];
            let mut right_values_gps = vec![ConsVar::default(); self.solver_param.number_of_edge_gp];
            let left_phis_gps_iter = self.basis.phi_edge_gp[indl].iter();
            let right_phis_gps_iter = self.basis.phi_edge_gp[indr].iter();
            self.compute_edges_values(edge, &mut left_values_gps, &mut right_values_gps);
            for ((left_value_gp, right_value_gp), (left_phi_gp, right_phi_gp)) in 
            left_values_gps.iter().zip(right_values_gps.iter()).zip(left_phis_gps_iter.zip(right_phis_gps_iter.rev())) {
                let num_flux = self.inviscid_flux.compute(left_value_gp, right_value_gp, &edge.normal, &self.flow_param.hcr);
                for (&flux, left_res_vars, right_res_vars) in izip!(&num_flux, &mut left_element.residual, &mut right_element.residual) {
                    for (&mut left_res, &mut right_res, &left_phi, &right_phi) in izip!(left_res_vars.iter_mut(), right_res_vars.iter_mut(), left_phi_gp, right_phi_gp) {
                        left_res -= flux * left_phi * 0.5 * edge.jacob_det;
                        right_res += flux * right_phi * 0.5 * edge.jacob_det;
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
        for (value_gp, phi_gp) in left_value_gp.iter_mut().zip(self.basis.phi_edge_gp[indl].iter()) {
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
        for (value_gp, phi_gp) in right_value_gp.iter_mut().zip(self.basis.phi_edge_gp[indr].iter().rev()) {
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