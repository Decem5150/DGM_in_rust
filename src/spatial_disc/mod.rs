pub mod flux;
pub mod quadrature;
pub mod gauss_point;
pub mod basis_function;
pub mod boundary;
pub use crate::mesh::{Element, Mesh};
pub use crate::solver::{SolverParameters, ConsVar, FlowParameters};
pub struct SpatialDisc<'a> {
    pub gauss_point: gauss_point::GaussPoints,
    pub inviscid_flux: Box<dyn flux::InvisFluxScheme<'a>>,
    // pub viscous_flux: Box<dyn flux::VisFluxScheme>,
    pub basis: basis_function::DubinerBasis,
    pub mesh: &'a Mesh<'a>,
    pub solver_param: &'a SolverParameters,
    pub flow_param: &'a FlowParameters,
}
impl<'a> SpatialDisc<'a> {
    pub fn integrate_over_cell(&self) {
        let nbasis = self.solver_param.number_of_basis_functions;
        let ngp = self.solver_param.number_of_cell_gp;
        for element in self.mesh.elements.iter() {
            let mut sol_gp = vec![ConsVar::default(); ngp];
            for (value_gp, phi_gp) in sol_gp.iter_mut().zip(self.basis.phi_cell_gp.iter()) {
                for (phi, sol) in phi_gp.iter().zip(element.solution.iter()) {
                    value_gp.density += sol.density * phi;
                    value_gp.x_momentum += sol.x_momentum * phi;
                    value_gp.y_momentum += sol.y_momentum * phi;
                    value_gp.energy += sol.energy * phi;
                }
            }
            for value_gp in sol_gp.iter() {
                let (f1, f2, f3, f4, g1, g2, g3, g4) = flux::flux(
                    value_gp.density, 
                    value_gp.x_momentum, 
                    value_gp.y_momentum, 
                    value_gp.energy, 
                    self.flow_param.hcr
                );
                
            }
        }
    }
    pub fn compute_fluxes(&self) {
        for edge in self.mesh.edges.iter() {
            let left_element = edge.elements[0];
            let right_element = edge.elements[1];
            let left_value = left_element.solution;
            let right_value = right_element.solution;
            self.inviscid_flux.compute(left_value, right_value, edge.invis_flux, edge.normal);
        }
    }
}