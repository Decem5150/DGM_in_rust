pub mod flux;
pub mod boundary_condition;
//pub mod limiter;
pub mod local_characteristics;
use ndarray::Array;
use ndarray::{Ix2, Ix3};
use ndarray::s;
use crate::basis_function::{DubinerBasis, GaussPoints};
use crate::mesh::Mesh;
use crate::solver::{FlowParameters, SolverParameters, MeshParameters};
pub enum InviscidFluxScheme {
    HLLC,
}
pub struct SpatialDisc<'a> {
    pub inviscid_flux_scheme: InviscidFluxScheme,
    // pub viscous_flux: Box<dyn flux::VisFluxScheme>,
    pub mesh: &'a Mesh,
    pub basis: &'a DubinerBasis,
    pub gauss_points: &'a GaussPoints,
    pub flow_param: &'a FlowParameters,
    pub mesh_param: &'a MeshParameters,
    pub solver_param: &'a SolverParameters,
}
impl<'a> SpatialDisc<'a> {
    pub fn compute_residuals(&self, residuals: &mut Array<f64, Ix3>, solutions: &Array<f64, Ix3>) {
        self.integrate_over_cell(residuals, solutions);
        self.integrate_over_edges(residuals, solutions);
        self.apply_bc(residuals, solutions);
        self.divide_residual_by_mass_mat_diag(residuals);
    }
    fn integrate_over_cell(&self, residuals: &mut Array<f64, Ix3>, solutions: &Array<f64, Ix3>) {
        let nbasis = self.solver_param.number_of_basis_functions;
        let ngp = self.solver_param.number_of_cell_gp;
        let neq = self.solver_param.number_of_equations;
        let nelem = self.mesh_param.number_of_elements;
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
                        residuals[[ielem, ivar, ibasis]] += (f[ivar] * self.mesh.elements[ielem].derivatives[[igp, ibasis]].get(&(1, 0)).unwrap()
                            + g[ivar] * self.mesh.elements[ielem].derivatives[[igp, ibasis]].get(&(0, 1)).unwrap())
                            * self.gauss_points.cell_weights[igp] * 0.5 * self.mesh.elements[ielem].jacob_det;
                    }
                }
            }
        }
    }
    fn integrate_over_edges(&self, residuals: &mut Array<f64, Ix3>, solutions: &Array<f64, Ix3>) {
        let nbasis = self.solver_param.number_of_basis_functions;
        let ngp = self.solver_param.number_of_edge_gp;
        let neq = self.solver_param.number_of_equations;
        for edge in self.mesh.edges.iter() {
            let ilelem = edge.in_cell_indices[0];
            let irelem = edge.in_cell_indices[1];
            let mut left_values_gps: Array<f64, Ix2> = Array::zeros((ngp, neq));
            let mut right_values_gps: Array<f64, Ix2> = Array::zeros((ngp, neq));
            let left_element = edge.ielements[0];
            let right_element = edge.ielements[1];
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
            for igp in 0..ngp {
                let num_flux = match self.inviscid_flux_scheme {
                    InviscidFluxScheme::HLLC => flux::hllc(
                        left_values_gps.slice(s![igp, ..]), 
                        right_values_gps.slice(s![igp, ..]), 
                        &edge.normal, 
                        &self.flow_param.hcr
                    ),
                };
                for ivar in 0..neq {
                    for ibasis in 0..nbasis {
                        residuals[[ilelem, ivar, ibasis]] -= num_flux[ivar] * self.basis.phis_edge_gps[[ilelem, igp, ibasis]] * 0.5 * edge.jacob_det;
                        residuals[[irelem, ivar, ibasis]] += num_flux[ivar] * self.basis.phis_edge_gps[[irelem, igp, ibasis]] * 0.5 * edge.jacob_det;
                    }
                }
            }
        }
    }
    fn divide_residual_by_mass_mat_diag(&self, residuals: &mut Array<f64, Ix3>) {
        let nbasis = self.solver_param.number_of_basis_functions;
        let neq = self.solver_param.number_of_equations;
        let nelem = self.mesh_param.number_of_elements;
        for ielem in 0..nelem {
            for ivar in 0..neq {
                for ibasis in 0..nbasis {
                    residuals[[ielem, ivar, ibasis]] /= self.mesh.elements[ielem].mass_mat_diag[ibasis];
                }
            }
        }
    }
}