pub mod flux;
pub mod boundary_condition;
pub mod limiter;
pub mod local_characteristics;
pub mod viscous;
use ndarray::{Array, Ix1, Ix4};
use ndarray::{Ix2, Ix3};
use ndarray::s;
use crate::basis_function::{DubinerBasis, GaussPoints};
use crate::debug_utils::check_for_nan;
use crate::mesh::{BoundaryType, EdgeTypeAndIndex, Mesh};
use crate::solver::{FlowParameters, SolverParameters, MeshParameters};
pub enum InviscidFluxScheme {
    HLLC,
}
pub struct SpatialDisc<'a> {
    // pub viscous_flux: Box<dyn flux::VisFluxScheme>,
    pub aux_vars: Array<f64, Ix4>,
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
        let ngp = self.solver_param.number_of_cell_gp;
        let neq = self.solver_param.number_of_equations;
        let nbasis = self.solver_param.number_of_basis_functions;
        let weights = &self.gauss_points.cell_weights;
        for &ielem in self.mesh.internal_element_indices.iter() {
            let mut global_lift_x: Array<f64, Ix2> = Array::zeros((neq, nbasis));
            let mut global_lift_y: Array<f64, Ix2> = Array::zeros((neq, nbasis));
            // solve for global lift
            let edge_weights = &self.gauss_points.edge_weights;
            for ivar in 0..neq {
                for ibasis in 0..nbasis {
                    // rhs
                    for (ilic, iedge) in self.mesh.elements[ielem].iedges.indexed_iter() {
                        let edge = &self.mesh.edges[*iedge];
                        let ilelem = ielem as usize;
                        let (iric, irelem, nx, ny) = {
                            if edge.in_cell_indices[0] as usize == ilic {
                                (edge.in_cell_indices[1] as usize, edge.ielements[1] as usize, edge.normal[0], edge.normal[1])
                            } else {
                                (edge.in_cell_indices[0] as usize, edge.ielements[0] as usize, -edge.normal[0], -edge.normal[1])
                            }
                        };
                        for igp in 0..self.solver_param.number_of_edge_gp {
                            let mut left_value = 0.0;
                            let mut right_value = 0.0;
                            for ibasis in 0..nbasis {
                                left_value += solutions[[ilelem, ivar, ibasis]] * self.basis.phis_edge_gps[[ilic, igp, ibasis]];
                                right_value += solutions[[irelem, ivar, ibasis]] * self.basis.phis_edge_gps[[iric, ngp - 1 - igp, ibasis]];
                            }
                            let jump_x = (left_value - right_value) * nx;
                            let jump_y = (left_value - right_value) * ny;
                            global_lift_x[[ivar, ibasis]] -= 0.5 * jump_x * self.basis.phis_edge_gps[[ilic, igp, ibasis]] * edge_weights[igp] * 0.5 * edge.jacob_det;
                            global_lift_y[[ivar, ibasis]] -= 0.5 * jump_y * self.basis.phis_edge_gps[[ilic, igp, ibasis]] * edge_weights[igp] * 0.5 * edge.jacob_det;
                        }
                            /*
                            EdgeTypeAndIndex::Boundary(iedge) => {
                                let edge = &self.mesh.boundary_edges[*iedge];
                                let ilelem = ielem;
                                let (nx, ny) = {
                                    if self.mesh.elements[ielem].ivertices[ilic] == edge.ivertices[0] {
                                        (edge.normal[0], edge.normal[1])
                                    } else {
                                        (-edge.normal[0], -edge.normal[1])
                                    }
                                };
                                match self.mesh.patches[edge.ipatch].boundary_type {
                                    BoundaryType::Wall => {
                                        for igp in 0..self.solver_param.number_of_edge_gp {
                                            let mut left_value = 0.0;
                                            for ibasis in 0..nbasis {
                                                if ivar == 1||ivar == 2 {
                                                    left_value += solutions[[ilelem, ivar, ibasis]] * self.basis.phis_edge_gps[[ilic, igp, ibasis]];
                                                }
                                            }
                                            let jump_x = left_value * nx;
                                            let jump_y = left_value * ny;
                                            global_lift_x[[ivar, ibasis]] -= 0.5 * jump_x * self.basis.phis_edge_gps[[ilic, igp, ibasis]] * edge_weights[igp] * 0.5 * edge.jacob_det;
                                            global_lift_y[[ivar, ibasis]] -= 0.5 * jump_y * self.basis.phis_edge_gps[[ilic, igp, ibasis]] * edge_weights[igp] * 0.5 * edge.jacob_det;
                                        }
                                    }
                                    BoundaryType::FarField => {

                                    }
                                }
                            }
                            */
                    }
                    global_lift_x[[ivar, ibasis]] /= self.basis.mass_mat_diag[ibasis] * 0.5 * self.mesh.elements[ielem].jacob_det;
                    global_lift_y[[ivar, ibasis]] /= self.basis.mass_mat_diag[ibasis] * 0.5 * self.mesh.elements[ielem].jacob_det;
                }
            }
            for igp in 0..ngp {
                let mut sol: Array<f64, Ix1> = Array::zeros(neq);
                for ivar in 0..solutions.shape()[1] {
                    for ibasis in 0..solutions.shape()[2] {
                        sol[ivar] += solutions[[ielem, ivar, ibasis]] * self.basis.phis_cell_gps[[igp, ibasis]];
                    }
                }
                let (f, g) = flux::flux(
                    &sol, 
                    self.flow_param.hcr
                );
                let mut fv: Array<f64, Ix1> = Array::zeros(neq);
                let mut gv: Array<f64, Ix1> = Array::zeros(neq);
                let viscosity = {
                    let pressure = (self.flow_param.hcr - 1.0) * (sol[3] - 0.5 * (sol[1] * sol[1] + sol[2] * sol[2]) / sol[0]);
                    let temperature = pressure / (self.flow_param.gas_constant * sol[0]);
                    Self::compute_dynamic_viscosity(temperature)
                };
                let (g11, g12, g21, g22) = self.viscous_jacobian(&sol, viscosity);
                for iflux in 0..neq {
                    for ivar in 0..neq {
                        let mut dsol_dx = 0.0;
                        let mut dsol_dy = 0.0;
                        let mut lift_x = 0.0;
                        let mut lift_y = 0.0;
                        for ibasis in 0..solutions.shape()[2] {
                            let dphi_dx = *self.mesh.elements[ielem].dphis_cell_gps[[igp, ibasis]].get(&(1, 0)).unwrap();
                            let dphi_dy = *self.mesh.elements[ielem].dphis_cell_gps[[igp, ibasis]].get(&(0, 1)).unwrap();
                            dsol_dx += solutions[[ielem, ivar, ibasis]] * dphi_dx;
                            dsol_dy += solutions[[ielem, ivar, ibasis]] * dphi_dy;
                            lift_x += global_lift_x[[ivar, ibasis]] * self.basis.phis_cell_gps[[igp, ibasis]];
                            lift_y += global_lift_y[[ivar, ibasis]] * self.basis.phis_cell_gps[[igp, ibasis]];
                        }
                        fv[iflux] += g11[[iflux, ivar]] * (dsol_dx + lift_x) + g12[[iflux, ivar]] * (dsol_dy + lift_y);
                        gv[iflux] += g21[[iflux, ivar]] * (dsol_dx + lift_x) + g22[[iflux, ivar]] * (dsol_dy + lift_y);
                    }
                }
                for ivar in 0..solutions.shape()[1] {
                    for ibasis in 0..solutions.shape()[2] {
                        let dphi_dx = *self.mesh.elements[ielem].dphis_cell_gps[[igp, ibasis]].get(&(1, 0)).unwrap();
                        let dphi_dy = *self.mesh.elements[ielem].dphis_cell_gps[[igp, ibasis]].get(&(0, 1)).unwrap();
                        residuals[[ielem, ivar, ibasis]] += ((f[ivar] - fv[ivar]) * dphi_dx + (g[ivar] - gv[ivar]) * dphi_dy)
                            * weights[igp] * 0.5 * self.mesh.elements[ielem].jacob_det;
                    }
                }
            }
        }
    }
    fn integrate_over_edges(&self, residuals: &mut Array<f64, Ix3>, solutions: &Array<f64, Ix3>) {
        let ngp = self.solver_param.number_of_edge_gp;
        let neq = self.solver_param.number_of_equations;
        let nbasis = self.solver_param.number_of_basis_functions;
        let weights = &self.gauss_points.edge_weights;
        for iedge in self.mesh.internal_edge_indices.iter() {
            let edge = &self.mesh.edges[*iedge];
            let mut left_local_lift_x: Array<f64, Ix2> = Array::zeros((neq, nbasis));
            let mut left_local_lift_y: Array<f64, Ix2> = Array::zeros((neq, nbasis));
            let mut right_local_lift_x: Array<f64, Ix2> = Array::zeros((neq, nbasis));
            let mut right_local_lift_y: Array<f64, Ix2> = Array::zeros((neq, nbasis));
            let ilic = edge.in_cell_indices[0] as usize;
            let iric = edge.in_cell_indices[1] as usize;
            let ilelem = edge.ielements[0] as usize;
            let irelem = edge.ielements[1] as usize;
            let left_dofs = solutions.slice(s![ilelem, .., ..]);
            let right_dofs = solutions.slice(s![irelem, .., ..]);
            let nx = edge.normal[0];
            let ny = edge.normal[1];
            // solve for local lift
            for ivar in 0..neq {
                for ibasis in 0..nbasis {
                    for igp in 0..ngp {
                        let mut left_value = 0.0;
                        let mut right_value = 0.0;
                        for ibasis in 0..nbasis {
                            left_value += left_dofs[[ivar, ibasis]] * self.basis.phis_edge_gps[[ilic, igp, ibasis]];
                            right_value += right_dofs[[ivar, ibasis]] * self.basis.phis_edge_gps[[iric, ngp - 1 - igp, ibasis]];
                        }
                        let jump_x = (left_value - right_value) * nx;
                        let jump_y = (left_value - right_value) * ny;
                        left_local_lift_x[[ivar, ibasis]] -= 0.5 * jump_x * self.basis.phis_edge_gps[[ilic, igp, ibasis]] * weights[igp] * 0.5 * edge.jacob_det;
                        left_local_lift_y[[ivar, ibasis]] -= 0.5 * jump_y * self.basis.phis_edge_gps[[ilic, igp, ibasis]] * weights[igp] * 0.5 * edge.jacob_det;
                        right_local_lift_x[[ivar, ibasis]] += 0.5 * jump_x * self.basis.phis_edge_gps[[iric, ngp - 1 - igp, ibasis]] * weights[ngp - 1 - igp] * 0.5 * edge.jacob_det;
                        right_local_lift_y[[ivar, ibasis]] += 0.5 * jump_y * self.basis.phis_edge_gps[[iric, ngp - 1 - igp, ibasis]] * weights[ngp - 1 - igp] * 0.5 * edge.jacob_det;
                    }
                    left_local_lift_x[[ivar, ibasis]] /= self.basis.mass_mat_diag[ibasis] * 0.5 * self.mesh.elements[ilelem].jacob_det;
                    left_local_lift_y[[ivar, ibasis]] /= self.basis.mass_mat_diag[ibasis] * 0.5 * self.mesh.elements[ilelem].jacob_det;
                    right_local_lift_x[[ivar, ibasis]] /= self.basis.mass_mat_diag[ibasis] * 0.5 * self.mesh.elements[irelem].jacob_det;
                    right_local_lift_y[[ivar, ibasis]] /= self.basis.mass_mat_diag[ibasis] * 0.5 * self.mesh.elements[irelem].jacob_det;
                }
            }
            // compute numerical flux
            for igp in 0..ngp {
                // compute inviscid flux
                let mut left_values: Array<f64, Ix1> = Array::zeros(neq);
                let mut right_values: Array<f64, Ix1> = Array::zeros(neq);
                for ivar in 0..neq {
                    for ibasis in 0..nbasis {
                        left_values[ivar] += left_dofs[[ivar, ibasis]] * self.basis.phis_edge_gps[[ilic, igp, ibasis]];  
                        right_values[ivar] += right_dofs[[ivar, ibasis]] * self.basis.phis_edge_gps[[iric, ngp - 1 - igp, ibasis]];
                    }
                }
                let num_flux = match self.solver_param.inviscid_flux_scheme {
                    InviscidFluxScheme::HLLC => flux::hllc(
                        &left_values, 
                        &right_values, 
                        &edge.normal, 
                        &self.flow_param.hcr
                        ),
                };
                // compute viscous flux
                let mut num_viscous_flux: Array<f64, Ix1> = Array::zeros(neq);
                // left contribution to viscous flux
                let mut fv: Array<f64, Ix1> = Array::zeros(neq);
                let mut gv: Array<f64, Ix1> = Array::zeros(neq);
                let viscosity = {
                    let pressure = (self.flow_param.hcr - 1.0) * (left_values[3] - 0.5 * (left_values[1] * left_values[1] + left_values[2] * left_values[2]) / left_values[0]);
                    let temperature = pressure / (self.flow_param.gas_constant * left_values[0]);
                    Self::compute_dynamic_viscosity(temperature)
                };
                let (g11, g12, g21, g22) = self.viscous_jacobian(&left_values, viscosity);
                for iflux in 0..neq {
                    for ivar in 0..neq {
                        let mut dsol_dx = 0.0;
                        let mut dsol_dy = 0.0;
                        let mut lift_x = 0.0;
                        let mut lift_y = 0.0;
                        for ibasis in 0..solutions.shape()[2] {
                            let dphi_dx = *self.mesh.elements[ilelem].dphis_edge_gps[[ilic, igp, ibasis]].get(&(1, 0)).unwrap();
                            let dphi_dy = *self.mesh.elements[ilelem].dphis_edge_gps[[ilic, igp, ibasis]].get(&(0, 1)).unwrap();
                            dsol_dx += left_dofs[[ivar, ibasis]] * dphi_dx;
                            dsol_dy += left_dofs[[ivar, ibasis]] * dphi_dy;
                            lift_x += left_local_lift_x[[ivar, ibasis]] * self.basis.phis_edge_gps[[ilic, igp, ibasis]];
                            lift_y += left_local_lift_y[[ivar, ibasis]] * self.basis.phis_edge_gps[[ilic, igp, ibasis]];
                        }
                        fv += g11[[iflux, ivar]] * (dsol_dx + lift_x) + g12[[iflux, ivar]] * (dsol_dy + lift_y);
                        gv += g21[[iflux, ivar]] * (dsol_dx + lift_x) + g22[[iflux, ivar]] * (dsol_dy + lift_y);
                    }
                }
                num_viscous_flux = num_viscous_flux + 0.5 * (nx * &fv + ny * &gv);
                // right contribution to viscous flux
                fv = Array::zeros(neq);
                gv = Array::zeros(neq);
                let viscosity = {
                    let pressure = (self.flow_param.hcr - 1.0) * (right_values[3] - 0.5 * (right_values[1] * right_values[1] + right_values[2] * right_values[2]) / right_values[0]);
                    let temperature = pressure / (self.flow_param.gas_constant * right_values[0]);
                    Self::compute_dynamic_viscosity(temperature)
                };
                let (g11, g12, g21, g22) = self.viscous_jacobian(&right_values, viscosity);
                for iflux in 0..neq {
                    for ivar in 0..neq {
                        let mut dsol_dx = 0.0;
                        let mut dsol_dy = 0.0;
                        let mut lift_x = 0.0;
                        let mut lift_y = 0.0;
                        for ibasis in 0..solutions.shape()[2] {
                            let dphi_dx = *self.mesh.elements[irelem].dphis_edge_gps[[iric, ngp - 1 - igp, ibasis]].get(&(1, 0)).unwrap();
                            let dphi_dy = *self.mesh.elements[irelem].dphis_edge_gps[[iric, ngp - 1 - igp, ibasis]].get(&(0, 1)).unwrap();
                            dsol_dx += right_dofs[[ivar, ibasis]] * dphi_dx;
                            dsol_dy += right_dofs[[ivar, ibasis]] * dphi_dy;
                            lift_x += right_local_lift_x[[ivar, ibasis]] * self.basis.phis_edge_gps[[iric, ngp - 1 - igp, ibasis]];
                            lift_y += right_local_lift_y[[ivar, ibasis]] * self.basis.phis_edge_gps[[iric, ngp - 1 - igp, ibasis]];
                        }
                        fv += g11[[iflux, ivar]] * (dsol_dx + lift_x) + g12[[iflux, ivar]] * (dsol_dy + lift_y);
                        gv += g21[[iflux, ivar]] * (dsol_dx + lift_x) + g22[[iflux, ivar]] * (dsol_dy + lift_y);
                    }
                }
                num_viscous_flux = num_viscous_flux + 0.5 * (nx * fv + ny * gv);

                for ivar in 0..residuals.shape()[1] {
                    for ibasis in 0..residuals.shape()[2] {
                        residuals[[ilelem, ivar, ibasis]] -= weights[igp] * (num_flux[ivar] + num_viscous_flux[ivar]) * self.basis.phis_edge_gps[[ilic, igp, ibasis]] * 0.5 * edge.jacob_det;
                        residuals[[irelem, ivar, ibasis]] += weights[ngp - 1 - igp] * (num_flux[ivar] + num_viscous_flux[ivar]) * self.basis.phis_edge_gps[[iric, ngp - 1 - igp, ibasis]] * 0.5 * edge.jacob_det;
                    }
                }
            }
        }
    }
    fn divide_residual_by_mass_mat_diag(&self, residuals: &mut Array<f64, Ix3>) {
        for ielem in 0..residuals.shape()[0] {
            for ivar in 0..residuals.shape()[1] {
                for ibasis in 0..residuals.shape()[2] {
                    residuals[[ielem, ivar, ibasis]] /= self.basis.mass_mat_diag[ibasis] * 0.5 * self.mesh.elements[ielem].jacob_det;
                }
            }
        }
    }
}