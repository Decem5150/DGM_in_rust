use ndarray::{Ix1, Ix2, Ix3};
use ndarray::{Array, Axis};
use ndarray::s;
use crate::mesh::{BoundaryType, Edge, Patch};
use super::{local_characteristics, SpatialDisc};
impl<'a> SpatialDisc<'a> {
    pub fn apply_bc(&self, residuals: &mut Array<f64, Ix3>, solutions: &Array<f64, Ix3>) {
        /*
        for patch in self.mesh.patches.iter() {
            match patch.boundary_type {
                BoundaryType::Wall => self.wall(residuals, solutions, patch),
                BoundaryType::FarField => self.far_field(residuals, solutions, patch),
                //_ => panic!("Boundary type not implemented"),
            }
        }
        */
        let ngp = self.solver_param.number_of_cell_gp;
        let neq = self.solver_param.number_of_equations;
        let nbasis = self.solver_param.number_of_basis_functions;
        let cell_weights = &self.gauss_points.cell_weights;
        let edge_weights = &self.gauss_points.edge_weights;
        for &ielem in self.mesh.boundary_element_indices.iter() {
            let element = &self.mesh.elements[ielem];
            let mut global_lift_x: Array<f64, Ix2> = Array::zeros((neq, nbasis));
            let mut global_lift_y: Array<f64, Ix2> = Array::zeros((neq, nbasis));
            for (ilic, iedge) in self.mesh.elements[ielem].iedges.indexed_iter() {
                let edge = &self.mesh.edges[*iedge];
                if edge.ipatch == -1 {
                    let ilelem = ielem as usize;
                    let (iric, irelem, nx, ny) = {
                        if edge.in_cell_indices[0] as usize == ilic {
                            (edge.in_cell_indices[1] as usize, edge.ielements[1] as usize, edge.normal[0], edge.normal[1])
                        } else {
                            (edge.in_cell_indices[0] as usize, edge.ielements[0] as usize, -edge.normal[0], -edge.normal[1])
                        }
                    };
                    let left_dofs = solutions.slice(s![ilelem, .., ..]);
                    let right_dofs = solutions.slice(s![irelem, .., ..]);
                    for ivar in 0..neq {
                        for ibasis in 0..nbasis {
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
                        }
                    }
                }
                else {
                    let (nx, ny) = {
                        if element.ivertices[ilic] == edge.ivertices[0] {
                            (edge.normal[0], edge.normal[1])
                        } else {
                            (-edge.normal[0], -edge.normal[1])
                        }
                    };
                    let mut local_lift_x: Array<f64, Ix2> = Array::zeros((neq, nbasis));
                    let mut local_lift_y: Array<f64, Ix2> = Array::zeros((neq, nbasis));
                    let dofs = solutions.slice(s![ielem, .., ..]);
                    let boundary_type = self.mesh.patches[edge.ipatch as usize].boundary_type;
                    match boundary_type {
                        BoundaryType::Wall => {
                            for igp in 0..ngp {
                                let left_values: Array<f64, Ix1> = self.compute_boundary_edges_values(edge, solutions, igp);
                                let pressure = (self.flow_param.hcr - 1.0) * (left_values[3] - 0.5 * (left_values[1] * left_values[1] + left_values[2] * left_values[2]) / left_values[0]);
                                let boundary_energy = pressure / (self.flow_param.hcr - 1.0);
                                let x_momentum_jump_x = left_values[1] * nx;
                                let x_momentum_jump_y = left_values[1] * ny;
                                let y_momentum_jump_x = left_values[2] * nx;
                                let y_momentum_jump_y = left_values[2] * ny;
                                let energy_jump_x = (left_values[3] - boundary_energy) * nx;
                                let energy_jump_y = (left_values[3] - boundary_energy) * ny;
                                for ibasis in 0..nbasis {
                                    local_lift_x[[1, ibasis]] -= 0.5 * x_momentum_jump_x * self.basis.phis_cell_gps[[igp, ibasis]] * cell_weights[igp] * 0.5 * element.jacob_det;
                                    local_lift_x[[2, ibasis]] -= 0.5 * x_momentum_jump_y * self.basis.phis_cell_gps[[igp, ibasis]] * cell_weights[igp] * 0.5 * element.jacob_det;
                                    local_lift_x[[3, ibasis]] -= 0.5 * energy_jump_x * self.basis.phis_cell_gps[[igp, ibasis]] * cell_weights[igp] * 0.5 * element.jacob_det;
                                    local_lift_y[[1, ibasis]] -= 0.5 * y_momentum_jump_x * self.basis.phis_cell_gps[[igp, ibasis]] * cell_weights[igp] * 0.5 * element.jacob_det;
                                    local_lift_y[[2, ibasis]] -= 0.5 * y_momentum_jump_y * self.basis.phis_cell_gps[[igp, ibasis]] * cell_weights[igp] * 0.5 * element.jacob_det;
                                    local_lift_y[[3, ibasis]] -= 0.5 * energy_jump_y * self.basis.phis_cell_gps[[igp, ibasis]] * cell_weights[igp] * 0.5 * element.jacob_det;
                                }
                            }
                            global_lift_x = global_lift_x + local_lift_x;
                            global_lift_y = global_lift_y + local_lift_y;
                            for ivar in 1..neq {
                                for ibasis in 0..nbasis {
                                    local_lift_x[[ivar, ibasis]] /= self.basis.mass_mat_diag[ibasis] * 0.5 * element.jacob_det;
                                    local_lift_y[[ivar, ibasis]] /= self.basis.mass_mat_diag[ibasis] * 0.5 * element.jacob_det;
                                }
                            }
                            for igp in 0..ngp {
                                let left_values: Array<f64, Ix1> = self.compute_boundary_edges_values(edge, solutions, igp);
                                let pressure = (self.flow_param.hcr - 1.0) * (left_values[3] - 0.5 * (left_values[1] * left_values[1] + left_values[2] * left_values[2]) / left_values[0]);
                                let boundary_inviscid_flux = [
                                    0.0,
                                    pressure * nx,
                                    pressure * ny,
                                    0.0,
                                ];
                                let mut boundary_viscous_flux: Array<f64, Ix1> = Array::zeros(neq);
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
                                            let dphi_dx = element.dphis_edge_gps[[ilic, igp, ibasis]].get(&(1, 0)).unwrap();
                                            let dphi_dy = element.dphis_edge_gps[[ilic, igp, ibasis]].get(&(0, 1)).unwrap();
                                            dsol_dx += dofs[[ivar, ibasis]] * dphi_dx;
                                            dsol_dy += dofs[[ivar, ibasis]] * dphi_dy;
                                            lift_x += local_lift_x[[ivar, ibasis]] * self.basis.phis_edge_gps[[ilic, igp, ibasis]];
                                            lift_y += local_lift_y[[ivar, ibasis]] * self.basis.phis_edge_gps[[ilic, igp, ibasis]];
                                        }
                                        fv += g11[[iflux, ivar]] * (dsol_dx + lift_x) + g12[[iflux, ivar]] * (dsol_dy + lift_y);
                                        gv += g21[[iflux, ivar]] * (dsol_dx + lift_x) + g22[[iflux, ivar]] * (dsol_dy + lift_y);
                                    }
                                }
                                boundary_viscous_flux = nx * &fv + ny * &gv;
                                for ivar in 0..residuals.shape()[1] {
                                    for ibasis in 0..residuals.shape()[2] {
                                        residuals[[ielem, ivar, ibasis]] -= edge_weights[igp] * (boundary_inviscid_flux[ivar] + boundary_viscous_flux[ivar]) * self.basis.phis_edge_gps[[ilic, igp, ibasis]] * 0.5 * edge.jacob_det;
                                    }
                                }
                            }
                        }
                        BoundaryType::FarField => {
                            let boundary_quantity = self.mesh.patches[edge.ipatch as usize].boundary_quantity.as_ref().unwrap();
                            let boundary_consvar = Array::from_vec(vec![
                                boundary_quantity.rho, 
                                boundary_quantity.rho * boundary_quantity.u, 
                                boundary_quantity.rho * boundary_quantity.v, 
                                boundary_quantity.p / (self.flow_param.hcr - 1.0) + 0.5 * boundary_quantity.rho * (boundary_quantity.u.powi(2) + boundary_quantity.v.powi(2))
                                ]);
                            for igp in 0..ngp {
                                let left_values: Array<f64, Ix1> = self.compute_boundary_edges_values(edge, solutions, igp);
                                for ibasis in 0..nbasis {
                                    let u = left_values[1] / left_values[0];
                                    let v = left_values[2] / left_values[0];
                                    let p = (hcr - 1.0) * (left_values[3] - 0.5 * (left_values[1] * left_values[1] + left_values[2] * left_values[2]) / left_values[0]);
                                    let c = (hcr * p / left_values[0]).sqrt();
                                    let vn = u * nx + v * ny;
                                    let (left_elem_lmatrix, left_elem_rmatrix) = local_characteristics::compute_eigenmatrix(left_values, nx, ny, hcr);
                                    let left_elem_eigenvalues = [vn - c, vn, vn, vn + c];
                                    let alpha: Array<f64, Ix2> = left_elem_lmatrix.dot(&left_values).insert_axis(Axis(1));
                                    let beta: Array<f64, Ix2> = left_elem_lmatrix.dot(&boundary_consvar.clone().insert_axis(Axis(1)));
                                    let mut q_new: Array<f64, Ix1> = Array::zeros(neq);
                                }
                            }
                        }
                    }

                }
            }
        }
    }
    fn wall(&self, residuals: &mut Array<f64, Ix3>, solutions: &Array<f64, Ix3>, patch: &Patch) {
        let nbasis = self.solver_param.number_of_basis_functions;
        let ngp = self.solver_param.number_of_edge_gp;
        let hcr = self.flow_param.hcr;
        for &iedge in patch.iedges.iter() {
            let ielem = self.mesh.boundary_edges[iedge].ielement;
            let in_cell_index = self.mesh.boundary_edges[iedge].in_cell_index;
            let (nx, ny) = {
                if self.mesh.elements[ielem].ivertices[in_cell_index] == self.mesh.boundary_edges[iedge].ivertices[0] {
                    (self.mesh.boundary_edges[iedge].normal[0], self.mesh.boundary_edges[iedge].normal[1])
                } else {
                    (-self.mesh.boundary_edges[iedge].normal[0], -self.mesh.boundary_edges[iedge].normal[1])
                }
            };
            let left_value: Array<f64, Ix2> = self.compute_boundary_edges_values(&self.mesh.boundary_edges[iedge], solutions);
            for igp in 0..ngp {
                let pressure = (hcr - 1.0) * (left_value[[igp, 3]] - 0.5 * (left_value[[igp, 1]] * left_value[[igp, 1]] + left_value[[igp, 2]] * left_value[[igp, 2]]) / left_value[[igp, 0]]);
                for ibasis in 0..nbasis {
                    residuals[[ielem, 1, ibasis]] -= pressure * self.gauss_points.edge_weights[igp] * self.basis.phis_edge_gps[[in_cell_index, igp, ibasis]] * nx * 0.5 * self.mesh.boundary_edges[iedge].jacob_det;
                    residuals[[ielem, 2, ibasis]] -= pressure * self.gauss_points.edge_weights[igp] * self.basis.phis_edge_gps[[in_cell_index, igp, ibasis]] * ny * 0.5 * self.mesh.boundary_edges[iedge].jacob_det;
                }
            }
        }
    }
    fn far_field(&self, residuals: &mut Array<f64, Ix3>, solutions: &Array<f64, Ix3>, patch: &Patch) {
        let nbasis = self.solver_param.number_of_basis_functions;
        let ngp = self.solver_param.number_of_edge_gp;
        let neq = self.solver_param.number_of_equations;
        let hcr = self.flow_param.hcr;
        let boundary_quantity = patch.boundary_quantity.as_ref().unwrap();
        let boundary_consvar = Array::from_vec(vec![
            boundary_quantity.rho, 
            boundary_quantity.rho * boundary_quantity.u, 
            boundary_quantity.rho * boundary_quantity.v, 
            boundary_quantity.p / (hcr - 1.0) + 0.5 * boundary_quantity.rho * (boundary_quantity.u.powi(2) + boundary_quantity.v.powi(2))
            ]);
        for &iedge in patch.iedges.iter() {
            let ielem = self.mesh.boundary_edges[iedge].ielement;
            let in_cell_index = self.mesh.boundary_edges[iedge].in_cell_index;
            let (nx, ny) = {
                if self.mesh.elements[ielem].ivertices[in_cell_index] == self.mesh.boundary_edges[iedge].ivertices[0] {
                    (self.mesh.boundary_edges[iedge].normal[0], self.mesh.boundary_edges[iedge].normal[1])
                } else {
                    (-self.mesh.boundary_edges[iedge].normal[0], -self.mesh.boundary_edges[iedge].normal[1])
                }
            };
            let left_values_gps: Array<f64, Ix2> = self.compute_boundary_edges_values(&self.mesh.boundary_edges[iedge], solutions);
            for igp in 0..ngp {
                for ibasis in 0..nbasis {
                    let u = left_values_gps[[igp, 1]] / left_values_gps[[igp, 0]];
                    let v = left_values_gps[[igp, 2]] / left_values_gps[[igp, 0]];
                    let p = (hcr - 1.0) * (left_values_gps[[igp, 3]] - 0.5 * (left_values_gps[[igp, 1]] * left_values_gps[[igp, 1]] + left_values_gps[[igp, 2]] * left_values_gps[[igp, 2]]) / left_values_gps[[igp, 0]]);
                    let c = (hcr * p / left_values_gps[[igp, 0]]).sqrt();
                    let vn = u * nx + v * ny;
                    let (left_elem_lmatrix, left_elem_rmatrix) = local_characteristics::compute_eigenmatrix(left_values_gps.slice(s![igp, ..]), nx, ny, hcr);
                    let left_elem_eigenvalues = [vn - c, vn, vn, vn + c];
                    let alpha: Array<f64, Ix2> = left_elem_lmatrix.dot(&left_values_gps.slice(s![igp, ..]).insert_axis(Axis(1)));
                    let beta: Array<f64, Ix2> = left_elem_lmatrix.dot(&boundary_consvar.clone().insert_axis(Axis(1)));
                    let mut q_new: Array<f64, Ix1> = Array::zeros(neq);
                    for i in 0..neq {
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
                    for ivar in 0..neq {
                        residuals[[ielem, ivar, ibasis]] -= flux[ivar] * self.gauss_points.edge_weights[igp] * self.basis.phis_edge_gps[[in_cell_index, igp, ibasis]] * 0.5 * self.mesh.boundary_edges[iedge].jacob_det;
                    }
                }
            }
        }
    }
    fn compute_boundary_edges_values(&self, edge: &Edge, solutions: &Array<f64, Ix3>, igp: usize) -> Array<f64, Ix1> {
        let nbasis = self.solver_param.number_of_basis_functions;
        let ngp = self.solver_param.number_of_edge_gp;
        let neq = self.solver_param.number_of_equations;
        let mut edges_values = Array::zeros(neq);
        let ielem = edge.ielements[0];
        let in_cell_index = edge.in_cell_indices[0] as usize;
        let internal_sol = solutions.slice(s![ielem, .., ..]);
        for ivar in 0..neq {
            for ibasis in 0..nbasis {
                edges_values[ivar] += internal_sol[[ivar, ibasis]] * self.basis.phis_edge_gps[[in_cell_index, igp, ibasis]];
            }
        }
        edges_values
    }
}
