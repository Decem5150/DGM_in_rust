use std::f64::EPSILON;

use ndarray::{s, Ix2, Ix3};
use ndarray::{Array, ArrayView, ArrayViewMut};
use crate::mesh::{BoundaryType, Edge, EdgeTypeAndIndex, Element};
use crate::spatial_disc::{local_characteristics, SpatialDisc};
impl SpatialDisc<'_> {
    pub fn modified_kxrcf(&self, solutions: &mut Array<f64, Ix3>, ielem: usize) -> bool {
        let nbasis = self.solver_param.number_of_basis_functions;
        let ngp = self.solver_param.number_of_edge_gp;
        let element = &self.mesh.elements[ielem];
        let mut density_jump = 0.0;
        let mut total_energy_jump = 0.0;
        let mut surface = 0.0;
        let mut min_density = 1.0e10;
        let mut min_total_energy = 1.0e10;
        for (in_cell_index, edge_type_and_index) in element.iedges.indexed_iter() {
            match edge_type_and_index {
                EdgeTypeAndIndex::Internal(iedge) => {
                    let edge = &self.mesh.edges[*iedge];
                    let ineighbour = element.ineighbours[in_cell_index].unwrap();
                    let (nx, ny, neighbour_in_cell_index) = {
                        if element.ivertices[in_cell_index] == edge.ivertices[0] {
                            (edge.normal[0], edge.normal[1], edge.in_cell_indices[0])
                        } else {
                            (-edge.normal[0], -edge.normal[1], edge.in_cell_indices[1])
                        }
                    };
                    for igp in 0..ngp {
                        let mut rho = 0.0;
                        let mut rho_e = 0.0;
                        for ibasis in 0..nbasis {
                            rho += solutions[[ielem, 0, ibasis]] * self.basis.phis_edge_gps[[in_cell_index, igp, ibasis]];
                            rho_e += solutions[[ielem, 3, ibasis]] * self.basis.phis_edge_gps[[in_cell_index, igp, ibasis]];
                        }
                        if rho < min_density {
                            min_density = rho;
                        }
                        if rho_e < min_total_energy {
                            min_total_energy = rho_e;
                        } 
                    }
                    let center_rho_velocity = {
                        let center_gp_index = ngp / 2;
                        let mut rho_u = 0.0;
                        let mut rho_v = 0.0;
                        for ibasis in 0..nbasis {
                            rho_u += solutions[[ielem, 1, ibasis]] * self.basis.phis_edge_gps[[in_cell_index, center_gp_index, ibasis]];
                            rho_v += solutions[[ielem, 2, ibasis]] * self.basis.phis_edge_gps[[in_cell_index, center_gp_index, ibasis]];
                        }
                        rho_u * nx + rho_v * ny
                    };
                    if center_rho_velocity < 0.0 {
                        let mut rho_jump = 0.0;
                        let mut rho_e_jump = 0.0;  
                        for igp in 0..ngp {
                            let mut int_rho_gp = 0.0;
                            let mut int_rho_e_gp = 0.0;
                            let mut ext_rho_gp = 0.0;
                            let mut ext_rho_e_gp = 0.0;
                            for ibasis in 0..nbasis {
                                int_rho_gp += solutions[[ielem, 0, ibasis]] * self.basis.phis_edge_gps[[in_cell_index, igp, ibasis]];
                                int_rho_e_gp += solutions[[ielem, 3, ibasis]] * self.basis.phis_edge_gps[[in_cell_index, igp, ibasis]];
                                ext_rho_gp += solutions[[ineighbour, 0, ibasis]] * self.basis.phis_edge_gps[[neighbour_in_cell_index, ngp - 1 - igp, ibasis]];
                                ext_rho_e_gp += solutions[[ineighbour, 3, ibasis]] * self.basis.phis_edge_gps[[neighbour_in_cell_index, ngp - 1 - igp, ibasis]];
                            }
                            rho_jump += 0.5 * edge.jacob_det * self.gauss_points.edge_weights[igp] * (int_rho_gp - ext_rho_gp);
                            rho_e_jump += 0.5 * edge.jacob_det * self.gauss_points.edge_weights[igp] * (int_rho_e_gp - ext_rho_e_gp);
                        }
                        density_jump += rho_jump;
                        total_energy_jump += rho_e_jump;
                        surface += edge.jacob_det;
                    }
                }
                EdgeTypeAndIndex::Boundary(_) => {},
            }
        }
        let h = element.circumradius;
        //dbg!(&surface);
        //dbg!(&density_jump.abs());
        let indicator_density = density_jump.abs() / (h.powf(1.5) * surface * min_density);
        let indicator_total_energy = total_energy_jump.abs() / (h.powf(1.5) * surface * min_total_energy);
        //dbg!(&indicator_density);
        (indicator_density >= 1.0) || (indicator_total_energy >= 1.0)
    }
    pub fn weno(&self, solutions: &mut Array<f64, Ix3>) {
        let neq = self.solver_param.number_of_equations;
        let nbasis = self.solver_param.number_of_basis_functions;
        let hcr = self.flow_param.hcr;
        let linear_weights = [0.998, 0.002];
        for (ielem, element) in self.mesh.elements.indexed_iter() {
            let mut neighbour_area = 0.0;
            if self.modified_kxrcf(solutions, ielem) {
                //dbg!(&ielem);
                let mut reconstructed_pol: Array<f64, Ix2> = Array::zeros((neq, nbasis));
                for (i, edge_type_and_index) in self.mesh.elements[ielem].iedges.indexed_iter() {
                    match edge_type_and_index {
                        EdgeTypeAndIndex::Internal(iedge) => {
                            let ineighbour = self.mesh.elements[ielem].ineighbours[i].unwrap();
                            neighbour_area += self.mesh.elements[ineighbour].jacob_det * 0.5;
                            let tgt_pol = solutions.slice(s![ielem, .., ..]);
                            let average_modified_nb_pol = {
                                let mut nb_pol = solutions.slice(s![ineighbour, .., ..]).to_owned();
                                for ivar in 0..neq {
                                    nb_pol[[ivar, 0]] = tgt_pol[[ivar, 0]];
                                }
                                nb_pol
                            };
                            let nx = self.mesh.edges[*iedge].normal[0];
                            let ny = self.mesh.edges[*iedge].normal[1];
                            let (lmatrix, rmatrix) = local_characteristics::compute_eigenmatrix(solutions.slice(s![ielem, .., 0]), nx, ny, hcr);
                            let proj_tgt_pol: Array<f64, Ix2> = lmatrix.dot(&tgt_pol);
                            let proj_nb_pol: Array<f64, Ix2> = lmatrix.dot(&average_modified_nb_pol);
                            let proj_pols: [Array<f64, Ix2>; 2] = [proj_tgt_pol, proj_nb_pol];
                            let mut nonlinear_weights = [0.0; 2];
                            {
                                let mut beta = [0.0; 2];
                                let k = self.solver_param.order_of_polynomials;
                                for (i, proj_pol) in proj_pols.iter().enumerate() {
                                    for l in 1..=k {
                                        for l_1 in 0..=l {
                                            let l_2 = l - l_1;
                                            let mut integral = 0.0;
                                            for (igp, cell_weight) in self.gauss_points.cell_weights.indexed_iter() {
                                                for ((_, ibasis), proj_pol_value) in proj_pol.indexed_iter() {
                                                    integral += 0.5 * element.jacob_det * cell_weight 
                                                        * (1.0 / Self::factorial(l) * proj_pol_value
                                                        * element.derivatives[[igp, ibasis]].get(&(l_1, l_2)).unwrap()).powf(2.0);
                                                }
                                            }
                                            /*
                                            for igp in 0..ngp {
                                                for ivar in 0..neq {
                                                    for ibasis in 0..nbasis {
                                                        integral += 0.5 * element.jacob_det * self.gauss_points.cell_weights[igp] 
                                                            * (1.0 / self.factorial(l) * proj_pol[[ivar, ibasis]] 
                                                            * element.derivatives[[igp, ibasis]].get(&(l_1, l_2)).unwrap()).powf(2.0);
                                                    }
                                                }
                                            }
                                            */
                                            beta[i] += element.jacob_det.powf((l - 1) as f64) * integral;                                            
                                        }
                                    }
                                    nonlinear_weights[i] = linear_weights[i] / (beta[i] + EPSILON).powf(2.0);
                                }
                            }
                            let mut weights = [0.0; 2];
                            for i in 0..2 {
                                weights[i] = nonlinear_weights[i] / (nonlinear_weights[0] + nonlinear_weights[1]);
                            }
                            let pol_new: Array<f64, Ix2> = weights[0] * &proj_pols[0] + weights[1] * &proj_pols[1];
                            reconstructed_pol = reconstructed_pol + rmatrix.dot(&pol_new) * self.mesh.elements[ineighbour].jacob_det * 0.5;
                        },
                        EdgeTypeAndIndex::Boundary(_) => {},
                    }
                }
                reconstructed_pol /= neighbour_area;
                //dbg!(&solutions.slice(s![ielem, .., ..]));
                //dbg!(&reconstructed_pol);
                solutions.slice_mut(s![ielem, .., ..]).assign(&reconstructed_pol);
            }
            else {}
        }
    }
    fn factorial(n: usize) -> f64 {
        let mut result = 1;
        for i in 1..=n {
            result *= i;
        }
        result as f64
    }
    /*
    pub fn hweno(&self, solutions: ArrayViewMut<f64, Ix3>) {
        let neq = self.solver_param.number_of_equations;
        let nbasis = self.solver_param.number_of_basis_functions;
        let hcr = self.flow_param.hcr;
        for &ielem in self.mesh.internal_elements.iter() {
            let orig_pols = [
                solutions.slice(s![ielem, .., ..]),
                solutions.slice(s![self.mesh.elements[ielem].ineighbours[0], .., ..]), 
                solutions.slice(s![self.mesh.elements[ielem].ineighbours[1], .., ..]), 
                solutions.slice(s![self.mesh.elements[ielem].ineighbours[2], .., ..])
                ];
            let proj_pol = Array::zeros((3, 4, neq, nbasis));
            if self.kxrcf(solutions, &self.mesh.elements[ielem]) {
        for (i, &iedge) in self.mesh.elements[ielem].iedges.indexed_iter() {
            let ineighbour = self.mesh.elements[ielem].ineighbours[iedge];
            let nx = self.mesh.edges[iedge].normal[0];
            let ny = self.mesh.edges[iedge].normal[1];
            let (lmatrix, rmatrix)= local_characteristics::compute_eigenmatrix(orig_pols[0], nx, ny, hcr);
            for l in 0..4 {
                proj_pol.slice_mut(s![i, l, .., ..]) = &lmatrix.dot(&orig_pols[l]);
            }
        }
        for (i, &iedge) in self.mesh.elements[ielem].iedges.indexed_iter() {
            for l in 1..=3 {
                for ivar in 0..neq {
                    let iuseful = {
                        let indices: Vec<usize> = vec![1, 2, 3];
                        indices = indices.remove(l - 1);
                        let dist: Vec<f64> = vec![
                            (proj_pol[[i, indices[0], ivar, 0]] - proj_pol[[i, 0, ivar, 0]]).abs(), 
                            (proj_pol[[i, indices[1], ivar, 0]] - proj_pol[[i, 0, ivar, 0]]).abs()
                            ];
                        if dist[0] < dist[1] {
                            indices[0]
                        }
                        else {
                            indices[1]
                        }
                    };
                }
            }
        }
    }
    */
}
