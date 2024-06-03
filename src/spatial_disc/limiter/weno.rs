use std::f64::EPSILON;

use ndarray::{s, Axis, Ix2, Ix3};
use ndarray::Array;
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
        for (in_cell_index, iedge) in element.iedges.indexed_iter() {
            let edge = &self.mesh.edges[*iedge];
            if edge.ielements[1] != -1 {
                let edge = &self.mesh.edges[*iedge];
                let ineighbour = element.ineighbours[in_cell_index] as usize;
                let (nx, ny, neighbour_in_cell_index) = {
                    if element.ivertices[in_cell_index] == edge.ivertices[0] {
                        (edge.normal[0], edge.normal[1], edge.in_cell_indices[0] as usize)
                    } else {
                        (-edge.normal[0], -edge.normal[1], edge.in_cell_indices[1] as usize)
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
                else {
                    continue;
                }
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
        for (ielem, element) in self.mesh.elements.indexed_iter() {
            if self.modified_kxrcf(solutions, ielem) {
                //dbg!(&ielem);
                let linear_weights = {
                    let mut weights = [0.997, 0.001, 0.001, 0.001];
                    let mut number_of_boundary_neighbours = 0;
                    for (i, iedge) in self.mesh.elements[ielem].iedges.indexed_iter() {
                        let edge_i = &self.mesh.edges[*iedge];
                        if edge_i.ielements[1] == -1 {
                            number_of_boundary_neighbours += 1;
                            weights[i] = 0.0;
                        }
                    }
                    match number_of_boundary_neighbours {
                        0 => weights[0] = 0.997,
                        1 => weights[1] = 0.998,
                        2 => weights[2] = 0.999,
                        _ => panic!("Invalid number of boundary neighbours!")
                    }
                    weights
                };
                let mut neighbour_area = 0.0;
                let tgt_pol = solutions.slice(s![ielem, .., ..]);
                let mut reconstructed_pol: Array<f64, Ix2> = Array::zeros((neq, nbasis));
                for (i, iedge) in self.mesh.elements[ielem].iedges.indexed_iter() {
                    let edge_i = &self.mesh.edges[*iedge];
                    if edge_i.ielements[1] != -1 {
                        let ineighbour = element.ineighbours[i] as usize;
                        neighbour_area += self.mesh.elements[ineighbour].jacob_det * 0.5;
                        let mut nonlinear_weights: Array<f64, Ix2> = Array::zeros((4, neq));
                        let mut proj_pols: Array<f64, Ix3> = Array::zeros((4, neq, nbasis));
                        let nx = self.mesh.edges[*iedge].normal[0];
                        let ny = self.mesh.edges[*iedge].normal[1];
                        let (lmatrix, rmatrix) = local_characteristics::compute_eigenmatrix(tgt_pol.slice(s![.., 0]), nx, ny, hcr);
                        proj_pols.slice_mut(s![0, .., ..]).assign(&lmatrix.dot(&tgt_pol));
                        for ivar in 0..neq {
                            let mut beta = 0.0;
                            for l in 1..=self.solver_param.order_of_polynomials {
                                for l_1 in 0..=l {
                                    let l_2 = l - l_1;
                                    let mut integral = 0.0;
                                    for (igp, cell_weight) in self.gauss_points.cell_weights.indexed_iter() {
                                        let mut value = 0.0;
                                        for ibasis in 0..nbasis {
                                            value += proj_pols[[0, ivar, ibasis]]
                                                * element.dphis_cell_gps[[igp, ibasis]].get(&(l_1, l_2)).unwrap();
                                        }
                                        integral += 0.5 * element.jacob_det * cell_weight * value.powf(2.0);
                                    }
                                    beta += element.jacob_det.powf((l - 1) as f64) * integral;                                            
                                }
                            }
                            nonlinear_weights[[0, ivar]] = linear_weights[0] / (beta + EPSILON).powf(2.0);
                        }
                        for (j, jedge) in element.iedges.indexed_iter() {
                            let edge_j = &self.mesh.edges[*jedge];
                            if edge_j.ielements[1] != -1 {
                                let jneighbour = element.ineighbours[j] as usize;
                                let average_modified_nb_pol = {
                                    let mut nb_pol = solutions.slice(s![jneighbour, .., ..]).to_owned();
                                    for ivar in 0..neq {
                                        nb_pol[[ivar, 0]] = tgt_pol[[ivar, 0]];
                                    }
                                    nb_pol
                                };
                                proj_pols.slice_mut(s![j + 1, .., ..]).assign(&lmatrix.dot(&average_modified_nb_pol));
                                for ivar in 0..neq {
                                    let mut beta = 0.0;
                                    for l in 1..=self.solver_param.order_of_polynomials {
                                        for l_1 in 0..=l {
                                            let l_2 = l - l_1;
                                            let mut integral = 0.0;
                                            for (igp, cell_weight) in self.gauss_points.cell_weights.indexed_iter() {
                                                let mut value = 0.0;
                                                for ibasis in 0..nbasis {
                                                    value += proj_pols[[j + 1, ivar, ibasis]]
                                                        * element.dphis_cell_gps[[igp, ibasis]].get(&(l_1, l_2)).unwrap();
                                                }
                                                integral += 0.5 * element.jacob_det * cell_weight * value.powf(2.0);
                                            }
                                            beta += element.jacob_det.powf((l - 1) as f64) * integral;                                            
                                        }
                                    }
                                    nonlinear_weights[[j + 1, ivar]] = linear_weights[j + 1] / (beta + EPSILON).powf(2.0);
                                }
                            }
                            else {
                                continue;
                            }
                        }
                        let mut weights: Array<f64, Ix2> = Array::zeros((4, neq));
                        let denominator = nonlinear_weights.sum_axis(Axis(0));
                        for ivar in 0..neq {
                            for j in 0..4 {
                                weights[[j, ivar]] = nonlinear_weights[[j, ivar]] / (denominator[ivar]);
                            }
                        }
                        let mut pol_new: Array<f64, Ix2> = Array::zeros((neq, nbasis));
                        for ivar in 0..neq {
                            for j in 0..4 {
                                for ibasis in 0..nbasis {
                                    pol_new[[ivar, ibasis]] += weights[[j, ivar]] * proj_pols[[j, ivar, ibasis]];
                                }
                            }
                        }
                        pol_new = rmatrix.dot(&pol_new);
                        reconstructed_pol = reconstructed_pol + pol_new * self.mesh.elements[ineighbour].jacob_det * 0.5;
                        /*
                        let ineighbour = self.mesh.elements[ielem].ineighbours[i] as usize;
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
                                                    * element.dphis_cell_gps[[igp, ibasis]].get(&(l_1, l_2)).unwrap()).powf(2.0);
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
                        */
                    }
                    else {
                        continue;
                    }
                }
                reconstructed_pol /= neighbour_area;
                //dbg!(&solutions.slice(s![ielem, .., ..]));
                //dbg!(&reconstructed_pol);
                solutions.slice_mut(s![ielem, .., ..]).assign(&reconstructed_pol);
            }
            else {
                continue;
            }
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
