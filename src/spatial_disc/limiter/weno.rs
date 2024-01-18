use ndarray::{Ix1, Ix2, Ix3, s};
use ndarray::{Array, ArrayView, ArrayViewMut};
use crate::mesh::{Element, Edge, NormalDirection};
use crate::spatial_disc::{local_characteristics, SpatialDisc};
impl SpatialDisc<'_> {
    pub fn modified_kxrcf(&self, solutions: &Array<f64, Ix3>) -> Vec<usize> {
        let nelem = self.solver_param.number_of_elements;
        let nbasis = self.solver_param.number_of_basis_functions;
        let ngp = self.solver_param.number_of_edge_gp;
        let neq = self.solver_param.number_of_equations;
        let total_energy_jump = Array::zeros(nelem);
        let mut indicator_density = 0.0;
        // let mut indicator_total_enthalpy = 0.0;
        let mut surface = 0.0;
        let mut min_density = {
            let mut rho = 0.0;
            for ibasis in 0..nbasis {
                rho += solutions[[ielem, 0, ibasis]] * self.basis.phis_edge_gps[[0, 0, ibasis]]; 
            }
            rho
        };
        for &iedge in element.iedges.iter() {
            let edge = &self.mesh.edges[iedge];
            let ilelem = edge.in_cell_index[0];
            let irelem = edge.in_cell_index[1];
            let ineighbour = element.ineighbours[ilelem].unwrap();
            for igp in 0..ngp {
                let mut rho = 0.0;
                let mut rho_u = 0.0;
                let mut rho_v = 0.0;
                for ibasis in 0..nbasis {
                    rho += solutions[[ielem, 0, ibasis]] * self.basis.phis_edge_gps[[ilelem, igp, ibasis]];
                    rho_u += solutions[[ielem, 1, ibasis]] * self.basis.phis_edge_gps[[ilelem, igp, ibasis]];
                    rho_v += solutions[[ielem, 2, ibasis]] * self.basis.phis_edge_gps[[ilelem, igp, ibasis]];
                }
                if rho < min_density {
                    min_density = rho;
                }
                let (nx, ny) = match element.normal_directions[iedge] {
                    NormalDirection::Inward => (-edge.normal[0], -edge.normal[1]),
                    NormalDirection::Outward => (edge.normal[0], edge.normal[1]),
                };
                let rho_vn = rho_u * nx + rho_v * ny;
                if rho_vn < 0.0 {
                    let mut rho_int = 0.0;
                    let mut rho_ext = 0.0;
                    for ibasis in 0..nbasis {
                        rho_int += solutions[[ielem, 0, ibasis]] * self.basis.phis_edge_gps[[ilelem, igp, ibasis]];
                        rho_ext += solutions[[ineighbour, 0, ibasis]] * self.basis.phis_edge_gps[[irelem, ngp - 1 - igp, ibasis]];
                    }
                    indicator_density += rho_int - rho_ext;
                    surface += edge.jacob_det;
                }
            }
        }
        let h = element.circumradius;
        indicator_density = indicator_density.abs() / (h.powf(1.5) * surface * min_density);
        indicator_density > 1.0
    }
    
    pub fn weno(&self, solutions: &mut Array<f64, Ix3>) {
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
            // let proj_pol = Array::zeros((3, 4, neq, nbasis));
            if self.modified_kxrcfkxrcf(solutions, &self.mesh.elements[ielem]) {
                for (i, &iedge) in self.mesh.elements[ielem].iedges.indexed_iter() {
                    let ineighbour = self.mesh.elements[ielem].ineighbours[i];
                    let nx = self.mesh.edges[iedge].normal[0];
                    let ny = self.mesh.edges[iedge].normal[1];
                    let (lmatrix, rmatrix)= local_characteristics::compute_eigenmatrix(orig_pols[0], nx, ny, hcr);
                    let proj_tgt_pol = lmatrix.dot(&orig_pols[0]);
                    let proj_neigh_pol = lmatrix.dot(&orig_pols[i + 1]);
                    for ivar in 0..neq {
                        let beta_tgt = ()self.mesh.elements[ielem].jacob_det
                    }

                }
                
            }
            else {}
        }
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
