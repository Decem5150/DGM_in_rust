use ndarray::{Ix1, Ix2, Ix3};
use ndarray::ArrayViewMut;
use crate::mesh::{Element, Edge, NormalDirection, EdgeType};
use super::Limiter;
impl<'a> Limiter<'a> {
    pub fn kxrcf(&self, solutions: ArrayViewMut<f64, Ix3>, element: &Element<'a>) -> bool {
        let nelem = self.solver_param.number_of_elements;
        let nbasis = self.solver_param.number_of_basis_functions;
        let ngp = self.solver_param.number_of_edge_gp;
        let neq = self.solver_param.number_of_equations;
        let ielem = element.solution_index;
        let mut indicator_density = 0.0;
        let mut indicator_total_enthalpy = 0.0;
        let mut surface = 0.0;
        for (iedge, edge_type) in element.edges.indexed_iter() {
            if let EdgeType::Internal(edge) = edge_type {
                let ilelem = edge.ind_in_left_elem;
                let irelem = edge.ind_in_right_elem;
                for igp in 0..ngp {
                    let mut rho_u = 0.0;
                    let mut rho_v = 0.0;
                    for ibasis in 0..nbasis {
                        rho_u += solutions[[ielem, 1, ibasis]] * self.basis.phis_edge_gps[[ilelem, igp, ibasis]];
                        rho_v += solutions[[ielem, 2, ibasis]] * self.basis.phis_edge_gps[[ilelem, igp, ibasis]];
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
                            rho_ext += solutions[[element.neighbours[iedge].solution_index, 0, ibasis]] * self.basis.phis_edge_gps[[irelem, ngp - 1 - igp, ibasis]];
                        }
                        indicator_density += rho_int - rho_ext;
                        surface += edge.jacob_det;
                    }
                    else {}
                }
            }
            else {}
        }
        let h = element.circumradius;
        let mut density_l2_norm = 0.0;
        for igp in 0..ngp {
            let mut rho = 0.0;
            for ibasis in 0..nbasis {
                rho += solutions[[ielem, 0, ibasis]] * self.basis.phis_edge_gps[[ielem, igp, ibasis]]; 
            }
            density_l2_norm += rho.powf(2.0) * self.gauss_points.cell_weights[igp] * 0.5 * element.jacob_det;
        }
        density_l2_norm = density_l2_norm.sqrt();
        indicator_density = indicator_density.abs() / (h.powf(nbasis as f64 * 0.5) * surface * density_l2_norm);
        if indicator_density > 1.0 {
            true
        }
        else {
            false
        }
    }
    pub fn hweno(&self, solutions: ArrayViewMut<f64, Ix3>) {
        for element in self.mesh.internal_elements.iter() {
            if self.kxrcf(solutions, element) {

            }
            else {}
        }
    }
}
