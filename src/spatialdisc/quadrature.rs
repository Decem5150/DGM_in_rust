use super::flux;
use super::basis_function;
use super::gauss_point;
use crate::mesh;
use crate::solver;
pub fn quadrature_triangle<'a>(res: &mut Vec<solver::ConsVar>, basis: &basis_function::DubinerBasis, quad: &gauss_point::Quadrature2D, element: &'a mesh::Element) {
    let mut gauss_values = vec![solver::ConsVar::default(); quad.points.len()];
    // Compute solutions at gauss points.
    for i in 0..quad.points.len() {
        for j in 0..basis.dof {
            gauss_values[i].density += basis.phi_gauss[i][j] * element.solution[j].density;
            gauss_values[i].x_momentum += basis.phi_gauss[i][j] * element.solution[j].x_momentum;
            gauss_values[i].y_momentum += basis.phi_gauss[i][j] * element.solution[j].y_momentum;
            gauss_values[i].energy += basis.phi_gauss[i][j] * element.solution[j].energy;
        }
    }
    // Compute cell integration term and add to residuals.
    for i in 0..quad.points.len() {
        let (f1, f2, f3, f4, g1, g2, g3, g4) = flux::flux(gauss_values[i].density, gauss_values[i].x_momentum, gauss_values[i].y_momentum, gauss_values[i].energy);
        for j in 0..basis.dof {
            res[j].density += quad.weights[i] * element.dphi_dx[i][j] * f1 + quad.weights[i] * element.dphi_dy[i][j] * g1;
            res[j].x_momentum += quad.weights[i] * element.dphi_dx[i][j] * f2 + quad.weights[i] * element.dphi_dy[i][j] * g2;
            res[j].y_momentum += quad.weights[i] * element.dphi_dx[i][j] * f3 + quad.weights[i] * element.dphi_dy[i][j] * g3;
            res[j].energy += quad.weights[i] * element.dphi_dx[i][j] * f4 + quad.weights[i] * element.dphi_dy[i][j] * g4;
        }
    }
}