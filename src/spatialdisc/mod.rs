pub mod flux;
pub mod quadrature;
pub mod gauss_point;
pub mod basis_function;
pub use crate::mesh::Element;
pub use crate::mesh::Mesh;
pub struct SpatialDisc{
    pub cell_gauss_points: Vec<[f64; 2]>,
    pub cell_gauss_weights: Vec<f64>,
    pub edge_gauss_points: Vec<f64>,
    pub edge_gauss_weights: Vec<f64>,
    pub compute_inviscid_flux: Box<dyn flux::InvisFluxScheme>,
    // pub compute_viscous_flux: Box<dyn flux::VisFluxScheme>,
}
impl<'a> SpatialDisc {
    pub fn integrate_over_cell(&self, element: &'a Element) -> (Vec<f64>, Vec<f64>, Vec<f64>, Vec<f64>) {
    }
    pub fn compute_fluxes(&self, mesh: &Mesh) {
        for edge in mesh.edges.iter() {
            let left_element = edge.elements[0];
            let right_element = edge.elements[1];
            let left_value = left_element.solution;
            let right_value = right_element.solution;
            let nx = edge.normal[0];
            let ny = edge.normal[1];
            self.compute_inviscid_flux.compute(left_value, right_value, edge.invis_flux, nx, ny);
        }
    }
}