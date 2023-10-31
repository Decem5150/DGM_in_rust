pub mod flux;
pub mod quadrature;
pub mod gauss_point;
pub mod basis_function;
pub mod boundary;
pub use crate::mesh::Element;
pub use crate::mesh::Mesh;
pub struct SpatialDisc<'a> {
    pub gauss_point: gauss_point::GaussPoints,
    pub inviscid_flux: Box<dyn flux::InvisFluxScheme<'a>>,
    // pub viscous_flux: Box<dyn flux::VisFluxScheme>,
}
impl<'a> SpatialDisc<'a> {
    pub fn integrate_over_cell(&self, element: &'a Element) -> (Vec<f64>, Vec<f64>, Vec<f64>, Vec<f64>) {
    }
    pub fn compute_fluxes(&self, mesh: &Mesh) {
        for edge in mesh.edges.iter() {
            let left_element = edge.elements[0];
            let right_element = edge.elements[1];
            let left_value = left_element.solution;
            let right_value = right_element.solution;
            self.inviscid_flux.compute(left_value, right_value, edge.invis_flux, edge.normal);
        }
    }
}