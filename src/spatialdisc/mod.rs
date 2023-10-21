pub mod flux;
pub mod quadrature;
pub mod gauss_point;
pub mod basis_function;
pub struct SpatialDisc<'a>{
    pub cell_gauss_points: Vec<[f64; 2]>,
    pub cell_gauss_weights: Vec<f64>,
    pub edge_gauss_points: Vec<f64>,
    pub edge_gauss_weights: Vec<f64>,
    pub compute_inviscid_flux: fn(left_value: ConsVar, right_value: ConsVar, flux: &'a mut InvisFlux, nx: f64, ny: f64), 
    pub compute_viscous_flux: fn(left_value: ConsVar, right_value: ConsVar, flux: &'a mut VisFlux, nx: f64, ny: f64),
    pub integrate_cell: fn(area: &'a f64, quad: gauss_point::Quadrature2D, element: &'a Element) -> (Vec<f64>, Vec<f64>, Vec<f64>, Vec<f64>),
    pub integrate_edge: fn(length: &'a f64, quad: gauss_point::Quadrature1D, edge: &'a Edge) -> (Vec<f64>, Vec<f64>, Vec<f64>, Vec<f64>),
}
impl<'a> SpatialDisc<'a> {
    pub fn integrate_over_cell(&self, area: &'a f64, weights: &'a Vec<f64>, points: &'a Vec<[f64; 2]>, element: &'a Element) -> (Vec<f64>, Vec<f64>, Vec<f64>, Vec<f64>) {
        (self.integrate_cell)();
    }
    pub fn compute_fluxes() {

    }
}