use super::spatialdisc;
use super::mesh;
#[derive(Default)]
#[derive(Clone)]
pub struct ConsVar {
    pub density: f64,
    pub x_momentum: f64,
    pub y_momentum: f64,
    pub energy: f64,
}
pub struct Solver<'a> {
    pub mesh: mesh::Mesh<'a>,
    pub residuals: Vec<Vec<ConsVar>>,
    pub spatial_disc: spatialdisc::SpatialDisc,
    pub temperal_disc: Option<fn()>,
    
}
impl<'a> Solver<'a> {
    pub fn compute_residuals(&mut self, u: &Vec<f64>, time: f64) {
        self.mesh.
        self.spatial_disc.compute_fluxes(self.residuals: &mut Vec<Vec<ConsVar>>, mesh: &mesh::Mesh<'a>);
    }
    pub fn time_step(&mut self, u: &Vec<f64>, time: f64) {
        
    }
    pub fn solve(&mut self, u: &Vec<f64>, time: f64) {
        self.compute_residuals(u, time);
    }
}