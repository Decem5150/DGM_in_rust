pub mod rk3;
pub use crate::solver::{SolverParameters, Solver, SolCoeff};
pub use crate::spatial_disc::SpatialDisc;
pub use crate::mesh::{Mesh, Element};
pub struct TemperalDisc<'a> {
    pub solver: &'a Solver<'a>,
    pub spatial_disc: &'a SpatialDisc<'a>,
    pub mesh: &'a Mesh<'a>,
    pub solver_param: &'a SolverParameters,
}
trait TimeScheme {
    fn step(&self, residuals: &mut Vec<SolCoeff>, time_step: f64);
}