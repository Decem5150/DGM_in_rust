pub mod rk3;
use ndarray::Array;
use ndarray::{Ix1, Ix2, Ix3};
use ndarray::array;
use ndarray::{ArrayView, ArrayViewMut};
pub use crate::solver::{SolverParameters, Solver, SolCoeff};
pub use crate::spatial_disc::SpatialDisc;
pub use crate::mesh::{Mesh, Element};
pub struct TemperalDisc<'a> {
    pub residuals: ArrayViewMut<'a, f64, Ix3>,
    pub solver: &'a Solver<'a>,
    pub spatial_disc: &'a SpatialDisc<'a>,
    pub mesh: &'a Mesh<'a>,
    pub solver_param: &'a SolverParameters,
    pub time_scheme: Box<dyn TimeScheme>,
}
impl<'a> TemperalDisc<'a> {
    pub fn set_residual_to_zero(mut residuals: ArrayViewMut<f64, Ix3>) {
        for residual in residuals.iter_mut() {
            *residual = 0.0;
        }
    }
}
trait TimeScheme {
    fn step(&mut self, solutions: ArrayViewMut<f64, Ix3>, residuals: ArrayViewMut<f64, Ix3>, time_step: f64);
}