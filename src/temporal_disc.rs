pub mod rk3;
use ndarray::Array;
use ndarray::{Ix1, Ix2, Ix3};
use ndarray::{ArrayView, ArrayViewMut};
pub use crate::solver::{SolverParameters, Solver, SolCoeff};
pub use crate::spatial_disc::SpatialDisc;
pub use crate::mesh::{Mesh, Element};
pub struct TemperalDisc<'a> {
    pub residuals: Array<f64, Ix3>,
    pub solver: &'a Solver<'a>,
    pub spatial_disc: &'a SpatialDisc<'a>,
    pub mesh: &'a Mesh<'a>,
    pub solver_param: &'a SolverParameters,
    pub current_time: f64,
    pub time_step: f64,
    pub time_scheme: Box<dyn TimeScheme>,
}
impl<'a> TemperalDisc<'a> {
    pub fn new (solver: &'a Solver<'a>, spatial_disc: &'a SpatialDisc<'a>, mesh: &'a Mesh<'a>, solver_param: &'a SolverParameters) -> Self {
        let nelem = solver_param.number_of_elements;
        let nbasis = solver_param.number_of_basis_functions;
        let neq = solver_param.number_of_equations;
        let mut residuals = Array::zeros((nelem, neq, nbasis));
        let time_scheme: Box<dyn TimeScheme> = Box::new(rk3::RungeKutta3rd::new(spatial_disc));
        TemperalDisc {
            residuals,
            solver,
            spatial_disc,
            mesh,
            solver_param,
            current_time: 0.0,
            time_step: 0.0,
            time_scheme,
        }
    }
    pub fn time_march(&mut self, solutions: ArrayViewMut<f64, Ix3>) {
        while self.current_time < self.solver_param.final_time {
            if self.current_time + self.time_step > self.solver_param.final_time {
                self.time_step = self.solver_param.final_time - self.current_time;
            }
            self.time_scheme.step(solutions, self.residuals.view_mut(), self.time_step);
            self.current_time += self.time_step;
        }
    }
    pub fn set_residual_to_zero(mut residuals: ArrayViewMut<f64, Ix3>) {
        for residual in residuals.iter_mut() {
            *residual = 0.0;
        }
    }
}
trait TimeScheme {
    fn step(&mut self, solutions: ArrayViewMut<f64, Ix3>, residuals: ArrayViewMut<f64, Ix3>, time_step: f64);
}