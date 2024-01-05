use ndarray::{Array, ArrayViewMut, Zip, ArrayView};
use ndarray::Ix3;
use crate::solver::Solver;
use crate::spatial_disc::SpatialDisc;

use super::TemperalDisc;
pub struct RungeKutta3rd {
    pub u1: Array<f64, Ix3>,
    pub u2: Array<f64, Ix3>,
}
impl RungeKutta3rd {
    pub fn new(nelem: usize, nbasis: usize, neq: usize) -> RungeKutta3rd {
        let u1 = Array::zeros((nelem, neq, nbasis));
        let u2 = Array::zeros((nelem, neq, nbasis));
        RungeKutta3rd {
            u1,
            u2,
        }
    }
}
impl RungeKutta3rd {
    pub fn step(&mut self, spatial_disc: &SpatialDisc, residuals: &mut Array<f64, Ix3>, solutions: &mut Array<f64, Ix3>, time_step: f64) {
        TemperalDisc::set_residual_to_zero(residuals);
        spatial_disc.compute_residuals(residuals, solutions);
        self.u1 = &*solutions + time_step * &*residuals;
        // limiter
        TemperalDisc::set_residual_to_zero(residuals);
        spatial_disc.compute_residuals(residuals, &self.u1);
        self.u2 = 0.75 * &*solutions + 0.25 * (&self.u1 + time_step * &*residuals);
        // limiter
        TemperalDisc::set_residual_to_zero(residuals);
        spatial_disc.compute_residuals(residuals, &self.u2);
        *solutions = 1.0 / 3.0 * &*solutions + 2.0 / 3.0 * (&self.u2 + time_step * &*residuals);
    }
}