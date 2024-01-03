use ndarray::{Array, ArrayViewMut, Zip};
use ndarray::Ix3;
use crate::spatial_disc::SpatialDisc;
use super::TemperalDisc;
pub struct RungeKutta3rd {
    pub u1: Array<f64, Ix3>,
    pub u2: Array<f64, Ix3>,
}
impl RungeKutta3rd {
    pub fn new (spatial_disc: &SpatialDisc) -> RungeKutta3rd {
        let nelem = spatial_disc.mesh_param.number_of_elements;
        let nbasis = spatial_disc.solver_param.number_of_basis_functions;
        let neq = spatial_disc.solver_param.number_of_equations;
        let u1 = Array::zeros((nelem, neq, nbasis));
        let u2 = Array::zeros((nelem, neq, nbasis));
        RungeKutta3rd {
            u1,
            u2,
        }
    }
    pub fn step(&mut self, spatial_disc: &SpatialDisc, mut u: ArrayViewMut<f64, Ix3>, mut residuals: ArrayViewMut<f64, Ix3>, time_step: f64) {
        TemperalDisc::set_residual_to_zero(residuals);
        spatial_disc.compute_residuals(residuals, u.view());
        Zip::from(u).and(&mut self.u1).and(residuals).for_each(|u0, u1, res| {
            *u1 = *u0 + time_step * (*res);
        });
        // limiter
        TemperalDisc::set_residual_to_zero(residuals);
        spatial_disc.compute_residuals(residuals, self.u1.view());
        Zip::from(u).and(&mut self.u1).and(&mut self.u2).and(residuals).for_each(|u0, u1, u2, res| {
            *u2 = 0.75 * (*u0) + 0.25 * (*u1 + time_step * (*res));
        });
        // limiter
        TemperalDisc::set_residual_to_zero(residuals);
        spatial_disc.compute_residuals(residuals, self.u2.view());
        Zip::from(u).and(&mut self.u2).and(residuals).for_each(|u0, u2, res| {
            *u0 = 1.0 / 3.0 * (*u0) + 2.0 / 3.0 * (*u2 + time_step * (*res));
        });
    }
}