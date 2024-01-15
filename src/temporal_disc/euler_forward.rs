use ndarray::Array;
use ndarray::Ix3;
use crate::debug_utils::check_for_nan;
use crate::spatial_disc::SpatialDisc;
use super::TemperalDisc;
impl<'a> TemperalDisc<'a> {
    pub fn euler_forward(&mut self, spatial_disc: &SpatialDisc, residuals: &mut Array<f64, Ix3>, solutions: &mut Array<f64, Ix3>, time_step: f64) {
        TemperalDisc::set_residual_to_zero(residuals);
        spatial_disc.compute_residuals(residuals, solutions);
        check_for_nan("residuals".to_string(), residuals);
        *solutions = &*solutions + time_step * &*residuals;
    }
}