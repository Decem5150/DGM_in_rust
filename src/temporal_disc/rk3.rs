use ndarray::Array;
use ndarray::Ix3;
use crate::debug_utils::check_for_nan;
use crate::spatial_disc::SpatialDisc;
use super::TemperalDisc;
impl<'a> TemperalDisc<'a> {
    pub fn rk3(&mut self, spatial_disc: &SpatialDisc, residuals: &mut Array<f64, Ix3>, solutions: &mut Array<f64, Ix3>, time_step: f64) {
        TemperalDisc::set_residual_to_zero(residuals);
        spatial_disc.compute_residuals(residuals, solutions);
        self.temporary_solutions[0] = &*solutions + time_step * &*residuals;
        //println!("finished rk3 first step");
        //println!("applying limiter");
        spatial_disc.apply_limiter(&mut self.temporary_solutions[0]);
        TemperalDisc::set_residual_to_zero(residuals);
        spatial_disc.compute_residuals(residuals, &self.temporary_solutions[0]);
        self.temporary_solutions[1] = 0.75 * &*solutions + 0.25 * (&self.temporary_solutions[0] + time_step * &*residuals);
        //println!("finished rk3 second step");
        //println!("applying limiter");
        spatial_disc.apply_limiter(&mut self.temporary_solutions[1]);
        TemperalDisc::set_residual_to_zero(residuals);
        spatial_disc.compute_residuals(residuals, &self.temporary_solutions[1]);
        *solutions = 1.0 / 3.0 * &*solutions + 2.0 / 3.0 * (&self.temporary_solutions[1] + time_step * &*residuals);
        //println!("finished rk3 third step");
        //println!("applying limiter");
        spatial_disc.apply_limiter(solutions);
    }
}