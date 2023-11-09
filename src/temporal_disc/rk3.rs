use crate::mesh::{Mesh, Element};
use crate::solver::{ConsVar, SolCoeff, Solver};
use crate::spatial_disc::SpatialDisc;
use super::{TemperalDisc, TimeScheme};
struct RungeKutta3rd<'a> {
    pub spatial_disc: &'a SpatialDisc<'a>,
}
impl<'a> TimeScheme for RungeKutta3rd<'a> {
    fn step(&self, residuals: &mut Vec<SolCoeff>, time_step: f64) {
        Solver::set_residual_to_zero(residuals);
        self.spatial_disc.compute_residuals();
        let u1 = elements.iter_mut()
            .map(|element| {
                element.solution.iter_mut()
                    .zip(&element.residuals)
                    .for_each(|(sol, res)| *sol += res * time_step);
                &mut element.solution
            })
    }
}