use ndarray::Array;
use ndarray::{Ix1, Ix2, Ix3};
use ndarray::array;
use ndarray::{ArrayView, ArrayViewMut};
use crate::mesh::{Mesh, Element};
use crate::solver::{ConsVar, SolCoeff, Solver};
use crate::spatial_disc::SpatialDisc;
use super::{TemperalDisc, TimeScheme};
struct RungeKutta3rd<'a> {
    pub spatial_disc: &'a SpatialDisc<'a>,
    pub u1:Array<f64, Ix3>,
    pub u2:Array<f64, Ix3>,
}
impl<'a> TimeScheme for RungeKutta3rd<'a> {
    fn step(&mut self, mut u: ArrayViewMut<f64, Ix3>, mut residuals: ArrayViewMut<f64, Ix3>, time_step: f64) {
        TemperalDisc::set_residual_to_zero(residuals);
        self.spatial_disc.compute_residuals(residuals, u);
        residuals *= time_step;
        self.u1 += u;

    }
}