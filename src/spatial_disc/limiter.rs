use ndarray::Ix3;
use ndarray::Array;
pub mod weno;

use super::SpatialDisc;
impl<'a> SpatialDisc<'a> {
    pub fn apply_limiter(&self, solutions: &mut Array<f64, Ix3>) {
        self.weno(solutions);
    }
}