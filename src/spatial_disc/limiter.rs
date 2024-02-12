use ndarray::{Ix1, Ix3};
use ndarray::Array;
pub mod weno;

use super::SpatialDisc;
pub enum LimiterType {
    WENO,
    None,
}
pub enum DetectorType {
    TVB,
    KXRCF,
    None,
}
impl<'a> SpatialDisc<'a> {
    pub fn apply_limiter(&self, solutions: &mut Array<f64, Ix3>) {
        self.weno(solutions);
    }
}