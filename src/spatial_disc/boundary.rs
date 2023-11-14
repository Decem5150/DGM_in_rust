use ndarray::{Array, Ix2};
use ndarray::{ArrayView, ArrayViewMut};
use crate::mesh::Element;
pub trait BoundaryCondition<'a> {
    fn apply(&self, left_values_gps: ArrayView<f64, Ix2>, right_values_gps: ArrayViewMut<f64, Ix2>, normal: [f64; 2], hcr: f64);
}
struct NoSlipWall;
impl<'a> BoundaryCondition<'a> for NoSlipWall {
    fn apply(&self, left_values_gps: ArrayView<f64, Ix2>, mut right_values_gps: ArrayViewMut<f64, Ix2>, normal: [f64; 2], hcr: f64) {
        for igp in 0..left_values_gps.shape()[0] {
            let pl = (hcr - 1.0) * (left_values_gps[[igp, 3]] - 0.5 * (left_values_gps[[igp, 1]] * left_values_gps[[igp, 1]] + left_values_gps[[igp, 2]] * left_values_gps[[igp, 2]]) / left_values_gps[[igp, 0]]);
            right_values_gps[[igp, 0]] = left_values_gps[[igp, 0]];
            right_values_gps[[igp, 1]] = -left_values_gps[[igp, 1]];
            right_values_gps[[igp, 2]] = -left_values_gps[[igp, 2]];
            right_values_gps[[igp, 3]] = pl / (hcr - 1.0) + 0.5 * (right_values_gps[[igp, 1]] * right_values_gps[[igp, 1]] + right_values_gps[[igp, 2]] * right_values_gps[[igp, 2]]) / right_values_gps[[igp, 0]];
        }
    }
}
struct FreeSlipWall;
impl<'a> BoundaryCondition<'a> for FreeSlipWall {
    fn apply(&self, left_values_gps: ArrayView<f64, Ix2>, mut right_values_gps: ArrayViewMut<f64, Ix2>, normal: [f64; 2], hcr: f64) {
        let nx = normal[0];
        let ny = normal[1];
        for igp in 0..left_values_gps.shape()[0] {
            let pl = (hcr - 1.0) * (left_values_gps[[igp, 3]] - 0.5 * (left_values_gps[[igp, 1]] * left_values_gps[[igp, 1]] + left_values_gps[[igp, 2]] * left_values_gps[[igp, 2]]) / left_values_gps[[igp, 0]]);
            let vnl = left_values_gps[[igp, 1]] * nx + left_values_gps[[igp, 2]] * ny;
            right_values_gps[[igp, 0]] = left_values_gps[[igp, 0]];
            right_values_gps[[igp, 1]] = left_values_gps[[igp, 1]] - 2.0 * nx * vnl;
            right_values_gps[[igp, 2]] = left_values_gps[[igp, 2]] - 2.0 * ny * vnl;
            right_values_gps[[igp, 3]] = pl / (hcr - 1.0) + 0.5 * (right_values_gps[[igp, 1]] * right_values_gps[[igp, 1]] + right_values_gps[[igp, 2]] * right_values_gps[[igp, 2]]) / right_values_gps[[igp, 0]];
        }
    }
}
struct FarField;
impl<'a> BoundaryCondition<'a> for FarField {
    fn apply(&self, left_values_gps: ArrayView<f64, Ix2>, mut right_values_gps: ArrayViewMut<f64, Ix2>, normal: [f64; 2], hcr: f64) {
        let nx = normal[0];
        let ny = normal[1];
        for igp in 0..left_values_gps.shape()[0] {
            
        }
    }
}