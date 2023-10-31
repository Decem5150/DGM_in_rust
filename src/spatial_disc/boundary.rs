use crate::solver::ConsVar;
use crate::mesh::Element;
pub trait BoundaryCondition<'a> {
    fn apply(&self, left_element: &'a mut Element, right_element: &'a mut Element, flux: &'a mut ConsVar, normal: [f64; 2], hcr: f64);
}
struct NoSlipWall;
impl<'a> BoundaryCondition<'a> for NoSlipWall {
    fn apply(&self, left_element: &'a mut Element, right_element: &'a mut Element, flux: &'a mut ConsVar, normal: [f64; 2], hcr: f64) {
        right_element.solution.density = left_element.solution.density;
        right_element.solution.x_momentum = - left_element.solution.x_momentum;
        right_element.solution.y_momentum = - left_element.solution.y_momentum;
        let p = (hcr - 1.0) * (left_element.solution.energy - 0.5 * (left_element.solution.x_momentum * left_element.solution.x_momentum + left_element.solution.y_momentum * left_element.solution.y_momentum) / left_element.solution.density);
        right_element.solution.energy = p / (hcr - 1.0) + 0.5 * (right_element.solution.x_momentum * right_element.solution.x_momentum + right_element.solution.y_momentum * right_element.solution.y_momentum) / right_element.solution.density;
    }
}
struct FreeSlipWall;
struct FarField;