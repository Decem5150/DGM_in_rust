use crate::mesh;
use crate::solver::ConsVar;
trait TimeScheme {
    fn compute(&self, mesh: &mesh::Mesh, time_step: f64, time: f64, u: &Vec<f64>, f: &dyn Fn(&Vec<f64>, &Vec<f64>, f64) -> Vec<f64>) -> Vec<f64>;
}
struct RungeKutta4th;
impl TimeScheme for RungeKutta4th {
    fn compute() {
        compute_residuals();
    }
}
pub fn runge_kutta_4th<'a>(mesh: &'a mesh::Mesh, delta_t: f64, time: f64, u: &'a Vec<f64>, f: &'a dyn Fn(&'a Vec<f64>, &'a Vec<f64>, f64) -> Vec<f64>) -> Vec<f64> {
    
}