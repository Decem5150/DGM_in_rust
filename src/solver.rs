use super::spatial_disc;
use super::mesh;
#[derive(Default)]
#[derive(Clone)]
pub struct ConsVar {
    pub density: f64,
    pub x_momentum: f64,
    pub y_momentum: f64,
    pub energy: f64,
}
pub struct SolverParameters {
    pub cfl: f64,
    pub final_time: f64,
    pub number_of_time_steps: i64,
    pub number_of_gauss_points: i64,
    pub number_of_elements: i64,
    pub number_of_basis_functions: i64,
    pub number_of_edges: i64,
    pub number_of_vertices: i64,
    pub number_of_patches: i64,
}
struct FlowParameters {
    pub hcr: f64,
    pub gas_constant: f64,
    pub viscosity: f64,
    pub prandtl_number: f64,
}
pub struct Solver<'a> {
    pub mesh: mesh::Mesh<'a>,
    pub residuals: Vec<Vec<ConsVar>>,
    pub spatial_disc: spatial_disc::SpatialDisc<'a>,
    pub temperal_disc: Option<fn()>,
    pub solver_param: SolverParameters, 
    pub flow_param: FlowParameters,
}
impl<'a> Solver<'a> {
    pub fn compute_residuals(&mut self, u: &Vec<f64>, time: f64) {
        self.mesh.
        self.spatial_disc.compute_fluxes(self.residuals: &mut Vec<Vec<ConsVar>>, mesh: &mesh::Mesh<'a>);
    }
    pub fn time_step(&mut self, u: &Vec<f64>, time: f64) {
        
    }
    pub fn solve(&mut self, u: &Vec<f64>, time: f64) {
        for step in 0..self.solver_param.number_of_time_steps {
            self.time_step(u, time);
        }
        self.compute_residuals(u, time);
    }
}