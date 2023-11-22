use ndarray::Array;
use ndarray::Ix3;
use super::spatial_disc;
use super::mesh;
use crate::temporal_disc;
pub struct SolverParameters {
    pub cfl: f64,
    pub final_time: f64,
    pub number_of_cell_gp: usize,
    pub number_of_edge_gp: usize,
    pub number_of_equations: usize,
    pub number_of_basis_functions: usize,
    pub number_of_elements: usize,
    pub number_of_edges: usize,
    pub number_of_vertices: usize,
    pub number_of_patches: usize,
}
pub struct FlowParameters {
    pub hcr: f64,
    pub gas_constant: f64,
    pub viscosity: f64,
    pub prandtl_number: f64,
}
pub struct Solver<'a> {
    pub mesh: mesh::Mesh<'a>,
    pub solutions: Array<f64, Ix3>,
    pub spatial_disc: spatial_disc::SpatialDisc<'a>,
    pub temperal_disc: temporal_disc::TemperalDisc<'a>,
    pub solver_param: SolverParameters, 
    pub flow_param: FlowParameters,
}
impl<'a> Solver<'a> {
    pub fn solve(&mut self) {
        self.temperal_disc.time_march(self.solutions.view_mut());
    }
}