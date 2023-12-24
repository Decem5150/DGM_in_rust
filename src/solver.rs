use std::cell::RefCell;
use std::rc::Rc;
use ndarray::Array;
use ndarray::Ix3;
use crate::{spatial_disc, temporal_disc, mesh, basis_function, gauss_point};
pub struct SolverParameters {
    pub cfl: f64,
    pub final_time: f64,
    pub number_of_cell_gp: usize,
    pub number_of_edge_gp: usize,
    pub number_of_equations: usize,
    pub number_of_basis_functions: usize,
}
pub struct MeshParameters {
    pub number_of_elements: usize,
    pub number_of_edges: usize,
    pub number_of_vertices: usize,
    pub number_of_patches: usize,
}
pub struct FlowParameters {
    pub hcr: f64,
    /*
    pub gas_constant: f64,
    pub viscosity: f64,
    pub prandtl_number: f64,
    */
}
pub struct Solver<'a> {
    pub solutions: Array<f64, Ix3>,
    pub mesh: Rc<mesh::Mesh>,
    pub spatial_disc: Rc<spatial_disc::SpatialDisc<'a>>,
    pub temperal_disc: Rc<temporal_disc::TemperalDisc<'a>>,
    pub basis: Rc<basis_function::DubinerBasis>,
    pub gauss_point: Rc<gauss_point::GaussPoints>,
    pub solver_param: Rc<SolverParameters>, 
    pub mesh_param: Rc<MeshParameters>,
    pub flow_param: Rc<FlowParameters>,
}
impl<'a> Solver<'a> {
    pub fn solve(&mut self) {
        self.temperal_disc.time_march(self.solutions.view_mut());
    }
}