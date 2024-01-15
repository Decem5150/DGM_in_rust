use ndarray::Array;
use ndarray::Ix3;
use crate::basis_function::DubinerBasis;
use crate::basis_function::GaussPoints;
use crate::io::initial_solution_parser::InitialSolution;
use crate::mesh::Mesh;
use crate::spatial_disc::SpatialDisc;
use crate::temporal_disc::TemperalDisc;
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
    pub residuals: Array<f64, Ix3>,
    pub solutions: Array<f64, Ix3>,
    pub spatial_disc: SpatialDisc<'a>,
    pub temporal_disc: TemperalDisc<'a>,
    pub mesh: &'a Mesh,
    pub basis: &'a DubinerBasis,
    pub gauss_points: &'a GaussPoints,
    pub flow_param: &'a FlowParameters,
    pub mesh_param: &'a MeshParameters,
    pub solver_param: &'a SolverParameters,
}
impl<'a> Solver<'a> {
    pub fn set_initial_solution(&mut self) {
        let initial_solution = InitialSolution::parse();
        let initial_q = {
            let density = initial_solution.density;
            let pressure = initial_solution.pressure;
            let sound_speed = (self.flow_param.hcr * pressure / density).sqrt();
            let x_velocity = initial_solution.x_mach_number * sound_speed;
            let y_velocity = initial_solution.y_mach_number * sound_speed;
            let total_energy = pressure / (self.flow_param.hcr - 1.0) + 0.5 * density * (x_velocity.powi(2) + y_velocity.powi(2));
            [
                density,
                density * x_velocity,
                density * y_velocity,
                total_energy
            ]
        };
        for ielement in 0..self.mesh_param.number_of_elements {
            for ieq in 0..self.solver_param.number_of_equations {
                for ibasis in 0..self.solver_param.number_of_basis_functions {
                    for igp in 0..self.solver_param.number_of_cell_gp {
                        self.solutions[[ielement, ieq, ibasis]] += 0.5 * self.mesh.elements[ielement].jacob_det * self.gauss_points.cell_weights[igp] * initial_q[ieq] * self.basis.phis_cell_gps[[igp, ibasis]];
                    }
                }
            }
        }
        for ielement in 0..self.mesh_param.number_of_elements {
            for ieq in 0..self.solver_param.number_of_equations {
                for ibasis in 0..self.solver_param.number_of_basis_functions {
                    self.solutions[[ielement, ieq, ibasis]] /= self.basis.mass_mat_diag[ibasis] * 0.5 * self.mesh.elements[ielement].jacob_det;
                }
            }
        }
    }
    pub fn solve(&mut self) {
        self.temporal_disc.time_march(&self.spatial_disc, &mut self.residuals, &mut self.solutions);
    }
}