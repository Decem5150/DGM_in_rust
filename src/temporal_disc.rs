pub mod rk3;
use ndarray::{Array, ArrayViewMut, Ix3, ArrayView};
use crate::solver::{SolverParameters, MeshParameters, FlowParameters};
use crate::spatial_disc::SpatialDisc;
use crate::mesh::Mesh;
use self::rk3::RungeKutta3rd;
pub enum TimeScheme {
    RK3(RungeKutta3rd),
}
pub struct TemperalDisc<'a> {
    pub residuals: Array<f64, Ix3>,
    pub mesh: &'a Mesh<'a>,
    pub spatial_disc: &'a SpatialDisc<'a>,
    pub solver_param: &'a SolverParameters,
    pub mesh_param: &'a MeshParameters,
    pub flow_param: &'a FlowParameters,
    pub current_time: f64,
    pub time_step: f64,
    pub time_scheme: TimeScheme,
}
impl<'a> TemperalDisc<'a> {
    pub fn new (mesh: &Mesh, spatial_disc: &SpatialDisc, solver_param: &SolverParameters, mesh_param: &MeshParameters, flow_param: &FlowParameters) -> Self {
        let nelem = mesh_param.number_of_elements;
        let nbasis = solver_param.number_of_basis_functions;
        let neq = solver_param.number_of_equations;
        let mut residuals = Array::zeros((nelem, neq, nbasis));
        let time_scheme = TimeScheme::RK3(RungeKutta3rd::new(spatial_disc));
        TemperalDisc {
            residuals,
            mesh,
            spatial_disc,
            solver_param,
            mesh_param,
            flow_param,
            current_time: 0.0,
            time_step: 0.0,
            time_scheme,
        }
    }
    pub fn time_march(&mut self, solutions: ArrayViewMut<f64, Ix3>) {
        while self.current_time < self.solver_param.final_time {
            if self.current_time + self.time_step > self.solver_param.final_time {
                self.time_step = self.solver_param.final_time - self.current_time;
            } else {
                self.compute_time_step(solutions.view());
            }
            match self.time_scheme {
                TimeScheme::RK3(ref mut rk3) => rk3.step(self.spatial_disc, solutions, self.residuals.view_mut(), self.time_step),
            }
            self.current_time += self.time_step;
        }
    }
    pub fn compute_time_step(&mut self, solutions: ArrayView<f64, Ix3>) {
        let mut time_step = 1.0e10;
        for ielem in 0..self.mesh_param.number_of_elements {
            let mut u_max = 0.0;
            for ibasis in 0..self.solver_param.number_of_basis_functions {
                let u = solutions[[ielem, 0, ibasis]];
                let v = solutions[[ielem, 1, ibasis]];
                let p = solutions[[ielem, 3, ibasis]];
                let hcr = self.flow_param.hcr;
                let a = (hcr * p).sqrt();
                let u = (u * u + v * v).sqrt();
                let c = (hcr * p).sqrt();
                let c = (c / solutions[[ielem, 0, ibasis]]).sqrt();
                let u = u + c;
                if u > u_max {
                    u_max = u;
                }
            }
            let dx = self.mesh.elements[ielem].circumradius;
            let dt = self.solver_param.cfl * dx / u_max;
            if dt < time_step {
                time_step = dt;
            }
        }
        self.time_step = time_step;
    }
    pub fn set_residual_to_zero(mut residuals: ArrayViewMut<f64, Ix3>) {
        for residual in residuals.iter_mut() {
            *residual = 0.0;
        }
    }
}