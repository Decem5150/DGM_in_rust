pub mod rk3;
use ndarray::{Array, Ix3};
use crate::{solver::{FlowParameters, MeshParameters, SolverParameters}, mesh::Mesh, basis_function::{DubinerBasis, GaussPoints}, spatial_disc::SpatialDisc};
pub enum TimeScheme {
    RK3,
}
pub struct TemperalDisc<'a> {
    pub temporary_solutions: Vec<Array<f64, Ix3>>,
    pub current_time: f64,
    pub current_step: usize,
    pub time_scheme: TimeScheme,
    pub mesh: &'a Mesh,
    pub basis: &'a DubinerBasis,
    pub gauss_points: &'a GaussPoints,
    pub flow_param: &'a FlowParameters,
    pub mesh_param: &'a MeshParameters,
    pub solver_param: &'a SolverParameters,
}
impl<'a> TemperalDisc<'a> {
    pub fn time_march(&mut self, spatial_disc: &SpatialDisc, residuals: &mut Array<f64, Ix3>, solutions: &mut Array<f64, Ix3>) {
        let mut time_step = 1.0e10;
        while self.current_time < self.solver_param.final_time {
            if self.current_time + time_step > self.solver_param.final_time {
                time_step = self.solver_param.final_time - self.current_time;
            } else {
                time_step = self.compute_time_step(solutions);
            }
            match self.time_scheme {
                TimeScheme::RK3 => self.rk3(spatial_disc, residuals, solutions, time_step),
            }
            self.current_time += time_step;
        }
    }
    fn compute_time_step(&self, solutions: &mut Array<f64, Ix3>) -> f64 {
        let nbasis = self.solver_param.number_of_basis_functions;
        let nelem = self.mesh_param.number_of_elements;
        let hcr = self.flow_param.hcr;
        let mut time_step = 1.0e10;
        for ielem in 0..nelem {
            let mut rho = 0.0;
            let mut rho_u = 0.0;
            let mut rho_v = 0.0;
            let mut rho_e = 0.0;
            for ibasis in 0..nbasis {
                rho += solutions[[ielem, 0, ibasis]];
                rho_u += solutions[[ielem, 1, ibasis]];
                rho_v += solutions[[ielem, 2, ibasis]];
                rho_e += solutions[[ielem, 3, ibasis]];
            }
            let rho = rho;
            let u = rho_u / rho;
            let v = rho_v / rho;
            let p = (hcr - 1.0) * (rho_e - 0.5 * rho * (u * u + v * v));
            let c = (hcr * p / rho).sqrt();
            let speed = u + c;
            let dx = self.mesh.elements[ielem].circumradius;
            let dt = self.solver_param.cfl * dx / speed;
            if dt < time_step {
                time_step = dt;
            }
        }
        time_step
    }
    pub fn set_residual_to_zero(residuals: &mut Array<f64, Ix3>) {
        for residual in residuals.iter_mut() {
            *residual = 0.0;
        }
    }
}