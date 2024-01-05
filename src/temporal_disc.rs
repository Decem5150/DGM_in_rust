pub mod rk3;
use core::time;

use ndarray::{Array, ArrayViewMut, Ix3, ArrayView};
use crate::{solver::{Solver, FlowParameters, MeshParameters, SolverParameters}, mesh::Mesh, basis_function::{DubinerBasis, GaussPoints}, spatial_disc::SpatialDisc};
use self::rk3::RungeKutta3rd;
pub enum TimeScheme {
    RK3(RungeKutta3rd),
}
pub struct TemperalDisc<'a> {
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
                TimeScheme::RK3(ref mut rk3) => rk3.step(spatial_disc, residuals, solutions, time_step),
            }
            self.current_time += time_step;
        }
    }
    fn compute_time_step(&self, solutions: &mut Array<f64, Ix3>) -> f64 {
        let nbasis = self.solver_param.number_of_basis_functions;
        let ngp = self.solver_param.number_of_cell_gp;
        let neq = self.solver_param.number_of_equations;
        let nelem = self.mesh_param.number_of_elements;
        let mut time_step = 1.0e10;
        for ielem in 0..nelem {
            let mut u_max = 0.0;
            for ibasis in 0..nbasis {
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
        time_step
    }
    pub fn set_residual_to_zero(residuals: &mut Array<f64, Ix3>) {
        for residual in residuals.iter_mut() {
            *residual = 0.0;
        }
    }
}