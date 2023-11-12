use ndarray::Array;
use ndarray::Ix3;
use ndarray::array;
use ndarray::{ArrayView, ArrayViewMut};
use super::spatial_disc;
use super::mesh;
use crate::spatial_disc::gauss_point::GaussPoints;
use crate::spatial_disc::basis_function::DubinerBasis;
use crate::temporal_disc;
#[derive(Default)]
#[derive(Clone)]
pub struct ConsVar {
    pub density: f64,
    pub x_momentum: f64,
    pub y_momentum: f64,
    pub energy: f64,
}
impl ConsVar {
    pub fn iter(&self) -> ConsVarIterator {
        ConsVarIterator {
            cons_var: self,
            index: 0,
        }
    }
}
impl<'a> IntoIterator for &'a ConsVar {
    type Item = &'a f64;
    type IntoIter = ConsVarIterator<'a>;
    fn into_iter(self) -> Self::IntoIter {
        ConsVarIterator {
            cons_var: &self,
            index: 0,
        }
    }
}
pub struct ConsVarIterator<'a> {
    pub cons_var: &'a ConsVar,
    pub index: usize,
}
impl<'a> Iterator for ConsVarIterator<'a> {
    type Item = &'a f64;
    fn next(&mut self) -> Option<Self::Item> {
        self.index += 1;
        match self.index {
            1 => Some(&self.cons_var.density),
            2 => Some(&self.cons_var.x_momentum),
            3 => Some(&self.cons_var.y_momentum),
            4 => Some(&self.cons_var.energy),
            _ => None,
        }
    }
}
#[derive(Default)]
#[derive(Clone)]
pub struct SolCoeff {
    pub density: Vec<f64>,
    pub x_momentum: Vec<f64>,
    pub y_momentum: Vec<f64>,
    pub energy: Vec<f64>,
}
impl SolCoeff {
    pub fn iter(&self) -> SolCoeffIterator {
        SolCoeffIterator {
            sol_coeff: self,
            index: 0,
        }
    }
    pub fn iter_mut(&mut self) -> SolCoeffIteratorMut {
        SolCoeffIteratorMut {
            sol_coeff: self,
            index: 0,
        }
    }  
}
impl<'a> IntoIterator for &'a SolCoeff {
    type Item = &'a Vec<f64>;
    type IntoIter = SolCoeffIterator<'a>;
    fn into_iter(self) -> Self::IntoIter {
        SolCoeffIterator {
            sol_coeff: self,
            index: 0,
        }
    }
}
impl<'a> IntoIterator for &'a mut SolCoeff {
    type Item = &'a mut Vec<f64>;
    type IntoIter = SolCoeffIteratorMut<'a>;
    fn into_iter(self) -> Self::IntoIter {
        SolCoeffIteratorMut {
            sol_coeff: self,
            index: 0,
        }
    }
}
pub struct SolCoeffIterator<'a> {
    pub sol_coeff: &'a SolCoeff,
    pub index: usize,
}
impl<'a> Iterator for SolCoeffIterator<'a> {
    type Item = &'a Vec<f64>;
    fn next(&mut self) -> Option<Self::Item> {
        self.index += 1;
        match self.index {
            1 => Some(&self.sol_coeff.density),
            2 => Some(&self.sol_coeff.x_momentum),
            3 => Some(&self.sol_coeff.y_momentum),
            4 => Some(&self.sol_coeff.energy),
            _ => None,
        }
    }
}
pub struct SolCoeffIteratorMut<'a> {
    pub sol_coeff: &'a mut SolCoeff,
    pub index: usize,
}
impl<'a> Iterator for SolCoeffIteratorMut<'a> {
    type Item = &'a mut Vec<f64>;
    fn next(&mut self) -> Option<Self::Item> {
        self.index += 1;
        match self.index {
            1 => Some(&mut self.sol_coeff.density),
            2 => Some(&mut self.sol_coeff.x_momentum),
            3 => Some(&mut self.sol_coeff.y_momentum),
            4 => Some(&mut self.sol_coeff.energy),
            _ => None,
        }
    }
}
pub struct SolverParameters {
    pub cfl: f64,
    pub final_time: f64,
    pub number_of_time_steps: usize,
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
    pub fn time_step(&mut self, u: &Vec<f64>, time: f64) {
        
    }
    pub fn solve(&mut self, u: &Vec<f64>, time: f64) {
        /*
        let mut residuals:Array<f64, Ix3> = Array::zeros((self.solver_param.number_of_elements, self.solver_param.number_of_equations, self.solver_param.number_of_basis_functions));
        let mut u1:Array<f64, Ix3> = Array::zeros((self.solver_param.number_of_elements, self.solver_param.number_of_equations, self.solver_param.number_of_basis_functions));
        let mut u2:Array<f64, Ix3> = Array::zeros((self.solver_param.number_of_elements, self.solver_param.number_of_equations, self.solver_param.number_of_basis_functions));
        let mut u3:Array<f64, Ix3> = Array::zeros((self.solver_param.number_of_elements, self.solver_param.number_of_equations, self.solver_param.number_of_basis_functions));
        */
        let u_vec = vec!(self.solutions.view_mut(), u1.view_mut(), u2.view_mut(), u3.view_mut());
        for step in 0..self.solver_param.number_of_time_steps {
            self.time_step(u, time);
        }
    }
}