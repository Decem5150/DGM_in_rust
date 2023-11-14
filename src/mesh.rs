use ndarray::Array;
use ndarray::{Ix1, Ix2, Ix3};
use ndarray::array;
use ndarray::{ArrayView, ArrayViewMut};
use crate::solver::{ConsVar, SolCoeff};
use crate::spatial_disc::gauss_point::GaussPoints;
use crate::spatial_disc::basis_function::DubinerBasis;
use crate::spatial_disc::boundary::BoundaryCondition;
pub struct Mesh<'a> { 
    pub elements: Array<Element<'a>, Ix1>,
    pub vertices: Array<Vertex, Ix1>,
    pub edges : Array<Edge<'a>, Ix1>,
    pub patches: Array<Patch<'a>, Ix1>,
}
pub struct Element<'a> {
    pub vertices: Array<&'a Vertex, Ix1>,
    pub edges: Array<&'a Edge<'a>, Ix1>,
    pub neighbours: Array<&'a Element<'a>, Ix1>,
    pub jacob_det: f64,
    pub mass_mat_diag: Array<f64, Ix1>,
    pub dphis_dx: Array<f64, Ix2>,
    pub dphis_dy: Array<f64, Ix2>,
    pub solution_index: usize
}
impl<'a> Element<'a> {
    pub fn compute_mass_mat(&mut self, quad: &GaussPoints, basis: &DubinerBasis) {
        self.mass_mat_diag.iter_mut().for_each(|mass| *mass = 0.0);
        for i in 0..basis.dof {
            for j in 0..quad.cell_points.len() {
                self.mass_mat_diag[i] += quad.cell_weights[j] * basis.phis_cell_gps[[j, i]] * basis.phis_cell_gps[[j, j]] * self.jacob_det;
            }
        }
    }
    pub fn compute_jacob_det(&mut self) {
        let x1 = self.vertices[0].x;
        let x2 = self.vertices[1].x;
        let x3 = self.vertices[2].x;
        let y1 = self.vertices[0].y;
        let y2 = self.vertices[1].y;
        let y3 = self.vertices[2].y;
        self.jacob_det = (x2 - x1) * (y3 - y1) - (x3 - x1) * (y2 - y1);
    }
    pub fn compute_derivatives(&mut self, basis: &DubinerBasis) {
        let x1 = self.vertices[0].x;
        let x2 = self.vertices[1].x;
        let x3 = self.vertices[2].x;
        let y1 = self.vertices[0].y;
        let y2 = self.vertices[1].y;
        let y3 = self.vertices[2].y;

        let dx_deta = x2 - x1;
        let dx_dxi = x3 - x1;
        let dy_deta = y2 - y1;
        let dy_dxi = y3 - y1;

        let deta_dx = dy_dxi / self.jacob_det;
        let deta_dy = -dx_dxi / self.jacob_det;
        let dxi_dx = -dy_deta / self.jacob_det;
        let dxi_dy = dx_deta / self.jacob_det;

        for i in 0..self.dphis_dx.shape()[0] {
            for j in 0..self.dphis_dx.shape()[1] {
                self.dphis_dx[[i, j]] = basis.dphis_dxi[[i, j]] * dxi_dx + basis.dphis_deta[[i, j]] * deta_dx;
                self.dphis_dy[[i, j]] = basis.dphis_dxi[[i, j]] * dxi_dy + basis.dphis_deta[[i, j]] * deta_dy;
            }
        }
    }
}
#[derive(Default)]
pub struct Vertex {
    pub x: f64,
    pub y: f64,
}
pub struct Edge<'a> {
    pub vertices: [&'a Vertex; 2],
    pub elements: [&'a Element<'a>; 2],
    pub jacob_det: f64,
    pub normal: [f64; 2],
    pub ind_in_left_elem: usize,
    pub ind_in_right_elem: usize,
}
impl<'a> Edge<'a> {
    pub fn compute_jacob_det(&mut self) {
        let x1 = self.vertices[0].x;
        let x2 = self.vertices[1].x;
        let y1 = self.vertices[0].y;
        let y2 = self.vertices[1].y;
        self.jacob_det = ((x2 - x1).powi(2) + (y2 - y1).powi(2)).sqrt();
    }
    pub fn compute_normal(&mut self) {
        let x1 = self.vertices[0].x;
        let x2 = self.vertices[1].x;
        let y1 = self.vertices[0].y;
        let y2 = self.vertices[1].y;
        self.normal[0] = (y2 - y1) / self.jacob_det;
        self.normal[1] = -(x2 - x1) / self.jacob_det;
    } 
}
pub struct BoundaryEdge<'a> {
    pub vertices: [&'a Vertex; 2],
    pub internal_element: &'a Element<'a>,
    pub jacob_det: f64,
    pub normal: [f64; 2],
    pub ind_in_internal_elem: usize,
}
pub struct Patch<'a> {
    pub boundary_edges: Array<&'a BoundaryEdge<'a>, Ix1>,
    pub bc: Box<dyn BoundaryCondition>,
}
impl<'a> Patch<'a> {
    pub fn apply_bc(&mut self) {
        self.bc.apply();
    }
}