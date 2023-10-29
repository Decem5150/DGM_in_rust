use crate::solver::{self, ConsVar};
use crate::spatialdisc::gauss_point::Quadrature2D;
use crate::spatialdisc::basis_function::DubinerBasis;

#[derive(Default)]
pub struct Mesh<'a> { 
    pub elements: Vec<Element<'a>>,
    pub vertices: Vec<Vertex>,
    pub edges : Vec<Edge<'a>>,
    pub patches: Vec<Patch<'a>>,
}
#[derive(Default)]
pub struct Element<'a> {
    pub vertices: Vec<&'a Vertex>,
    pub edges: Vec<&'a Edge<'a>>,
    pub neighbours: Vec<&'a Element<'a>>,
    pub solution: Vec<ConsVar>,
    pub jacob_det: f64,
    pub mass_mat_diag: Vec<f64>,
    pub dphi_dx: Vec<Vec<f64>>,
    pub dphi_dy: Vec<Vec<f64>>,
}
impl<'a> Element<'a> {
    pub fn compute_jacob_det(&mut self) {
        let x1 = self.vertices[0].x;
        let x2 = self.vertices[1].x;
        let x3 = self.vertices[2].x;
        let y1 = self.vertices[0].y;
        let y2 = self.vertices[1].y;
        let y3 = self.vertices[2].y;
        self.jacob_det = (x2 - x1) * (y3 - y1) - (x3 - x1) * (y2 - y1);
    }
    pub fn compute_mass_mat(&mut self, quad: &Quadrature2D, basis: &DubinerBasis) {
        for value in self.mass_mat_diag.iter_mut() {
            *value = 0.0;
        }
        for i in 0..basis.dof {
            for j in 0..quad.points.len() {
                self.mass_mat_diag[i] += quad.weights[j] * basis.phi_gauss[j][i] * basis.phi_gauss[j][j] * self.jacob_det;
            }
        }
    }
}
#[derive(Default)]
pub struct Vertex {
    pub x: f64,
    pub y: f64,
}
#[derive(Default)]
pub struct Edge<'a> {
    pub vertices: [&'a Vertex; 2],
    pub elements: [&'a Element<'a>; 2],
    pub invis_flux: Vec<ConsVar>,
    //pub vis_flux: Vec<ConsVar>,
    pub jacob_det: f64,
    pub normal: [f64; 2],
    pub index_in_left: usize,
    pub index_in_right: usize,
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
#[derive(Default)]
pub struct Patch<'a> {
    pub edges: Vec<&'a Edge<'a>>,
    pub bc: Box<dyn BoundaryCondition>,
}
impl<'a> Patch<'a> {
    pub fn apply_bc(&mut self) {
        bc.apply();
    }
}