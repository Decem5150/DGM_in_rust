use std::cell::RefCell;
use std::rc::Weak;
use ndarray::{Array, ArrayView};
use ndarray::{Ix1, Ix2};
use crate::spatial_disc::SpatialDisc;
use super::spatial_disc::gauss_point::GaussPoints;
use super::spatial_disc::basis_function::DubinerBasis;
pub struct Mesh<'a> { 
    pub elements: Array<Element, Ix1>,
    pub edges : Array<Edge, Ix1>,
    pub internal_elements: Array<usize, Ix1>,
    pub boundary_elements: Array<usize, Ix1>,
    pub boundary_edges: Array<BoundaryEdge, Ix1>,
    pub vertices: Array<Vertex, Ix1>,
    pub patches: Array<Patch, Ix1>,
    pub spatial_disc: RefCell<Weak<SpatialDisc<'a>>>,
}
impl<'a> Mesh<'a> {
    pub fn compute_mass_mat(&mut self) {
        let spatial_disc = self.spatial_disc.borrow().upgrade().unwrap();
        for element in self.elements.iter_mut() {
            element.mass_mat_diag.iter_mut().for_each (|mass| *mass = 0.0);
            for ibasis in 0..spatial_disc.basis.dof {
                for igp in 0..spatial_disc.gauss_point.cell_gp_number {
                    element.mass_mat_diag[ibasis] += self.spatial_disc.borrow().gauss_point.cell_weights[igp] * self.spatial_disc.borrow().basis.phis_cell_gps[[igp, ibasis]] * self.spatial_disc.borrow().basis.phis_cell_gps[[igp, ibasis]] * element.jacob_det;
                }
            }
        }
    }
    pub fn compute_jacob_det(&mut self) {
        for element in self.elements.iter_mut() {
            let x1 = self.vertices[element.ivertices[0]].x;
            let x2 = self.vertices[element.ivertices[1]].x;
            let x3 = self.vertices[element.ivertices[2]].x;
            let y1 = self.vertices[element.ivertices[0]].y;
            let y2 = self.vertices[element.ivertices[1]].y;
            let y3 = self.vertices[element.ivertices[2]].y;
            element.jacob_det = (x2 - x1) * (y3 - y1) - (x3 - x1) * (y2 - y1);
        }
    }
    pub fn compute_derivatives(&mut self) {
        let spatial_disc = self.spatial_disc.borrow().upgrade().unwrap();
        for element in self.elements.iter_mut() {
            let x1 = self.vertices[element.ivertices[0]].x;
            let x2 = self.vertices[element.ivertices[1]].x;
            let x3 = self.vertices[element.ivertices[2]].x;
            let y1 = self.vertices[element.ivertices[0]].y;
            let y2 = self.vertices[element.ivertices[1]].y;
            let y3 = self.vertices[element.ivertices[2]].y;
            let dx_deta = x2 - x1;
            let dx_dxi = x3 - x1;
            let dy_deta = y2 - y1;
            let dy_dxi = y3 - y1;
            let deta_dx = dy_dxi / element.jacob_det;
            let deta_dy = -dx_dxi / element.jacob_det;
            let dxi_dx = -dy_deta / element.jacob_det;
            let dxi_dy = dx_deta / element.jacob_det;
            for ibasis in 0..spatial_disc.basis.dof {
                for igp in 0..spatial_disc.gauss_point.cell_gp_number {
                    element.dphis_dx[[ibasis, igp]] = spatial_disc.basis.dphis_dxi[[igp, ibasis]] * dxi_dx + spatial_disc.basis.dphis_deta[[igp, ibasis]] * deta_dx;
                    element.dphis_dy[[ibasis, igp]] = spatial_disc.basis.dphis_dxi[[igp, ibasis]] * dxi_dy + spatial_disc.basis.dphis_deta[[igp, ibasis]] * deta_dy;
                }
            }
        }
    }
}
pub struct Vertex {
    pub x: f64,
    pub y: f64,
}
pub enum NormalDirection {
    Inward,
    Outward,
}
pub struct Element {
    pub ivertices: Array<usize, Ix1>,
    pub iedges: Array<usize, Ix1>,
    pub ineighbours: Array<usize, Ix1>,
    pub mass_mat_diag: Array<f64, Ix1>,
    pub dphis_dx: Array<f64, Ix2>,
    pub dphis_dy: Array<f64, Ix2>,
    pub normal_directions: Array<NormalDirection, Ix1>,
    pub jacob_det: f64,
    pub circumradius: f64,
}
impl Element {

}
pub enum EdgeType<'a> {
    Internal(&'a Edge),
    Boundary(&'a BoundaryEdge),
}
pub struct Edge {
    pub ielements: [usize; 2],
    pub ivertices: [usize; 2],
    pub jacob_det: f64,
    pub normal: [f64; 2],
    pub in_cell_index: [usize; 2],
    //pub hcr: f64,
}
impl Edge {
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
pub struct BoundaryEdge {
    pub ivertices: [usize; 2],
    pub ielement: usize,
    pub jacob_det: f64,
    pub normal: [f64; 2],
    pub in_cell_index: usize,
    //pub hcr: f64,
}
pub enum BoundaryType {
    Wall,
    FarField,
}
pub struct BoundaryQuantity {
    pub rho: f64,
    pub u: f64,
    pub v: f64,
    pub p: f64,
}
pub struct Patch {
    pub iedges: Array<usize, Ix1>,
    pub boundary_type: BoundaryType,
    pub boundary_quantity: Option<BoundaryQuantity>,
}