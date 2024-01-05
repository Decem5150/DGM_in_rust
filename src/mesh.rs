use std::collections::HashMap;
use ndarray::{Array, ArrayView};
use ndarray::{Ix1, Ix2};
use crate::gauss_point::GaussPoints;
use crate::basis_function::DubinerBasis;
pub struct Mesh { 
    pub elements: Array<Element, Ix1>,
    pub edges: Array<Edge, Ix1>,
    pub indices_internal_elements: Array<usize, Ix1>,
    pub indices_boundary_elements: Array<usize, Ix1>,
    pub boundary_edges: Array<BoundaryEdge, Ix1>,
    pub vertices: Array<Vertex, Ix1>,
    pub patches: Array<Patch, Ix1>,
}
impl Mesh {
    pub fn compute_mass_mat(&mut self, basis: &DubinerBasis, gauss_points: &GaussPoints, ngp: usize, nbasis: usize) {
        for element in self.elements.iter_mut() {
            element.mass_mat_diag.iter_mut().for_each (|mass| *mass = 0.0);
            for ibasis in 0..nbasis {
                for igp in 0..ngp {
                    element.mass_mat_diag[ibasis] += gauss_points.cell_weights[igp] * basis.phis_cell_gps[[igp, ibasis]] * basis.phis_cell_gps[[igp, ibasis]] * element.jacob_det;
                }
            }
        }
    }
    pub fn compute_normal(&mut self) {
        for edge in self.edges.iter_mut() {
            let x1 = self.vertices[edge.ivertices[0]].x;
            let x2 = self.vertices[edge.ivertices[1]].x;
            let y1 = self.vertices[edge.ivertices[0]].y;
            let y2 = self.vertices[edge.ivertices[1]].y;
            edge.normal[0] = (y2 - y1) / edge.jacob_det;
            edge.normal[1] = -(x2 - x1) / edge.jacob_det;
        }
    }
    pub fn compute_jacob_det(&mut self) {
        for edge in self.edges.iter_mut() {
            let x1 = self.vertices[edge.ivertices[0]].x;
            let x2 = self.vertices[edge.ivertices[1]].x;
            let y1 = self.vertices[edge.ivertices[0]].y;
            let y2 = self.vertices[edge.ivertices[1]].y;
            edge.jacob_det = ((x2 - x1).powi(2) + (y2 - y1).powi(2)).sqrt();
        }
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
    pub fn compute_circle_circumradius(&mut self) {
        for element in self.elements.iter_mut() {
            let x1 = self.vertices[element.ivertices[0]].x;
            let x2 = self.vertices[element.ivertices[1]].x;
            let x3 = self.vertices[element.ivertices[2]].x;
            let y1 = self.vertices[element.ivertices[0]].y;
            let y2 = self.vertices[element.ivertices[1]].y;
            let y3 = self.vertices[element.ivertices[2]].y;
            let a = ((x2 - x1).powi(2) + (y2 - y1).powi(2)).sqrt();
            let b = ((x3 - x2).powi(2) + (y3 - y2).powi(2)).sqrt();
            let c = ((x1 - x3).powi(2) + (y1 - y3).powi(2)).sqrt();
            let s = (a + b + c) / 2.0;
            element.circumradius = a * b * c / (4.0 * (s * (s - a) * (s - b) * (s - c)).sqrt());
        }
    }
    /*
    pub fn compute_derivatives(&mut self) {
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
            for igp in 0..self.gauss_points.cell_gp_number {
                for ibasis in 0..self.basis.dof {
                    element.derivatives[[igp, ibasis]].insert((1, 0), dxi_dx * self.basis.derivatives[[igp, ibasis]].get(&(1, 1)).unwrap() + deta_dx * self.basis.derivatives[[igp, ibasis]].get(&(0, 2)).unwrap());
                    element.derivatives[[igp, ibasis]].insert((0, 1), dxi_dy * self.basis.derivatives[[igp, ibasis]].get(&(1, 0)).unwrap() + deta_dy * self.basis.derivatives[[igp, ibasis]].get(&(0, 1)).unwrap());
                    element.derivatives[[igp, ibasis]].insert((2, 0), dxi_dx * self.basis.derivatives[[igp, ibasis]].get(&(2, 1)).unwrap() + deta_dx * self.basis.derivatives[[igp, ibasis]].get(&(0, 3)).unwrap());
                }
            }
        }
    }
    */
    pub fn set_neighbours(&mut self) {
        for element in self.elements.iter_mut() {
            for (in_cell_index, iedge) in element.iedges.iter().enumerate() {
                match *iedge {
                    EdgeTypeAndIndex::Boundary(_) => (),
                    EdgeTypeAndIndex::Internal(iedge) => {
                        let edge = &self.edges[iedge];
                        match element.normal_directions[in_cell_index] {
                            NormalDirection::Inward => {
                                element.ineighbours[in_cell_index] = Some(edge.ielements[0]);
                            }
                            NormalDirection::Outward => {
                                element.ineighbours[in_cell_index] = Some(edge.ielements[1]);
                            }
                        }
                    }
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
    pub iedges: Array<EdgeTypeAndIndex, Ix1>,
    pub ineighbours: Array<Option<usize>, Ix1>,
    pub mass_mat_diag: Array<f64, Ix1>,
    pub derivatives: Array<HashMap<(usize, usize), f64>, Ix2>,
    pub normal_directions: Array<NormalDirection, Ix1>,
    pub jacob_det: f64,
    pub circumradius: f64,
}
#[derive(Clone)]
pub enum EdgeTypeAndIndex {
    Boundary(usize),
    Internal(usize),
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
    pub in_cell_indices: [usize; 2],
    //pub hcr: f64,
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