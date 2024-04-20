use std::collections::HashMap;
use ndarray::{Array, Ix3};
use ndarray::{Ix1, Ix2};
use crate::gauss_point::GaussPoints;
use crate::basis_function::DubinerBasis;
pub struct Mesh { 
    pub elements: Array<Element, Ix1>,
    pub edges: Array<Edge, Ix1>,
    pub internal_edge_indices: Array<usize, Ix1>,
    pub internal_element_indices: Array<usize, Ix1>,
    pub boundary_element_indices: Array<usize, Ix1>,
    pub vertices: Array<Vertex, Ix1>,
    pub patches: Array<Patch, Ix1>,
}
impl Mesh {
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
    pub fn compute_circumradius(&mut self) {
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
    pub fn compute_minimum_height(&mut self) {
        for element in self.elements.iter_mut() {
            let surface = 0.5 * element.jacob_det;
            let a = {
                let x1 = self.vertices[element.ivertices[0]].x;
                let x2 = self.vertices[element.ivertices[1]].x;
                let y1 = self.vertices[element.ivertices[0]].y;
                let y2 = self.vertices[element.ivertices[1]].y;
                ((x2 - x1).powi(2) + (y2 - y1).powi(2)).sqrt()
            };
            let b = {
                let x2 = self.vertices[element.ivertices[1]].x;
                let x3 = self.vertices[element.ivertices[2]].x;
                let y2 = self.vertices[element.ivertices[1]].y;
                let y3 = self.vertices[element.ivertices[2]].y;
                ((x3 - x2).powi(2) + (y3 - y2).powi(2)).sqrt()
            };
            let c = {
                let x3 = self.vertices[element.ivertices[2]].x;
                let x1 = self.vertices[element.ivertices[0]].x;
                let y3 = self.vertices[element.ivertices[2]].y;
                let y1 = self.vertices[element.ivertices[0]].y;
                ((x1 - x3).powi(2) + (y1 - y3).powi(2)).sqrt()
            };
            element.minimum_height = (2.0 * surface / a).min(2.0 * surface / b).min(2.0 * surface / c);
        }
    }
    pub fn compute_dphi(&mut self, basis: &DubinerBasis, cell_ngp: usize, edge_ngp: usize, nbasis: usize) {
        for element in self.elements.iter_mut() {
            let x1 = self.vertices[element.ivertices[0]].x;
            let x2 = self.vertices[element.ivertices[1]].x;
            let x3 = self.vertices[element.ivertices[2]].x;
            let y1 = self.vertices[element.ivertices[0]].y;
            let y2 = self.vertices[element.ivertices[1]].y;
            let y3 = self.vertices[element.ivertices[2]].y;
            let dx_dxi = x2 - x1;
            let dx_deta = x3 - x1;
            let dy_dxi = y2 - y1;
            let dy_deta = y3 - y1;
            let dxi_dx = dy_deta / element.jacob_det;
            let dxi_dy = -dx_deta / element.jacob_det;
            let deta_dx = -dy_dxi / element.jacob_det;
            let deta_dy = dx_dxi / element.jacob_det;
            for igp in 0..cell_ngp {
                for ibasis in 0..nbasis {
                    let dphi_dxi = *basis.dphis_cell_gps[[igp, ibasis]].get(&(1, 0)).unwrap();
                    let dphi_deta = *basis.dphis_cell_gps[[igp, ibasis]].get(&(0, 1)).unwrap();
                    let dphi_dxidxi = *basis.dphis_cell_gps[[igp, ibasis]].get(&(2, 0)).unwrap();
                    let dphi_dxideta = *basis.dphis_cell_gps[[igp, ibasis]].get(&(1, 1)).unwrap();
                    let dphi_detadeta = *basis.dphis_cell_gps[[igp, ibasis]].get(&(0, 2)).unwrap();
                    element.dphis_cell_gps[[igp, ibasis]].insert((1, 0), dxi_dx * dphi_dxi + deta_dx * dphi_deta);
                    element.dphis_cell_gps[[igp, ibasis]].insert((0, 1), dxi_dy * dphi_dxi + deta_dy * dphi_deta);
                    element.dphis_cell_gps[[igp, ibasis]].insert((2, 0), dxi_dx.powf(2.0) * dphi_dxidxi + 2.0 * dxi_dx * deta_dx * dphi_dxideta + deta_dx.powf(2.0) * dphi_detadeta);
                    element.dphis_cell_gps[[igp, ibasis]].insert((1, 1), dxi_dx * dxi_dy * dphi_dxidxi + (dxi_dx * deta_dy + deta_dx * dxi_dy) * dphi_dxideta + deta_dx * deta_dy * dphi_detadeta);
                    element.dphis_cell_gps[[igp, ibasis]].insert((0, 2), dxi_dy.powf(2.0) * dphi_dxidxi + 2.0 * dxi_dy * deta_dy * dphi_dxideta + deta_dy.powf(2.0) * dphi_detadeta);
                }
            }
            for iedge in 0..3 {
                for igp in 0..edge_ngp {
                    for ibasis in 0..nbasis {
                        let dphi_dxi = *basis.dphis_edge_gps[[iedge, igp, ibasis]].get(&(1, 0)).unwrap();
                        let dphi_deta = *basis.dphis_edge_gps[[iedge, igp, ibasis]].get(&(0, 1)).unwrap();
                        let dphi_dxidxi = *basis.dphis_edge_gps[[iedge, igp, ibasis]].get(&(2, 0)).unwrap();
                        let dphi_dxideta = *basis.dphis_edge_gps[[iedge, igp, ibasis]].get(&(1, 1)).unwrap();
                        let dphi_detadeta = *basis.dphis_edge_gps[[iedge, igp, ibasis]].get(&(0, 2)).unwrap();
                        element.dphis_edge_gps[[iedge, igp, ibasis]].insert((1, 0), dxi_dx * dphi_dxi + deta_dx * dphi_deta);
                        element.dphis_edge_gps[[iedge, igp, ibasis]].insert((0, 1), dxi_dy * dphi_dxi + deta_dy * dphi_deta);
                        element.dphis_edge_gps[[iedge, igp, ibasis]].insert((2, 0), dxi_dx.powf(2.0) * dphi_dxidxi + 2.0 * dxi_dx * deta_dx * dphi_dxideta + deta_dx.powf(2.0) * dphi_detadeta);
                        element.dphis_edge_gps[[iedge, igp, ibasis]].insert((1, 1), dxi_dx * dxi_dy * dphi_dxidxi + (dxi_dx * deta_dy + deta_dx * dxi_dy) * dphi_dxideta + deta_dx * deta_dy * dphi_detadeta);
                        element.dphis_edge_gps[[iedge, igp, ibasis]].insert((0, 2), dxi_dy.powf(2.0) * dphi_dxidxi + 2.0 * dxi_dy * deta_dy * dphi_dxideta + deta_dy.powf(2.0) * dphi_detadeta);
                    }
                }
            }
        }
    }
    pub fn set_neighbours(&mut self) {
        for element in self.elements.iter_mut() {
            for (in_cell_index, iedge) in element.iedges.iter().enumerate() {
                let edge = &self.edges[*iedge];
                if edge.ielements[1] == -1 {
                    continue;
                }
                else if element.ivertices[in_cell_index] == edge.ivertices[0] {
                    element.ineighbours[in_cell_index] = edge.ielements[1];
                } else {
                    element.ineighbours[in_cell_index] = edge.ielements[0];
                }
            }
        }
    }
}
pub struct Vertex {
    pub x: f64,
    pub y: f64,
}
#[derive(Debug)]
pub struct Element {
    pub ivertices: Array<usize, Ix1>,
    pub iedges: Array<usize, Ix1>,
    pub ineighbours: Array<isize, Ix1>,
    pub dphis_cell_gps: Array<HashMap<(usize, usize), f64>, Ix2>,
    pub dphis_edge_gps: Array<HashMap<(usize, usize), f64>, Ix3>,
    pub jacob_det: f64,
    pub circumradius: f64,
    pub minimum_height: f64,
}
#[derive(Clone, Debug)]
pub enum EdgeTypeAndIndex {
    Boundary(usize),
    Internal(usize),
}
#[derive(Debug)]
pub struct Edge {
    pub ielements: [isize; 2],
    pub ivertices: [usize; 2],
    pub jacob_det: f64,
    pub normal: [f64; 2],
    pub in_cell_indices: [isize; 2],
    pub ipatch: isize
    //pub hcr: f64,
}
#[derive(Clone)]
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