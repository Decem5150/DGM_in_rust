use ndarray::{Ix1, Ix2, Ix3};
use ndarray::Array;
use ndarray::ArrayViewMut;
use super::LimiterScheme;
use crate::mesh::{Mesh, Element, Edge, BoundaryEdge};
pub struct HWENO<'a> {
    detectors: Array<f64, Ix2>, //Density and total energy-based
    mesh: &'a Mesh<'a>,
}
impl<'a> HWENO<'a> {
    fn kxrcf(&self, solutions: ArrayViewMut<f64, Ix3>) {
        for ielem in 0..self.mesh.elements.len() {

        }
    }
}
impl<'a> LimiterScheme for HWENO<'a> {
    fn apply(&self, solutions: ArrayViewMut<f64, Ix3>){
        
    }
}