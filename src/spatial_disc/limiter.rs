use ndarray::{Ix1, Ix3};
use ndarray::Array;
use ndarray::{ArrayView, ArrayViewMut};
pub mod weno;
use crate::mesh::{Mesh, Element, Edge, BoundaryEdge};
use crate::solver::{SolverParameters, FlowParameters};
use crate::basis_function::{DubinerBasis, GaussPoints};
pub enum LimiterType {
    WENO,
    None,
}
pub enum DetectorType {
    TVB,
    KXRCF,
    None,
}
pub struct Limiter<'a> {
    pub limiter_type: LimiterType,
    pub detector_type: DetectorType,
    pub mesh: &'a Mesh<'a>,
    pub is_troubled: Array<bool, Ix1>,
    pub basis: &'a DubinerBasis<'a>,
    pub gauss_points: &'a GaussPoints,
    pub solver_param: &'a SolverParameters,
    pub flow_param: &'a FlowParameters,
}
impl<'a> Limiter<'a> {
    pub fn apply_limiter(&self, solutions: ArrayViewMut<f64, Ix3>) {
        match self.limiter_type {
            LimiterType::HWENO => self.hweno(solutions),
            LimiterType::WENO => panic!("WENO limiter not implemented"),
            LimiterType::None => (),  
        }
    }
}