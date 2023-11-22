use ndarray::{Ix1, Ix3};
use ndarray::Array;
use ndarray::ArrayViewMut;
pub mod hweno;
use crate::mesh::{Mesh, Element, Edge, BoundaryEdge};
pub struct Limiter{
    pub limiter: Box<dyn LimiterScheme>,

}
pub trait LimiterScheme {
    fn apply(&self, solutions: ArrayViewMut<f64, Ix3>);
}