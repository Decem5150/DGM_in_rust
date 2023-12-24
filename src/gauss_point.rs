use ndarray::{Array, stack, Axis};
use ndarray::{Ix1, Ix2, Ix3};
use ndarray::array;
pub struct GaussPoints {
    pub cell_gp_number: usize,
    pub edge_gp_number: usize,
    pub cell_points: Array<f64, Ix2>,
    pub cell_weights: Array<f64, Ix1>,
    pub edge_points: Array<f64, Ix3>,
    pub edge_weights: Array<f64, Ix1>,
}
impl GaussPoints {
    pub fn new(cell_gp: usize, edge_gp: usize) -> GaussPoints {
        let (cell_points, cell_weights) = gauss_points_triangle(cell_gp);
        let (interval_points, edge_weights) = gauss_points_interval(edge_gp);

        let ip11 = interval_points.mapv(|x| (x + 1.0) / 2.0);
        let ip12 = interval_points.mapv(|x| 0.0);
        let ip1 = stack(Axis(0), &[ip11.view(), ip12.view()]).unwrap().t().to_owned();
        let ip21 = interval_points.mapv(|x| (1.0 - x) / 2.0);
        let ip22 = interval_points.mapv(|x| (1.0 + x) / 2.0);
        let ip2 = stack(Axis(0), &[ip21.view(), ip22.view()]).unwrap().t().to_owned();
        let ip31 = interval_points.mapv(|x| 0.0);
        let ip32 = interval_points.mapv(|x| (x + 1.0) / 2.0);
        let ip3 = stack(Axis(0), &[ip31.view(), ip32.view()]).unwrap().t().to_owned();
        
        let edge_points = stack(Axis(0), &[ip1.view(), ip2.view(), ip3.view()]).unwrap();
        
        GaussPoints {
            cell_gp_number: cell_gp,
            edge_gp_number: edge_gp,
            cell_points,
            cell_weights,
            edge_points,
            edge_weights,
        }
    }
}
pub fn gauss_points_interval(number_of_points: usize) -> (Array<f64, Ix1>, Array<f64, Ix1>) {
    let (gauss_points, gauss_weights) = match number_of_points {
        3 => {
            let points = array![-0.7745966692414834, 0.0, 0.7745966692414834];
            let weights = array![0.5555555555555556, 0.8888888888888889, 0.5555555555555556];
            (points, weights) 
        }
        4 => {
            let points = array![-0.8611363115940526, -0.3399810435848563, 0.3399810435848563, 0.8611363115940526];
            let weights = array![0.3478548451374538, 0.6521451548625461, 0.6521451548625461, 0.3478548451374538];
            (points, weights)
        }
        5 => {
            let points = array![-0.9061798459386640, -0.5384693101056831, 0.0, 0.5384693101056831, 0.9061798459386640];
            let weights = array![0.2369268850561891, 0.4786286704993665, 0.5688888888888889, 0.4786286704993665, 0.2369268850561891];
            (points, weights)
        }
        _ => {
            panic!("Number of points not supported");
        }        
    };
    (gauss_points, gauss_weights)
}
pub fn gauss_points_triangle(number_of_points: usize) -> (Array<f64, Ix2>, Array<f64, Ix1>) {
    let (gauss_points, gauss_weights) = match number_of_points {
        4 => {
            let points = array![[0.33333333333, 0.33333333333], [0.2, 0.6], [0.2, 0.2], [0.6, 0.2]];
            let weights = array![-0.28125, 0.26041666667, 0.26041666667, 0.26041666667];
            (points, weights)
        }
        _ => {
            panic!("Number of points not supported");
        }     
    };
    (gauss_points, gauss_weights)
}