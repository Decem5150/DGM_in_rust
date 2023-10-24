pub struct GaussPoints {
    pub cell_points: Vec<[f64; 2]>,
    pub cell_weights: Vec<f64>,
    pub edge_points: [Vec<f64>; 3],
    pub edge_weights: Vec<f64>,
}
impl GaussPoints {
    pub fn new(number_of_points: i64) -> GaussPoints {
        let (cell_points, cell_weights) = gauss_points_triangle(number_of_points);
        let (interval_points, edge_weights) = gauss_points_interval(number_of_points);
        GaussPoints {
            cell_points,
            cell_weights,
            edge_points: [interval_points.iter().map(|&x| ((x + 1.0) / 2.0, 0.0)).collect(), 
                          interval_points.iter().map(|&x| ((1.0 - x) / 2.0, (1.0 + x) / 2.0)).collect(), 
                          interval_points.iter().map(|&x| (0.0, (x + 1.0) / 2.0)).collect()],
            edge_weights,
        }
    }
}
pub fn gauss_points_interval(number_of_points: i64) -> (Vec<f64>, Vec<f64>) {
    let (gauss_points, gauss_weights) = match number_of_points {
        3 => {
            let points: Vec<f64> = vec!{-0.7745966692414834, 0.0, 0.7745966692414834};
            let weights: Vec<f64> = vec!{0.5555555555555556, 0.8888888888888889, 0.5555555555555556};
            (points, weights) 
        }
        4 => {
            let points: Vec<f64> = vec!{-0.8611363115940526, -0.3399810435848563, 0.3399810435848563, 0.8611363115940526};
            let weights: Vec<f64> = vec!{0.3478548451374538, 0.6521451548625461, 0.6521451548625461, 0.3478548451374538};
            (points, weights)
        }
        5 => {
            let points: Vec<f64> = vec!{-0.9061798459386640, -0.5384693101056831, 0.0, 0.5384693101056831, 0.9061798459386640};
            let weights: Vec<f64> = vec!{0.2369268850561891, 0.4786286704993665, 0.5688888888888889, 0.4786286704993665, 0.2369268850561891};
            (points, weights)
        }
        _ => {
            panic!("Number of points not supported");
        }        
    };
    (gauss_points, gauss_weights)
}
pub fn gauss_points_triangle(number_of_points: i64) -> (Vec<[f64; 2]>, Vec<f64>) {
    let (gauss_points, gauss_weights) = match number_of_points {
        4 => {
            let points: Vec<[f64; 2]> = vec![[0.33333333333, 0.33333333333], [0.2, 0.6], [0.2, 0.2], [0.6, 0.2]];
            let weights: Vec<f64> = vec![-0.28125, 0.26041666667, 0.26041666667, 0.26041666667];
            (points, weights)
        }
        _ => {
            panic!("Number of points not supported");
        }     
    };
    (gauss_points, gauss_weights)
}