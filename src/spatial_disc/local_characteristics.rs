use ndarray::Array;
use ndarray::{ArrayView, ArrayViewMut};
use ndarray::{Ix1, Ix2};
pub fn compute_eigenmatrix(q: ArrayView<f64, Ix1>, nx: f64, ny: f64, hcr: f64) -> (Array<f64, Ix2>, Array<f64, Ix2>) {
    let u = q[1] / q[0];
    let v = q[2] / q[0];
    let p = (hcr - 1.0) * (q[3] - 0.5 * (q[1] * q[1] + q[2] * q[2]) / q[0]);
    let c = (hcr * p / q[0]).sqrt();
    let h = (q[3] + p) / q[0];
    let b1 = (hcr - 1.0) / c.powf(2.0);
    let b2 = b1 * (u.powf(2.0) + v.powf(2.0)) / 2.0;
    let mut lmatrix: Array<f64, Ix2> = Array::zeros((4, 4));
    let mut rmatrix: Array<f64, Ix2> = Array::zeros((4, 4));

    lmatrix[[0, 0]] = (b2 + (u * nx + v * ny) / c) / 2.0;
    lmatrix[[0, 1]] = -(b1 * u + nx / c) / 2.0;
    lmatrix[[0, 2]] = -(b1 * v + ny / c) / 2.0;
    lmatrix[[0, 3]] = b1 / 2.0;
    lmatrix[[1, 0]] = ny * u - nx * v;
    lmatrix[[1, 1]] = -ny;
    lmatrix[[1, 2]] = nx;
    lmatrix[[1, 3]] = 0.0;
    lmatrix[[2, 0]] = 1.0 - b2;
    lmatrix[[2, 1]] = b1 * u;
    lmatrix[[2, 2]] = b1 * v;
    lmatrix[[2, 3]] = -b1;
    lmatrix[[3, 0]] = (b2 - (u * nx + v * ny) / c) / 2.0;
    lmatrix[[3, 1]] = -(b1 * u - nx / c) / 2.0;
    lmatrix[[3, 2]] = -(b1 * v - ny / c) / 2.0;
    lmatrix[[3, 3]] = b1 / 2.0;

    rmatrix[[0, 0]] = 1.0;
    rmatrix[[0, 1]] = 0.0;
    rmatrix[[0, 2]] = 1.0;
    rmatrix[[0, 3]] = 1.0;
    rmatrix[[1, 0]] = u - c * nx;
    rmatrix[[1, 1]] = -ny;
    rmatrix[[1, 2]] = u;
    rmatrix[[1, 3]] = u + c * nx;
    rmatrix[[2, 0]] = v - c * ny;
    rmatrix[[2, 1]] = nx;
    rmatrix[[2, 2]] = v;
    rmatrix[[2, 3]] = v + c * ny;
    rmatrix[[3, 0]] = h - c * (u * nx + v * ny);
    rmatrix[[3, 1]] = - u * ny + v * nx;
    rmatrix[[3, 2]] = (u.powf(2.0) + v.powf(2.0)) / 2.0;
    rmatrix[[3, 3]] = h + c * (u * nx + v * ny);

    (lmatrix, rmatrix)
}

    