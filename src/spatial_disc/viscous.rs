use ndarray::{Array, ArrayView, Ix1, Ix2};

use super::SpatialDisc;

impl<'a> SpatialDisc<'a> {
    pub fn compute_dynamic_viscosity(temperature: f64) -> f64 {
        // Sutherland's law
        let mu_ref = 1.716e-5;
        let sutherland_constant = 110.4;
        let t_ref = 273.15;
        mu_ref * (temperature / t_ref).powf(1.5) * (t_ref + sutherland_constant) / (temperature + sutherland_constant)
    }
    pub fn viscous_jacobian(&self, q: &Array<f64, Ix1>, mu: f64) -> (Array<f64, Ix2>, Array<f64, Ix2>, Array<f64, Ix2>, Array<f64, Ix2>) {
        let neq = self.solver_param.number_of_equations;
        let mut g11: Array<f64, Ix2>= Array::zeros((neq, neq));
        let mut g12: Array<f64, Ix2>= Array::zeros((neq, neq));
        let mut g21: Array<f64, Ix2>= Array::zeros((neq, neq));
        let mut g22: Array<f64, Ix2>= Array::zeros((neq, neq));
        let hcr = self.flow_param.hcr;
        let pr = self.flow_param.prandtl_number;
        let u = q[1] / q[0];
        let v = q[2] / q[0];
        let e = q[3] / q[0];

        g11[[1, 0]] = - 4.0 / 3.0 * u;
        g11[[1, 1]] = 4.0 / 3.0;
        g11[[2, 0]] = - v;
        g11[[2, 2]] = 1.0;
        g11[[3, 0]] = - (4.0 / 3.0 * u.powf(2.0) + v.powf(2.0) + hcr / pr * (e - u.powf(2.0) - v.powf(2.0)));
        g11[[3, 1]] = (4.0 / 3.0 - hcr / pr) * u;
        g11[[3, 2]] = (1.0 - hcr / pr) * v;
        g11[[3, 3]] = hcr / pr;

        g12[[1, 0]] = 2.0 / 3.0 * v;
        g12[[1, 2]] = - 2.0 / 3.0;
        g12[[2, 0]] = - u;
        g12[[2, 1]] = 1.0;
        g12[[3, 0]] = - 1.0 / 3.0 * u * v;
        g12[[3, 1]] = v;
        g12[[3, 2]] = - 2.0 / 3.0 * u;

        g21[[1, 0]] = - v;
        g21[[1, 2]] = 1.0;
        g21[[2, 0]] = 2.0 / 3.0 * u;
        g21[[2, 1]] = - 2.0 / 3.0;
        g21[[3, 0]] = - 1.0 / 3.0 * u * v;
        g21[[3, 1]] = - 2.0 / 3.0 * v;
        g21[[3, 2]] = u;

        g22[[1, 0]] = - u;
        g22[[1, 1]] = 1.0;
        g22[[2, 0]] = - 4.0 / 3.0 * v;
        g22[[2, 2]] = 4.0 / 3.0;
        g22[[3, 0]] = - (u.powf(2.0) + 4.0 / 3.0 * v.powf(2.0) + hcr / pr * (e - u.powf(2.0) - v.powf(2.0)));
        g22[[3, 1]] = (1.0 - hcr / pr) * u;
        g22[[3, 2]] = (4.0 / 3.0 - hcr / pr) * v;
        g22[[3, 3]] = hcr / pr;

        g11 *= mu / q[0];
        g12 *= mu / q[0];
        g21 *= mu / q[0];
        g22 *= mu / q[0];

        (g11, g12, g21, g22)
    }
}