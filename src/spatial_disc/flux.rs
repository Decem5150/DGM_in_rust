use crate::solver::ConsVar;
pub trait InvisFluxScheme<'a> {
    fn compute(&self, left_value: ConsVar, right_value: ConsVar, flux: &'a mut ConsVar, normal: [f64; 2], hcr: f64);
}
pub struct HLLC;
impl<'a> InvisFluxScheme<'a> for HLLC {
    fn compute(&self, left_value: ConsVar, right_value: ConsVar, flux: &'a mut ConsVar, normal: [f64; 2], hcr: f64) {
        let nx = normal[0];
        let ny = normal[1];
        let q1l = left_value.density;
        let q2l = left_value.x_momentum;
        let q3l = left_value.y_momentum;
        let q4l = left_value.energy;
        let q1r = right_value.density;
        let q2r = right_value.x_momentum;
        let q3r = right_value.y_momentum;
        let q4r = right_value.energy;
        let ul = q2l / q1l;
        let ur = q2r / q1r;
        let vl = q3l / q1l;
        let vr = q3r / q1r;
        let vnl = ul * nx + vl * ny;
        let vnr = ur * nx + vr * ny;
        let p1l = (hcr - 1.0) * (q4l - 0.5 * (q2l * q2l + q3l * q3l) / q1l);
        let p1r = (hcr - 1.0) * (q4r - 0.5 * (q2r * q2r + q3r * q3r) / q1r);
        let c1l = (hcr * p1l / q1l).sqrt();
        let c1r = (hcr * p1r / q1r).sqrt();
        let f1l = q2l;
        let f2l = q2l * vnl + p1l * nx;
        let f3l = q3l * vnl + p1l * ny;
        let f4l = vnl * (q4l + p1l);
        let f1r = q2r;
        let f2r = q2r * vnr + p1r * nx;
        let f3r = q3r * vnr + p1r * ny;
        let f4r = vnr * (q4r + p1r);
        let sl = vnl - c1l;
        let sr = vnr + c1r;
        if sl >= 0.0 {
            flux.density_flux = f1l;
            flux.x_momentum_flux = f2l;
            flux.y_momentum_flux = f3l;
            flux.energy_flux = f4l;
            return;
        } else if sr <= 0.0 {
            flux.density_flux = f1r;
            flux.x_momentum_flux = f2r;
            flux.y_momentum_flux = f3r;
            flux.energy_flux = f4r;
            return;
        } else {
            let s_star = (p1r - p1l + q1l * vnl * (sl - vnl) - q1r * vnr * (sr - vnr)) / (q1l * (sl - vnl) - q1r * (sr - vnr));
            if s_star >= 0.0 {
                flux.density_flux = f1l + sl * (q1l - q1l * (sl - vnl)) / (sl - s_star);
                flux.x_momentum_flux = f2l + sl * (q2l - q2l * (sl - vnl) + (s_star - vnl) * (q1l - q1l * (sl - vnl))) / (sl - s_star);
                flux.y_momentum_flux = f3l + sl * (q3l - q3l * (sl - vnl)) / (sl - s_star);
                flux.energy_flux = f4l + sl * (q4l - q4l * (sl - vnl) + (s_star - vnl) * (q1l - q1l * (sl - vnl))) / (sl - s_star);
                return;
            }
            else {
                flux.density_flux = f1r + sr * (q1r - q1r * (sr - vnr)) / (sr - s_star);
                flux.x_momentum_flux = f2r + sr * (q2r - q2r * (sr - vnr) + (s_star - vnr) * (q1r - q1r * (sr - vnr))) / (sr - s_star);
                flux.y_momentum_flux = f3r + sr * (q3r - q3r * (sr - vnr)) / (sr - s_star);
                flux.energy_flux = f4r + sr * (q4r - q4r * (sr - vnr) + (s_star - vnr) * (q1r - q1r * (sr - vnr))) / (sr - s_star);
                return;
            }
        } 
    }
}
pub fn flux(q1: f64, q2: f64, q3: f64, q4: f64, hcr: f64) -> (f64, f64, f64, f64, f64, f64, f64, f64) {
    let hcr: f64 = 1.4;
    let u = q2 / q1;
    let v = q3 / q1;
    let p = (hcr - 1.0) * (q4 - 0.5 * (q2 * q2 + q3 * q3) / q1);
    let f1 = q2;
    let f2 = q2 * u + p;
    let f3 = q3 * u;
    let f4 = u * (q4 + p);
    let g1 = q3;
    let g2 = q2 * v;
    let g3 = q3 * v + p;
    let g4 = v * (q4 + p);
    (f1, f2, f3, f4, g1, g2, g3, g4)
}
pub fn hllc_flux<'a>(left_value: ConsVar, right_value: ConsVar, flux: &'a mut ConsVar, normal: [f64; 2]) {
    let nx = normal[0];
    let ny = normal[1];
    let hcr: f64 = 1.4;
    let q1l = left_value.density;
    let q2l = left_value.x_momentum;
    let q3l = left_value.y_momentum;
    let q4l = left_value.energy;
    let q1r = right_value.density;
    let q2r = right_value.x_momentum;
    let q3r = right_value.y_momentum;
    let q4r = right_value.energy;
    let ul = q2l / q1l;
    let ur = q2r / q1r;
    let vl = q3l / q1l;
    let vr = q3r / q1r;
    let vnl = ul * nx + vl * ny;
    let vnr = ur * nx + vr * ny;
    let p1l = (hcr - 1.0) * (q4l - 0.5 * (q2l * q2l + q3l * q3l) / q1l);
    let p1r = (hcr - 1.0) * (q4r - 0.5 * (q2r * q2r + q3r * q3r) / q1r);
    let c1l = (hcr * p1l / q1l).sqrt();
    let c1r = (hcr * p1r / q1r).sqrt();
    let f1l = q2l;
    let f2l = q2l * vnl + p1l * nx;
    let f3l = q3l * vnl + p1l * ny;
    let f4l = vnl * (q4l + p1l);
    let f1r = q2r;
    let f2r = q2r * vnr + p1r * nx;
    let f3r = q3r * vnr + p1r * ny;
    let f4r = vnr * (q4r + p1r);
    let sl = vnl - c1l;
    let sr = vnr + c1r;
    if sl >= 0.0 {
        flux.density_flux = f1l;
        flux.x_momentum_flux = f2l;
        flux.y_momentum_flux = f3l;
        flux.energy_flux = f4l;
        return;
    } else if sr <= 0.0 {
        flux.density_flux = f1r;
        flux.x_momentum_flux = f2r;
        flux.y_momentum_flux = f3r;
        flux.energy_flux = f4r;
        return;
    } else {
        let s_star = (p1r - p1l + q1l * vnl * (sl - vnl) - q1r * vnr * (sr - vnr)) / (q1l * (sl - vnl) - q1r * (sr - vnr));
        if s_star >= 0.0 {
            flux.density_flux = f1l + sl * (q1l - q1l * (sl - vnl)) / (sl - s_star);
            flux.x_momentum_flux = f2l + sl * (q2l - q2l * (sl - vnl) + (s_star - vnl) * (q1l - q1l * (sl - vnl))) / (sl - s_star);
            flux.y_momentum_flux = f3l + sl * (q3l - q3l * (sl - vnl)) / (sl - s_star);
            flux.energy_flux = f4l + sl * (q4l - q4l * (sl - vnl) + (s_star - vnl) * (q1l - q1l * (sl - vnl))) / (sl - s_star);
            return;
        }
        else {
            flux.density_flux = f1r + sr * (q1r - q1r * (sr - vnr)) / (sr - s_star);
            flux.x_momentum_flux = f2r + sr * (q2r - q2r * (sr - vnr) + (s_star - vnr) * (q1r - q1r * (sr - vnr))) / (sr - s_star);
            flux.y_momentum_flux = f3r + sr * (q3r - q3r * (sr - vnr)) / (sr - s_star);
            flux.energy_flux = f4r + sr * (q4r - q4r * (sr - vnr) + (s_star - vnr) * (q1r - q1r * (sr - vnr))) / (sr - s_star);
            return;
        }
    }
}
