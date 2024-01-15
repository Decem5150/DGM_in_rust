use ndarray::Ix1;
use ndarray::ArrayView;
/*
pub fn hr_slau2(left_value: ArrayView<f64, Ix1>, right_value: ArrayView<f64, Ix1>, normal: &[f64; 2], hcr: &f64) -> [f64; 4] {
    let mut flux = [0.0f64; 4];
    let nx = normal[0];
    let ny = normal[1];
    let ql = left_value;
    let qr = right_value;
    let ul = ql[1] / ql[0];
    let ur = qr[1] / qr[0];
    let vl = ql[2] / ql[0];
    let vr = qr[2] / qr[0];
    let vnl = ul * nx + vl * ny;
    let vnr = ur * nx + vr * ny;
    let pl = (hcr - 1.0) * (ql[3] - 0.5 * (ql[1] * ql[1] + ql[2] * ql[2]) / ql[0]);
    let pr = (hcr - 1.0) * (qr[3] - 0.5 * (qr[1] * qr[1] + qr[2] * qr[2]) / qr[0]);
    let vnl = ul * nx + vl * ny;
    let vnr = ur * nx + vr * ny;
    let c_ave = 0.5 * ((hcr * pl / ql[0]).sqrt() + (hcr * pr / qr[0]).sqrt());
    let mach_left = vnl / c_ave;
    let mach_right = vnr / c_ave;
    let mach_tmp = (1.0 / c_ave * (0.5 * (ul * ul + ur * ur + vl * vl + vr * vr)).sqrt()).min(1.0);
    let chi = (1.0 - mach_tmp).powf(2.0);
}
*/
pub fn hllc(left_value: ArrayView<f64, Ix1>, right_value: ArrayView<f64, Ix1>, normal: &[f64; 2], hcr: &f64) -> [f64; 4] {
    let mut flux = [0.0_f64; 4];
    let nx = normal[0];
    let ny = normal[1];
    let ql = [
        left_value[0], 
        left_value[1] * nx + left_value[2] * ny, 
        - left_value[1] * ny + left_value[2] * nx, 
        left_value[3]
        ];
    let qr = [
        right_value[0], 
        right_value[1] * nx + right_value[2] * ny, 
        - right_value[1] * ny + right_value[2] * nx, 
        right_value[3]
        ];
    let ul = ql[1] / ql[0];
    let ur = qr[1] / qr[0];
    let vl = ql[2] / ql[0];
    let vr = qr[2] / qr[0];
    /* 
    let vnl = ul * nx + vl * ny;
    let vnr = ur * nx + vr * ny;
    let vtl = -ul * ny + vl * nx;
    let vtr = -ur * ny + vr * nx;
    */
    let pl = (hcr - 1.0) * (ql[3] - 0.5 * (ql[1] * ql[1] + ql[2] * ql[2]) / ql[0]);
    let pr = (hcr - 1.0) * (qr[3] - 0.5 * (qr[1] * qr[1] + qr[2] * qr[2]) / qr[0]);
    let cl = (hcr * pl / ql[0]).sqrt();
    let cr = (hcr * pr / qr[0]).sqrt();
    let mut fl = [0.0f64; 4];
    let mut fr = [0.0f64; 4];
    fl[0] = ql[1];
    fl[1] = ql[1] * ul + pl;
    fl[2] = ql[2] * ul;
    fl[3] = ul * (ql[3] + pl);
    fr[0] = qr[1];
    fr[1] = qr[1] * ur + pr;
    fr[2] = qr[2] * ur;
    fr[3] = ur * (qr[3] + pr);
    let p_star = {
        let zeta = (hcr - 1.0) / (2.0 * hcr);
        ((cl + cr - (zeta - 1.0) / 2.0 * (ur - ul)) / (cl / pl.powf(zeta) + cr / pr.powf(zeta))).powf(1.0 / zeta)
    };
    let sl = {
        let ql_lower = {
            if p_star <= pl {
                1.0
            }
            else {
                (1.0 + (hcr + 1.0) / (2.0 * hcr) * (p_star / pl - 1.0)).sqrt()
            }
        };
        ul - cl * ql_lower
    };
    let sr = {
        let qr_lower = {
            if p_star <= pr {
                1.0
            }
            else {
                (1.0 + (hcr + 1.0) / (2.0 * hcr) * (p_star / pr - 1.0)).sqrt()
            }
        };
        ur + cr * qr_lower
    };
    if sl >= 0.0 {
        flux[0] = fl[0];
        flux[1] = fl[1] * nx - fl[2] * ny;
        flux[2] = fl[1] * ny + fl[2] * nx;
        flux[3] = fl[3];
        for &value in flux.iter() {
            if value.is_nan() {
                panic!("Found NaN value in num_flux!");
            }
        }
        flux
    } else if sr <= 0.0 {
        flux[0] = fr[0];
        flux[1] = fr[1] * nx - fr[2] * ny;
        flux[2] = fr[1] * ny + fr[2] * nx;
        flux[3] = fr[3];
        for &value in flux.iter() {
            if value.is_nan() {
                panic!("Found NaN value in num_flux!");
            }
        }
        flux
    } else {
        let s_star = (pr - pl + ql[0] * ul * (sl - ul) - qr[0] * ur * (sr - ur)) / (ql[0] * (sl - ul) - qr[0] * (sr - ur));
        if s_star >= 0.0 {
            let mut q_star = [0.0_f64; 4];
            q_star[0] = ql[0] * (sl - ul) / (sl - s_star);
            q_star[1] = ql[0] * (sl - ul) / (sl - s_star) * s_star;
            q_star[2] = ql[0] * (sl - ul) / (sl - s_star) * vl;
            q_star[3] = ql[0] * (sl - ul) / (sl - s_star) * (ql[3] / ql[0] + (s_star - ul) * (s_star + pl / (ql[0] * (sl - ul))));

            flux[0] = fl[0] + sl * (q_star[0] - ql[0]);
            (flux[1], flux[2]) = {
                let flux_normal = fl[1] + sl * (q_star[1] - ql[1]);
                let flux_tangent = fl[2] + sl * (q_star[2] - ql[2]);
                (flux_normal * nx - flux_tangent * ny, flux_normal * ny + flux_tangent * nx)
            };
            flux[3] = fl[3] + sl * (q_star[3] - ql[3]);
            for &value in flux.iter() {
                if value.is_nan() {
                    panic!("Found NaN value in num_flux!");
                }
            }
            flux
        }
        else {
            let mut q_star = [0.0_f64; 4];
            q_star[0] = qr[0] * (sr - ur) / (sr - s_star);
            q_star[1] = qr[0] * (sr - ur) / (sr - s_star) * s_star;
            q_star[2] = qr[0] * (sr - ur) / (sr - s_star) * vr;
            q_star[3] = qr[0] * (sr - ur) / (sr - s_star) * (qr[3] / qr[0] + (s_star - ur) * (s_star + pr / (qr[0] * (sr - ur))));

            flux[0] = fr[0] + sr * (q_star[0] - qr[0]);
            (flux[1], flux[2]) = {
                let flux_normal = fr[1] + sr * (q_star[1] - qr[1]);
                let flux_tangent = fr[2] + sr * (q_star[2] - qr[2]);
                (flux_normal * nx - flux_tangent * ny, flux_normal * ny + flux_tangent * nx)
            };
            flux[3] = fr[3] + sr * (q_star[3] - qr[3]);
            for &value in flux.iter() {
                if value.is_nan() {
                    panic!("Found NaN value in num_flux!");
                }
            }
            flux
        }
    }
}
pub fn flux(q: ArrayView<f64, Ix1>, hcr: f64) -> ([f64; 4], [f64; 4]) {
    let mut f = [0.0f64; 4];
    let mut g = [0.0f64; 4];
    let u = q[1] / q[0];
    let v = q[2] / q[0];
    let p = (hcr - 1.0) * (q[3] - 0.5 * (q[1] * q[1] + q[2] * q[2]) / q[0]);
    f[0] = q[1];
    f[1] = q[1] * u + p;
    f[2] = q[2] * u;
    f[3] = u * (q[3] + p);
    g[0] = q[2];
    g[1] = q[1] * v;
    g[2] = q[2] * v + p;
    g[3] = v * (q[3] + p);
    (f, g)
}