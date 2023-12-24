use ndarray::Ix1;
use ndarray::ArrayView;

pub fn hllc(left_value: ArrayView<f64, Ix1>, right_value: ArrayView<f64, Ix1>, normal: &[f64; 2], hcr: &f64) -> [f64; 4] {
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
    let cl = (hcr * pl / ql[0]).sqrt();
    let cr = (hcr * pr / qr[0]).sqrt();
    let mut fl = [0.0f64; 4];
    let mut fr = [0.0f64; 4];
    fl[0] = ql[1];
    fl[1] = ql[1] * ul + pl * nx;
    fl[2] = ql[2] * ul + pl * ny;
    fl[3] = ul * (ql[3] + pl);
    fr[0] = qr[1];
    fr[1] = qr[1] * ur + pr * nx;
    fr[2] = qr[2] * ur + pr * ny;
    fr[3] = ur * (qr[3] + pr);
    let sl = vnl - cl;
    let sr = vnr + cr;
    if sl >= 0.0 {
        flux[0] = fl[0];
        flux[1] = fl[1];
        flux[2] = fl[2];
        flux[3] = fl[3];
        flux
    } else if sr <= 0.0 {
        flux[0] = fr[0];
        flux[1] = fr[1];
        flux[2] = fr[2];
        flux[3] = fr[3];
        flux
    } else {
        let s_star = (pr - pl + ql[0] * vnl * (sl - vnl) - qr[0] * vnr * (sr - vnr)) / (ql[0] * (sl - vnl) - qr[0] * (sr - vnr));
        if s_star >= 0.0 {
            flux[0] = fl[0] + sl * (ql[0] - ql[0] * (sl - vnl)) / (sl - s_star);
            flux[1] = fl[1] + sl * (ql[1] - ql[1] * (sl - vnl) + (s_star - vnl) * (ql[0] - ql[0] * (sl - vnl))) / (sl - s_star);
            flux[2] = fl[2] + sl * (ql[2] - ql[2] * (sl - vnl)) / (sl - s_star);
            flux[3] = fl[3] + sl * (ql[3] - ql[3] * (sl - vnl) + (s_star - vnl) * (ql[0] - ql[0] * (sl - vnl))) / (sl - s_star);
            flux
        }
        else {
            flux[0] = fr[0] + sr * (qr[0] - qr[0] * (sr - vnr)) / (sr - s_star);
            flux[1] = fr[1] + sr * (qr[1] - qr[1] * (sr - vnr) + (s_star - vnr) * (qr[0] - qr[0] * (sr - vnr))) / (sr - s_star);
            flux[2] = fr[2] + sr * (qr[2] - qr[2] * (sr - vnr)) / (sr - s_star);
            flux[3] = fr[3] + sr * (qr[3] - qr[3] * (sr - vnr) + (s_star - vnr) * (qr[0] - qr[0] * (sr - vnr))) / (sr - s_star);
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