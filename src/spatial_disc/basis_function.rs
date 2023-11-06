pub use super::gauss_point::GaussPoints;
pub struct DubinerBasis<'a> {
    pub dof: usize,
    pub phi_cell_gp: Vec<Vec<f64>>,
    pub phi_edge_gp: [Vec<Vec<f64>>; 3],
    pub dphi_dxi: Vec<Vec<f64>>,
    pub dphi_deta: Vec<Vec<f64>>,
    pub mass_mat: Vec<f64>,
    pub gauss_points: &'a GaussPoints,
}
impl DubinerBasis<'_> {
    fn compute_phi_cell(&mut self, cell_points: Vec<[f64; 2]>, gp_number: usize) {
        for i in 0..gp_number {
            let xi = cell_points[i][0];
            let eta = cell_points[i][1];
            self.phi_cell_gp[i][0] = 1.0;
            self.phi_cell_gp[i][1] = xi;
            self.phi_cell_gp[i][2] = eta;
            self.phi_cell_gp[i][3] = xi * (2.0 * xi - 1.0);
            self.phi_cell_gp[i][4] = 4.0 * xi * eta;
            self.phi_cell_gp[i][5] = eta * (2.0 * eta - 1.0);
        }
    }
    fn compute_phi_edge(&mut self) {
        for i in 0..3 {
            for (phi_gp, point) in self.phi_edge_gp[i].iter_mut().zip(self.gauss_points.edge_points[i].iter()) {
                let xi = point[0];
                let eta = point[1];
                phi_gp[0] = 1.0;
                phi_gp[1] = xi;
                phi_gp[2] = eta;
                phi_gp[3] = xi * (2.0 * xi - 1.0);
                phi_gp[4] = 4.0 * xi * eta;
                phi_gp[5] = eta * (2.0 * eta - 1.0);
            }
        }
    }
    fn compute_derivatives(&mut self, gauss_points: Vec<[f64; 2]>, gp_number: usize) {
        for i in 0..gp_number {
            let xi = gauss_points[i][0];
            let eta = gauss_points[i][1];
            self.dphi_dxi[i][0] = 0.0;
            self.dphi_dxi[i][1] = 1.0;
            self.dphi_dxi[i][2] = 0.0;
            self.dphi_dxi[i][3] = 4.0 * xi - 1.0;
            self.dphi_dxi[i][4] = 4.0 * eta;
            self.dphi_dxi[i][5] = 0.0;
            self.dphi_deta[i][0] = 0.0;
            self.dphi_deta[i][1] = 0.0;
            self.dphi_deta[i][2] = 1.0;
            self.dphi_deta[i][3] = 0.0;
            self.dphi_deta[i][4] = 4.0 * xi;
            self.dphi_deta[i][5] = 4.0 * eta - 1.0;
        }
    }
}