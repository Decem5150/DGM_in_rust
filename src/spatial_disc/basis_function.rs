pub struct DubinerBasis {
    pub dof: usize,
    pub phi_cell_gp: Vec<Vec<f64>>,
    pub phi_edge_gp: Vec<Vec<f64>>,
    pub dphi_dxi: Vec<Vec<f64>>,
    pub dphi_deta: Vec<Vec<f64>>,
    pub mass_mat: Vec<f64>,
}
impl DubinerBasis {
    fn compute_phi_gauss(&mut self, gauss_points: Vec<[f64; 2]>, gp_number: usize) {
        for i in 0..gp_number {
            let xi = gauss_points[i][0];
            let eta = gauss_points[i][1];
            self.phi_cell_gp[i][0] = 1.0;
            self.phi_cell_gp[i][1] = xi;
            self.phi_cell_gp[i][2] = eta;
            self.phi_cell_gp[i][3] = xi * (2.0 * xi - 1.0);
            self.phi_cell_gp[i][4] = 4.0 * xi * eta;
            self.phi_cell_gp[i][5] = eta * (2.0 * eta - 1.0);
        }
    }
    fn compute_phi_edge(&mut self, edge_gauss_points: Vec<f64>, gp_number: usize) {
        for i in 0..gp_number {
            let x = edge_gauss_points[i];
            self.phi_edge_gp[i][0] = 1.0;
            self.phi_edge_gp[i][1] = x;
            self.phi_edge_gp[i][2] = 0.5 * (3.0 * x.powi(2) - 1.0);
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
    fn compute_mass_mat(&mut self) {
        
    }
}