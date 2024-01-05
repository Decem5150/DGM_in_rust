use std::fs;
use serde::Deserialize;
use serde_json;
#[derive(Deserialize, Debug)]
pub struct Parameters {
    pub mesh_file: String,
    pub cfl: f64,
    pub final_time: f64,
    pub number_of_cell_gp: usize,
    pub number_of_edge_gp: usize,
    pub number_of_equations: usize,
    pub number_of_basis_functions: usize,
    pub hcr: f64,
    pub inviscid_flux_scheme: String,
    pub time_scheme: String,
}
impl Parameters {
    pub fn parse() -> Self {
        let data = fs::read_to_string("inputfiles/parameters.json").expect("Failed to read parameters.json");
        let parameters: Self = serde_json::from_str(&data).expect("Failed to parse parameters");
        parameters
    }
}
