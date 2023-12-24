use std::fs;
use serde::Deserialize;
use serde_json;
#[derive(Deserialize)]
pub struct Parameters {
    pub cfl: f64,
    pub final_time: f64,
    pub number_of_cell_gp: usize,
    pub number_of_edge_gp: usize,
    pub number_of_equations: usize,
    pub number_of_basis_functions: usize,
    pub hcr: f64,
}
pub struct ParametersParser;
impl ParametersParser {
    pub fn parse(&self) -> Parameters {
        let data = fs::read_to_string("inputfiles/parameters.json").expect("Failed to read parameters.json");
        let parameters: Parameters = serde_json::from_str(&data).expect("Failed to parse parameters");
        parameters
    }
}