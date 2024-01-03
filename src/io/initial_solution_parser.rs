use std::fs;
use serde::Deserialize;
use serde_json;
#[derive(Deserialize)]
pub struct InitialSolution {
    pub density: f64,
    pub u: f64,
    pub v: f64,
    pub pressure: f64,
}
impl InitialSolution {
    pub fn parse() -> Self {
        let data = fs::read_to_string("inputfiles/initial_solution.json").expect("Failed to read initial_solution.json");
        let initial_solution: Self = serde_json::from_str(&data).expect("Failed to parse initial solution");
        initial_solution
    }
}