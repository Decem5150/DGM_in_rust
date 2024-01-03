use std::fs;
use serde::Deserialize;
use serde_json;
#[derive(Deserialize)]
pub struct FreeStreamCondition {
    pub mach_number: f64,
    pub angle_of_attack: f64,
    pub density: f64,
    pub pressure: f64,
    pub hcr: f64,
}
impl FreeStreamCondition {
    pub fn parse() -> Self {
        let data = fs::read_to_string("inputfiles/freestream_condition.json").expect("Failed to read freestream_condition.json");
        let freestream_condition: Self = serde_json::from_str(&data).expect("Failed to parse freestream condition");
        freestream_condition
    }
}