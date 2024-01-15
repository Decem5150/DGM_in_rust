use std::fs;
use serde::Deserialize;
use serde_json;
#[derive(Deserialize)]
pub struct FreeStreamCondition {
    pub angle_of_attack: f64,
    pub mach_number: f64,
    pub freestream_density: f64,
    pub freestream_pressure: f64,
    pub freestream_temperature: f64,
    pub hcr: f64,
}
impl FreeStreamCondition {
    pub fn parse() -> Self {
        let data = fs::read_to_string("inputfiles/freestream.json").expect("Failed to read freestream_condition.json");
        let freestream_condition: Self = serde_json::from_str(&data).expect("Failed to parse freestream condition");
        freestream_condition
    }
}