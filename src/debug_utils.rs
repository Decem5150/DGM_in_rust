use std::fs::File;

use csv::Writer;
use ndarray::{Array, ArrayBase, Data, Ix1};

pub fn check_for_nan<S,D>(name: String, array: &ArrayBase<S, D>)
where
    S: Data<Elem = f64>, 
    D: ndarray::Dimension,
{
    for &value in array.iter() {
        if value.is_nan() {
            panic!("Found NaN value in {}!", name);
        }
    }
}
pub fn write_to_csv(solutions: Array<f64, Ix1>, current_step: usize) {
    let file_name = format!("outputfiles/debug/step={}.csv", current_step);
    let file = File::create(&file_name).unwrap();
    let mut writer = Writer::from_writer(file);
    writer.write_record(&["point 1", "point 2", "point 3", "point 4", "point5"]).unwrap();
    writer.write_record(&[solutions[0].to_string(), solutions[1].to_string(), solutions[2].to_string(), solutions[3].to_string(), solutions[4].to_string()]).unwrap();
    writer.flush().unwrap();
}