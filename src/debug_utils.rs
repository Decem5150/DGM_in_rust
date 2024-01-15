use ndarray::{ArrayBase, Data};

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