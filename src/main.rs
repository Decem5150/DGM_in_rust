mod spatial_disc;
mod temporal_disc;
mod solver;
mod mesh;
mod basis_function;
mod gauss_point;
mod io;
mod initialization;
fn main() {
    let mut solver = initialization::Initializer::initialize_solver();
    solver.set_initial_solution();

}
