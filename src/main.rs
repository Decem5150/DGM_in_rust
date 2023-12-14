mod spatial_disc;
mod temporal_disc;
mod solver;
mod mesh;
mod basis_function;
mod gauss_point;
mod io;
mod initialization;
fn main() {
    let mut solver = solver::Solver::default();
    let mut mesh = mesh::Mesh::default();
}
