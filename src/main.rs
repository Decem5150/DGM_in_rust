mod spatial_disc;
mod temporal_disc;
mod solver;
mod mesh;
mod basis_function;
mod gauss_point;
mod io;
mod initialization;
fn main() {
    let (mut mesh, mut basis_function, mut gauss_points, flow_param, mesh_param, solver_param) = initialization::Initializer::initialize_mesh_basis_gauss_params();
    let mut solver = initialization::Initializer::initialize_solver(&mesh, &basis_function, &gauss_points, &flow_param, &mesh_param, &solver_param);
    solver.set_initial_solution();

}
