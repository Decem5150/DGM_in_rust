use io::write_results::write_to_vtk;

mod spatial_disc;
mod temporal_disc;
mod solver;
mod mesh;
mod basis_function;
mod gauss_point;
mod io;
mod initialization;
fn main() {
    let (mesh, basis_function, gauss_points, flow_param, mesh_param, solver_param) = initialization::initialize_mesh_basis_gauss_params();
    let mut solver = initialization::initialize_solver(&mesh, &basis_function, &gauss_points, &flow_param, &mesh_param, &solver_param);
    //dbg!(&basis_function);
    solver.set_initial_solution();
    //solver.solve();
    write_to_vtk(&solver.solutions, &mesh, &basis_function, &gauss_points, &flow_param, &mesh_param, &solver_param);
}
