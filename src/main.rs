use io::write_results::write_to_vtk;

mod spatial_disc;
mod temporal_disc;
mod solver;
mod mesh;
mod basis_function;
mod gauss_point;
mod io;
mod initialization;
mod debug_utils;
fn main() {
    let (mesh, basis_function, gauss_points, flow_param, mesh_param, solver_param) = initialization::initialize_mesh_basis_gauss_params();
    let mut solver = initialization::initialize_solver(&mesh, &basis_function, &gauss_points, &flow_param, &mesh_param, &solver_param);
    //dbg!(&basis_function.phis_cell_gps);
    solver.set_initial_solution();
    /* 
    dbg!(&mesh.elements[2759]);
    dbg!(&mesh.boundary_edges[7]);
    dbg!(&mesh.edges[4742]);
    dbg!(&mesh.edges[4743]);
    */
    solver.solve();
    write_to_vtk(&solver.solutions, &mesh, &basis_function, &gauss_points, &flow_param, &mesh_param, &solver_param);
}
