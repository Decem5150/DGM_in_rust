use crate::mesh;
pub fn runge_kutta_4<'a>(mesh: &'a mesh::Mesh, time_step: f64, time: f64, u: &'a Vec<f64>, f: &'a dyn Fn(&'a Vec<f64>, &'a Vec<f64>, f64) -> Vec<f64>) -> Vec<f64> {
    let mut k1 = f(u, &time, time_step);
    let mut k2 = f(&mesh.add_vec(u, &mesh.scalar_mult(&k1, 0.5)), &time, time_step);
    let mut k3 = f(&mesh.add_vec(u, &mesh.scalar_mult(&k2, 0.5)), &time, time_step);
    let mut k4 = f(&mesh.add_vec(u, &k3), &time, time_step);
    let mut u_new = mesh.add_vec(u, &mesh.scalar_mult(&k1, 1.0/6.0));
    u_new = mesh.add_vec(&u_new, &mesh.scalar_mult(&k2, 1.0/3.0));
    u_new = mesh.add_vec(&u_new, &mesh.scalar_mult(&k3, 1.0/3.0));
    u_new = mesh.add_vec(&u_new, &mesh.scalar_mult(&k4, 1.0/6.0));
    u_new
}