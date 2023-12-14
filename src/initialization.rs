use ndarray::{Array, Ix1};
use crate::io::mesh_parser::{EdgeBlock, ElementBlock};
use crate::mesh::{BoundaryEdge, Vertex, Edge, Element, Patch};
struct EdgeBuilder {
    pub ielements: Option<[usize; 2]>,
    pub ivertices: Option<[usize; 2]>,
}
struct ElementBuilder {
    pub ivertices: Array<usize, Ix1>,
    pub iedges: Array<usize, Ix1>,
    pub ineighbours: Array<usize, Ix1>,
}
pub struct MeshBuilder {

}
impl MeshBuilder {
    pub fn initialize_mesh(vertices_vec: Vec<Vertex>, boundary_edge_blocks: Vec<EdgeBlock>, element_blocks: Vec<ElementBlock>) -> (Array<Vertex, Ix1>, Array<Edge, Ix1>, Array<BoundaryEdge, Ix1>, Array<Element, Ix1>, Array<Patch, Ix1>) {
        let vertices = Array::from(vertices_vec);
        let mut boundary_edges_vec = Vec::new();
        let mut elements_vec = Vec::new();
        for block in boundary_edge_blocks.into_iter() {
            let physical_name = block.physical_name;
            boundary_edges_vec.extend(
        }
    }
}