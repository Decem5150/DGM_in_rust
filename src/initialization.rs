use std::collections::HashMap;

use ndarray::{Array, Ix1};
use crate::io::mesh_parser::{EdgeBlock, ElementBlock};
use crate::mesh::{BoundaryEdge, BoundaryType, Vertex, Edge, Element, Patch};
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
        let mut edges_vec = Vec::new();
        let mut elements_vec = Vec::new();
        let mut patches_vec = Vec::new();
        let mut boundary_edge_map = HashMap::new();
        let mut edge_map = HashMap::new();
        for block in boundary_edge_blocks.into_iter() {
            let physical_name = block.physical_name;
            let start_index = boundary_edges_vec.len();
            boundary_edges_vec.reserve(block.node_pairs.len());
            for node_pair in block.node_pairs.into_iter() {
                let ivertices = node_pair;
                let normal = [0.0, 0.0];
                let jacob_det = 0.0;
                let ielement = 0;
                let in_cell_index = 0;
                let (node_1, node_2) = {
                    if ivertices[0] < ivertices[1] {
                        (ivertices[0], ivertices[1])
                    } else {
                        (ivertices[1], ivertices[0])
                    }
                };
                boundary_edges_vec.push(BoundaryEdge {
                    ivertices,
                    normal,
                    jacob_det,
                    ielement,
                    in_cell_index,
                });
                boundary_edge_map.insert((node_1, node_2), boundary_edges_vec.len() - 1);
            }
            let end_index = boundary_edges_vec.len() - 1;
            let boundary_type = match physical_name.as_str() {
                "wall" => BoundaryType::Wall,
                "farfield" => BoundaryType::FarField,
                _ => panic!("Boundary type not implemented"),
            };
            let patch = Patch {
                iedges: (start_index..end_index).collect(),
                boundary_type,
                boundary_quantity: None,
            };
        }
        for block in element_blocks.into_iter() {
            let start_index = elements_vec.len();
            edges_vec.reserve(block.node_ids.len());
            for &node_id in block.node_ids.iter() {
                let mut node_pairs = vec![]
                let node_pair_0 = {
                    if node_id[0] < node_id[1] {
                        (node_id[0], node_id[1])
                    } else {
                        (node_id[1], node_id[0])
                    }
                };
                let node_pair_1 = {
                    if node_id[1] < node_id[2] {
                        (node_id[1], node_id[2])
                    } else {
                        (node_id[2], node_id[1])
                    }
                };
                let node_pair_2 = {
                    if node_id[2] < node_id[0] {
                        (node_id[2], node_id[0])
                    } else {
                        (node_id[0], node_id[2])
                    }
                };
                if edge_map.get(node_pair_0
            }
            let end_index = elements_vec.len();
            let patch = Patch {
                iedges: (start_index..end_index).collect(),
                boundary_type: BoundaryType::Internal,
                boundary_quantity: None,
            };
            patches_vec.push(patch);
        }
    }
}