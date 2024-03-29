use std::collections::HashMap;
use ndarray::{Array, Ix1};
use crate::basis_function::{DubinerBasis, GaussPoints};
use crate::io::boundary_condition_parser::FreeStreamCondition;
use crate::io::mesh_parser::{EdgeBlock, ElementBlock, GmshParser};
use crate::io::param_parser::Parameters;
use crate::mesh::{BoundaryEdge, BoundaryType, Vertex, Edge, EdgeTypeAndIndex, Element, Patch, Mesh, BoundaryQuantity};
use crate::solver::{SolverParameters, Solver, FlowParameters, MeshParameters};
use crate::spatial_disc::{InviscidFluxScheme, SpatialDisc};
use crate::temporal_disc::{TimeScheme, TemperalDisc};
pub fn initialize_mesh_basis_gauss_params() -> (Mesh, DubinerBasis, GaussPoints, FlowParameters, MeshParameters, SolverParameters) {
    let parameters = Parameters::parse();
    let mut solver_parameters = SolverParameters {
        cfl: 0.0,
        final_time: parameters.final_time,
        final_step: parameters.final_step,
        number_of_cell_gp: parameters.number_of_cell_gp,
        number_of_edge_gp: parameters.number_of_edge_gp,
        number_of_equations: parameters.number_of_equations,
        number_of_basis_functions: parameters.number_of_basis_functions,
        order_of_polynomials: parameters.order_of_polynomials,
        inviscid_flux_scheme: {
            match parameters.inviscid_flux_scheme.to_lowercase().as_str() {
                "hllc" => InviscidFluxScheme::HLLC,
                _ => panic!("Inviscid flux scheme not implemented"),
            }
        }
    };
    solver_parameters.cfl = {
        let p = solver_parameters.order_of_polynomials as f64;
        1.0 / ((2.0 * p + 1.0) * (1.0 + 4.0 / (p + 2.0).powf(2.0)))
    };
    dbg!(&solver_parameters.cfl);
    let flow_parameters = FlowParameters {
        hcr: parameters.hcr,
    };
    let gauss_points = GaussPoints::new(solver_parameters.number_of_cell_gp, solver_parameters.number_of_edge_gp);
    let basis = DubinerBasis::new(solver_parameters.number_of_basis_functions, &gauss_points);
    let (vertices_vec, boundary_edge_blocks, element_blocks) = GmshParser::parse(parameters.mesh_file.clone());
    let (vertices, edges, boundary_edges, elements, patches, indices_internal_elements, indices_boundary_elements) = process_mesh(&solver_parameters, vertices_vec, boundary_edge_blocks, element_blocks);
    let mesh_parameters = MeshParameters {
        number_of_elements: elements.len(),
        number_of_edges: edges.len(),
        number_of_vertices: vertices.len(),
        number_of_patches: patches.len(),
    };
    let mut mesh = Mesh {
        vertices,
        edges,
        boundary_edges,
        elements,
        patches,
        indices_internal_elements,
        indices_boundary_elements,
    };
    mesh.compute_jacob_det();
    mesh.compute_normal();
    //mesh.compute_mass_mat(&basis, &gauss_points, solver_parameters.number_of_cell_gp, solver_parameters.number_of_basis_functions);
    mesh.compute_derivatives(&basis, solver_parameters.number_of_cell_gp, solver_parameters.number_of_basis_functions);
    mesh.compute_circumradius();
    mesh.compute_minimum_height();
    mesh.set_neighbours();
    (mesh, basis, gauss_points, flow_parameters, mesh_parameters, solver_parameters)
}
pub fn initialize_solver<'a>(mesh: &'a Mesh, basis: &'a DubinerBasis, gauss_points: &'a GaussPoints, flow_parameters: &'a FlowParameters, mesh_parameters: &'a MeshParameters, solver_parameters: &'a SolverParameters) -> Solver<'a> {
    let time_scheme = TimeScheme::RK3;
    let residuals = Array::zeros((mesh_parameters.number_of_elements, solver_parameters.number_of_equations, solver_parameters.number_of_basis_functions));
    let solutions = Array::zeros((mesh_parameters.number_of_elements, solver_parameters.number_of_equations, solver_parameters.number_of_basis_functions));
    let spatial_disc = SpatialDisc {
        mesh: &mesh,
        basis: &basis,
        gauss_points: &gauss_points,
        flow_param: &flow_parameters,
        mesh_param: &mesh_parameters,
        solver_param: &solver_parameters,
    };
    let nelem = mesh_parameters.number_of_elements;
    let neq = solver_parameters.number_of_equations;
    let nbasis = solver_parameters.number_of_basis_functions;
    let temporary_solutions = vec![Array::zeros((nelem, neq, nbasis)); 2];
    let temporal_disc = TemperalDisc {
        temporary_solutions,
        current_time: 0.0,
        current_step: 0,
        time_scheme,
        mesh: &mesh,
        basis: &basis,
        gauss_points: &gauss_points,
        flow_param: &flow_parameters,
        mesh_param: &mesh_parameters,
        solver_param: &solver_parameters,
    };
    let solver = Solver {
        residuals,
        solutions,
        spatial_disc,
        temporal_disc,
        mesh: &mesh,
        basis: &basis,
        gauss_points: &gauss_points,
        flow_param: &flow_parameters,
        mesh_param: &mesh_parameters,
        solver_param: &solver_parameters,
    };
    solver
}
pub fn process_mesh(solver_param: &SolverParameters, vertices_vec: Vec<Vertex>, boundary_edge_blocks: Vec<EdgeBlock>, element_blocks: Vec<ElementBlock>) -> (Array<Vertex, Ix1>, Array<Edge, Ix1>, Array<BoundaryEdge, Ix1>, Array<Element, Ix1>, Array<Patch, Ix1>, Array<usize, Ix1>, Array<usize, Ix1>) {
    let mut boundary_edges_vec = Vec::new();
    let mut edges_vec = Vec::new();
    let mut elements_vec = Vec::new();
    let mut patches_vec = Vec::new();
    let mut edge_map: HashMap<[usize; 2], EdgeTypeAndIndex> = HashMap::new();
    let mut indices_internal_elements_vec = Vec::new();
    let mut indices_boundary_elements_vec = Vec::new();
    let free_stream_condition = FreeStreamCondition::parse();
    for block in boundary_edge_blocks.into_iter() {
        let physical_name = block.physical_name;
        let (boundary_type, boundary_quantity) = match physical_name.to_lowercase().as_str() {
            "wall" => (BoundaryType::Wall, None),
            "farfield" => {
                let sound_speed = (free_stream_condition.hcr * free_stream_condition.freestream_pressure / free_stream_condition.freestream_density).sqrt();
                (BoundaryType::FarField, Some(BoundaryQuantity {
                    rho: free_stream_condition.freestream_density,
                    u: free_stream_condition.mach_number * sound_speed * free_stream_condition.angle_of_attack.cos(),
                    v: free_stream_condition.mach_number * sound_speed * free_stream_condition.angle_of_attack.sin(),
                    p: free_stream_condition.freestream_pressure,
                }))
            },
            _ => panic!("Boundary type not implemented"),
        };
        let start_index = boundary_edges_vec.len();
        boundary_edges_vec.reserve(block.node_pairs.len());
        for node_pair in block.node_pairs.into_iter() {
            let ivertices = node_pair;
            let normal = [0.0, 0.0];
            let jacob_det = 0.0;
            let ielement = 0;
            let in_cell_index = 0;
            let (node_1, node_2) = (ivertices[0], ivertices[1]);
            boundary_edges_vec.push(BoundaryEdge {
                ivertices,
                normal,
                jacob_det,
                ielement,
                in_cell_index,
                boundary_type: boundary_type.clone(),
            });
            edge_map.insert([node_1, node_2], EdgeTypeAndIndex::Boundary(boundary_edges_vec.len() - 1));
        }
        let end_index = boundary_edges_vec.len() - 1;
        let patch = Patch {
            iedges: (start_index..=end_index).collect(),
            boundary_type,
            boundary_quantity: boundary_quantity,
        };
        patches_vec.push(patch);
    }
    for block in element_blocks.into_iter() {
        edges_vec.reserve(block.node_ids.len() * 3);
        elements_vec.reserve(block.node_ids.len());
        // create edges from connectivities of elements
        for &node_id in block.node_ids.iter() {
            let mut node_pairs = [[0, 0]; 3];
            for i in 0..3 {
                node_pairs[i] = [node_id[i], node_id[(i + 1) % 3]];
            }
            for node_pair in node_pairs.iter() {
                if edge_map.get(node_pair).is_none() && edge_map.get(&[node_pair[1], node_pair[0]]).is_none() {
                    let ivertices = *node_pair;
                    let normal = [0.0, 0.0];
                    let jacob_det = 0.0;
                    let ielements = [0, 0];
                    let in_cell_indices = [0, 0];
                    edges_vec.push(Edge {
                        ivertices,
                        normal,
                        jacob_det,
                        ielements,
                        in_cell_indices,
                    });
                    edge_map.insert(*node_pair, EdgeTypeAndIndex::Internal(edges_vec.len() - 1));
                }
            }
        }
        // assign edges to elements
        for node_id in block.node_ids.iter() {
            let mut node_pairs = [[0, 0]; 3];
            for i in 0..3 {
                node_pairs[i] = [node_id[i], node_id[(i + 1) % 3]];
            }
            let mut ivertices = Vec::new();
            let mut iedges = Vec::new();
            for node_pair in node_pairs.iter() {
                ivertices.push(node_pair[0]);
                let edge_type_and_index = edge_map.get(node_pair).or_else(|| {
                    let mut swapped_node_pair = node_pair.clone();  // Swap only if the first get fails
                    swapped_node_pair.swap(0, 1);
                    edge_map.get(&swapped_node_pair)
                }).unwrap();
                iedges.push(edge_type_and_index.clone());
            }
            // store the indices of the elements in the edges
            for (in_cell_index, iedge_) in iedges.iter().enumerate() {
                match iedge_ {
                    EdgeTypeAndIndex::Boundary(iedge) => {
                        boundary_edges_vec[*iedge].ielement = elements_vec.len();
                        boundary_edges_vec[*iedge].in_cell_index = in_cell_index;
                    },
                    EdgeTypeAndIndex::Internal(iedge) => {
                        if ivertices[in_cell_index] == edges_vec[*iedge].ivertices[0] {
                            edges_vec[*iedge].ielements[0] = elements_vec.len();
                            edges_vec[*iedge].in_cell_indices[0] = in_cell_index;
                        } else {
                            edges_vec[*iedge].ielements[1] = elements_vec.len();
                            edges_vec[*iedge].in_cell_indices[1] = in_cell_index;
                        }
                    }
                }
            }
            let is_internal_element = iedges.iter().all(|iedge| matches!(iedge, EdgeTypeAndIndex::Internal(_)));
            let ivertices = Array::from(ivertices);
            let iedges = Array::from(iedges);
            let ineighbours = Array::from(vec![None; 3]);
            let ngp = solver_param.number_of_cell_gp;
            let nbasis = solver_param.number_of_basis_functions;
            let derivatives: Array<HashMap<(usize, usize), f64>, _> = Array::from_shape_fn((ngp, nbasis), |_| {
                HashMap::new()
            });
            let jacob_det = 0.0;
            let circumradius = 0.0;
            let minimum_height = 0.0;
            elements_vec.push(Element {
                ivertices,
                iedges,
                ineighbours,
                derivatives,
                jacob_det,
                circumradius,
                minimum_height,
            });
            if is_internal_element {
                indices_internal_elements_vec.push(elements_vec.len() - 1);
            } else {
                indices_boundary_elements_vec.push(elements_vec.len() - 1);
            }
        }
    }
    let vertices = Array::from(vertices_vec);
    let edges = Array::from(edges_vec);
    let boundary_edges = Array::from(boundary_edges_vec);
    let elements = Array::from(elements_vec);
    let patches = Array::from(patches_vec);
    let indices_internal_elements = Array::from(indices_internal_elements_vec);
    let indices_boundary_elements = Array::from(indices_boundary_elements_vec);
    (vertices, edges, boundary_edges, elements, patches, indices_internal_elements, indices_boundary_elements)
}