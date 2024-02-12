use ndarray::{Array, Ix3, s};
use vtkio::{model::{VertexNumbers, Cells, CellType, DataSet, UnstructuredGridPiece, Attributes, Version, ByteOrder, DataArray, Attribute, ElementType}, IOBuffer, Vtk};

use crate::{mesh::Mesh, basis_function::{DubinerBasis, GaussPoints}, solver::{FlowParameters, MeshParameters, SolverParameters}};

pub fn write_to_vtk(solutions: &Array<f64, Ix3>, mesh: &Mesh, basis_function: &DubinerBasis, gauss_points: &GaussPoints, flow_param: &FlowParameters, mesh_param: &MeshParameters, solver_param: &SolverParameters) {
    let points = IOBuffer::F64(mesh.vertices.iter().flat_map(|v| vec![v.x, v.y, 0.0]).collect());
    let cell_verts = VertexNumbers::XML {
        connectivity: mesh.elements.iter().flat_map(|elem| elem.ivertices.iter().map(|&index| index as u64)).collect(),
        offsets: {
            let mut offsets: Vec<u64> = Vec::new();
            let mut current_offset: u64 = 0;
            for elem in mesh.elements.iter() {
                // Get the number of vertices for the element
                let vertex_count = elem.ivertices.len() as u64;
                current_offset += vertex_count;
                offsets.push(current_offset);
            }
            offsets
        }
    };
    let cells = Cells {
        cell_verts,
        types: vec![CellType::Triangle; mesh.elements.len()]
    };
    let attributes = {
        let density = {
            let slice = solutions.slice(s![.., 0, 0]);
            slice.to_vec()
        };
        let x_momentum = {
            let slice = solutions.slice(s![.., 1, 0]);
            slice.to_vec()
        };
        let y_momentum = {
            let slice = solutions.slice(s![.., 2, 0]);
            slice.to_vec()
        };
        let energy = {
            let slice = solutions.slice(s![.., 3, 0]);
            slice.to_vec()
        };
        let x_velocity = x_momentum.iter().zip(&density).map(|(&x_momentum, &density)| x_momentum / density).collect::<Vec<f64>>();
        let y_velocity = y_momentum.iter().zip(&density).map(|(&y_momentum, &density)| y_momentum / density).collect::<Vec<f64>>();
        let velocity = x_velocity.iter().zip(&y_velocity).flat_map(|(&x_velocity, &y_velocity)| vec![x_velocity, y_velocity, 0.0]).collect::<Vec<f64>>();
        let pressure = energy.iter().zip(&density)
            .zip(&x_velocity)
            .zip(&y_velocity)
            .map(|(((&energy, &density), &x_vel), &y_vel)| (flow_param.hcr - 1.0) * (energy - 0.5 * density * (x_vel.powi(2) + y_vel.powi(2)))).collect::<Vec<f64>>();
        let density_attribute = Attribute::DataArray(DataArray {
            name: "Density".to_string(),
            elem: ElementType::Scalars { num_comp: 1, lookup_table: None },
            data: IOBuffer::F64(density),
        });
        let velocity_attribute = Attribute::DataArray(DataArray {
            name: "Velocity".to_string(),
            elem: ElementType::Vectors,
            data: IOBuffer::F64(velocity),
        });
        let pressure_attribute = Attribute::DataArray(DataArray {
            name: "Pressure".to_string(),
            elem: ElementType::Scalars { num_comp: 1, lookup_table: None },
            data: IOBuffer::F64(pressure),
        });
        Attributes {
            point: vec![],
            cell: vec![density_attribute, velocity_attribute, pressure_attribute],
        }
    };
    let dataset = DataSet::inline(UnstructuredGridPiece {
        points,
        cells,
        data: attributes,
    });
    let vtk = Vtk {
        version: Version::XML { major: 4, minor: 2 },
        title: String::new(),
        byte_order: ByteOrder::BigEndian,
        file_path: None,
        data: dataset,
    };
    vtk.export_ascii("outputfiles/airfoil solution at step=10.vtk")
        .expect(&format!("Failed to save file: {:?}", "outputfiles/solution.vtk"));
}