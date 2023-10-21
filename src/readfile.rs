use std::fs::File;
use std::io::{self, BufRead};
use std::path::Path;
fn read_gmsh_mesh(filename: &str) -> io::Result<(Vec<Vertex>, Vec<Element>, Vec<Edge>)> {
    let path = Path::new(filename);
    let file = File::open(&path)?;
    let reader = io::BufReader::new(file);
    let mut lines = reader.lines();

    let mut vertices: Vec<Vertex> = Vec::new();
    let mut elements: Vec<Element> = Vec::new();
    let mut edges: Vec<Edge> = Vec::new();

    while let Some(Ok(line)) = lines.next() {
        if line == "$Nodes" {
            break;
        }
    }
    // Read the meta-information line
    let meta_info_line = lines.next().unwrap()?;
    let parts: Vec<&str> = meta_info_line.split_whitespace().collect();
    let total_blocks: usize = parts[0].parse().unwrap();
    let total_nodes: usize = parts[1].parse().unwrap();
    vertices.resize(total_nodes, Vertex::default());
    // Read the nodes
    for i in 0..total_blocks {
        let block_info_line = lines.next().unwrap()?;
        let parts: Vec<&str> = block_info_line.split_whitespace().collect();
        let num_of_nodes: usize = parts[3].parse().unwrap();
        for _ in 0..num_of_nodes {
            lines.next();
        }
        for _ in 0..num_of_nodes {
            let node_line = lines.next().unwrap()?;
            let parts: Vec<&str> = node_line.split_whitespace().collect();
            let node_id: usize = parts[0].parse().unwrap();
            let x: f64 = parts[1].parse().unwrap();
            let y: f64 = parts[2].parse().unwrap();
            vertices[node_id - 1] = Vertex { x, y };
        }
    }
    while let Some(Ok(line)) = lines.next() {
        if line == "$Elements" {
            break;
        }
    }
    // Read the meta-information line
    let meta_info_line = lines.next().unwrap()?;
    let parts: Vec<&str> = meta_info_line.split_whitespace().collect();
    let total_blocks: usize = parts[0].parse().unwrap();
    let total_elements: usize = parts[1].parse().unwrap();
    elements.resize(total_elements, Element::default());
    // Read the elements
    for i in 0..total_blocks {
        let block_info_line = lines.next().unwrap()?;
        let parts: Vec<&str> = block_info_line.split_whitespace().collect();
        let entity_dim: usize = parts[0].parse().unwrap();
        let entity_tag: usize = parts[1].parse().unwrap();
        let element_type: usize = parts[2].parse().unwrap();
        let num_of_elements: usize = parts[3].parse().unwrap();
        if element_type == 1 {
            
        }
        if element_type == 2 {
            for _ in 0..num_of_elements {
                let element_line = lines.next().unwrap()?;
                let parts_in_block: Vec<&str> = element_line.split_whitespace().collect();
                let element_id: usize = parts_in_block[0].parse().unwrap();
                let mut node_tags: Vec<usize> = vec![0;3];
                for j in 0..num_of_tags {
                    tags[j] = parts_in_block[j + 1].parse().unwrap();
                }
                let mut element = Element{vertices: node_tags, edges: vec![0;3], neighbours: vec![0;3], area: 0.0};
                elements[element_id - 1] = element;
            }
        }
        

    }
    for line in reader.lines() {
        let line = line?;
        /*
        if line.starts_with("$Nodes") {
            in_nodes_section = true;
            continue;
        } else if line.starts_with("$EndNodes") {
            in_nodes_section = false;
            continue;
        } else if line.starts_with("$Elements") {
            in_elements_section = true;
            continue;
        } else if line.starts_with("$EndElements") {
            in_elements_section = false;
            continue;
        }
        */
}