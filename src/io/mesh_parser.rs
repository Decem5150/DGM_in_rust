use std::collections::HashMap;
use std::{path::Path, fs::File, io};
use std::io::BufRead;
use crate::mesh::Vertex;
pub struct EdgeBlock {
    pub node_pairs: Vec<[usize; 2]>,
    pub physical_name: String,
}
pub struct ElementBlock {
    pub node_ids: Vec<[usize; 3]>,
}
enum CurrentState {
    StandingBy,
    ReadingPhysicalNames,
    ReadingFirstLineOfPhysicalNames,
    ReadingEntities,
    ReadingFirstLineOfEntities,
    ReadingNodes,
    ReadingFirstLineOfNodes,
    ReadingFirstLineOfElements,
    ReadingEdgeBlocks,
    ReadingElementBlocks,
    ReadingFirstLineOfBlocks,
}
pub struct GmshParser;
impl GmshParser {
    pub fn parse(filename: String) -> (Vec<Vertex>, Vec<EdgeBlock>, Vec<ElementBlock>) {
        let path = Path::new(&filename);
        let file = File::open(path).expect("Failed to open file");
        let reader = io::BufReader::new(file);
        let mut physical_names: HashMap<usize, String> = HashMap::new();
        let mut entities_to_physical_tags: HashMap<usize, usize> = HashMap::new();
        let mut nodes: Vec<Vertex> = Vec::new();
        let mut edges: Vec<EdgeBlock> = Vec::new();
        let mut elements: Vec<ElementBlock> = Vec::new();
        let mut current_state = CurrentState::StandingBy;
        let mut line_count: usize = 0;
        let mut max_line_count: usize = 0;
        for line in reader.lines() {
            let line = line.expect("Failed to read line");
            // Check the current section of the file
            if line.starts_with("$End") {
                current_state = CurrentState::StandingBy;
                continue;
            } else if line.starts_with("$PhysicalNames") {
                current_state = CurrentState::ReadingFirstLineOfPhysicalNames;
                continue;
            } else if line.starts_with("$Entities") {
                current_state = CurrentState::ReadingFirstLineOfEntities;
                continue;
            } else if line.starts_with("$Nodes") {
                current_state = CurrentState::ReadingFirstLineOfNodes;
                continue;
            } else if line.starts_with("$Elements") {
                current_state = CurrentState::ReadingFirstLineOfElements;
                continue;
            } 
            match current_state {
                CurrentState::StandingBy => continue,
                CurrentState::ReadingFirstLineOfPhysicalNames => {
                    current_state = CurrentState::ReadingPhysicalNames;
                    continue;
                },
                CurrentState::ReadingPhysicalNames => {
                    let values: Vec<&str> = line.split_whitespace().collect();
                    if values.len() == 3 {
                        let physical_tag = values[1].parse::<usize>().unwrap();
                        let physical_name = values[2].trim_matches('"').to_string();
                        physical_names.insert(physical_tag, physical_name);
                    }
                    continue;
                },
                CurrentState::ReadingFirstLineOfEntities => {
                    current_state = CurrentState::ReadingEntities;
                    continue;
                },
                CurrentState::ReadingEntities => {
                    let values: Vec<&str> = line.split_whitespace().collect();
                    if values.len() >= 11 {
                        let entity_tag = values[0].parse::<usize>().unwrap();
                        let physical_tag = values[8].parse::<usize>().unwrap();
                        entities_to_physical_tags.insert(entity_tag, physical_tag);
                    }
                    continue;
                },
                CurrentState::ReadingFirstLineOfNodes => {
                    let values: Vec<&str> = line.split_whitespace().collect();
                    let num_of_nodes = values[1].parse::<usize>().unwrap();
                    nodes.reserve(num_of_nodes);
                    current_state = CurrentState::ReadingNodes;
                    continue;
                },
                CurrentState::ReadingNodes => {
                    let values: Vec<&str> = line.split_whitespace().collect();
                    if values.len() == 3 {
                        let x = values[0].parse::<f64>().unwrap();
                        let y = values[1].parse::<f64>().unwrap();
                        nodes.push(Vertex { x, y });
                    }
                    continue;
                },
                CurrentState::ReadingFirstLineOfElements => {
                    current_state = CurrentState::ReadingFirstLineOfBlocks;
                    continue;
                },
                CurrentState::ReadingFirstLineOfBlocks => {
                    let values: Vec<&str> = line.split_whitespace().collect();
                    let entity_dim = values[0].parse::<usize>().unwrap();
                    let entity_tag = values[1].parse::<usize>().unwrap();
                    let num_of_elements = values[3].parse::<usize>().unwrap();
                    let &physical_tag = entities_to_physical_tags.get(&entity_tag).unwrap();
                    let physical_name = physical_names.get(&physical_tag).unwrap().clone();
                    if entity_dim == 1 {
                        let mut edge_block = EdgeBlock{physical_name, node_pairs: Vec::new()};
                        edge_block.node_pairs.reserve(num_of_elements);
                        edges.push(edge_block);
                        current_state = CurrentState::ReadingEdgeBlocks;
                    }
                    if entity_dim == 2 {
                        let mut element_block = ElementBlock{node_ids: Vec::new()};
                        element_block.node_ids.reserve(num_of_elements);
                        elements.push(element_block);
                        current_state = CurrentState::ReadingElementBlocks;
                    }
                    line_count = 0;
                    max_line_count = num_of_elements;
                    continue;
                },
                CurrentState::ReadingEdgeBlocks => {
                    line_count += 1;
                    let values: Vec<&str> = line.split_whitespace().collect();
                    let edge_block = edges.last_mut().unwrap();
                    let node_id_1 = values[1].parse::<usize>().unwrap();
                    let node_id_2 = values[2].parse::<usize>().unwrap();
                    edge_block.node_pairs.push([node_id_1 - 1, node_id_2 - 1]);
                    if line_count == max_line_count {
                        current_state = CurrentState::ReadingFirstLineOfBlocks;
                    }
                    continue;
                },
                CurrentState::ReadingElementBlocks => {
                    line_count += 1;
                    let values: Vec<&str> = line.split_whitespace().collect();
                    let element_block = elements.last_mut().unwrap();
                    let node_id_1 = values[1].parse::<usize>().unwrap();
                    let node_id_2 = values[2].parse::<usize>().unwrap();
                    let node_id_3 = values[3].parse::<usize>().unwrap();
                    element_block.node_ids.push([node_id_1 - 1, node_id_2 - 1, node_id_3 - 1]);
                    if line_count == max_line_count {
                        current_state = CurrentState::ReadingFirstLineOfBlocks;
                    }
                    continue;
                },
            }
        }
        (nodes, edges, elements)
    }
}