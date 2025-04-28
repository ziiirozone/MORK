#[derive(Debug, Clone, PartialEq)]
enum NodeState {
    NotVisited,
    InPath(usize),
    InSCC(usize),
}

pub fn scc(graph: &Vec<Vec<bool>>) -> Vec<Vec<usize>> {
    #[allow(non_snake_case)]
    let mut S: Vec<usize> = Vec::new();
    #[allow(non_snake_case)]
    let mut B: Vec<usize> = Vec::new();
    let s = graph.len();
    if s == 0 {
        return Vec::new();
    }
    let mut scc_count: usize = 0;
    let mut current_node: usize;
    let mut state: Vec<NodeState> = vec![NodeState::NotVisited; s];
    let mut tree: Vec<Vec<usize>> = vec![(0..s).collect()];
    let mut l: usize;
    let mut path: Vec<usize> = Vec::new();
    loop {
        l = tree.len() - 1;
        match tree[l].pop() {
            Some(j) => {
                // next node to search
                match state[j] {
                    NodeState::NotVisited => {
                        current_node = j;
                        path.push(current_node);
                        S.push(current_node);
                        state[current_node] = NodeState::InPath(S.len() - 1);
                        B.push(S.len() - 1);
                        tree.push(
                            (0..s)
                                .filter(|&j| graph[current_node][j] && current_node != j)
                                .collect(),
                        );
                    }
                    NodeState::InPath(j1) => {
                        while j1 < B[B.len() - 1] {
                            B.pop();
                        }
                    }
                    _ => {}
                }
            }
            None => {
                // all edges of current node were searched
                tree.pop();
                match path.pop() {
                    Some(j) => {
                        // assign to the scc
                        current_node = j;
                        if let NodeState::InPath(j1) = state[current_node] {
                            if j1 == B[B.len() - 1] {
                                B.pop();
                                while j1 + 1 <= S.len() {
                                    state[S.pop().unwrap()] = NodeState::InSCC(scc_count);
                                }
                                scc_count += 1;
                            }
                        }
                    }
                    None => {
                        break;
                    } // no more path to search
                }
            }
        }
    }
    let mut scc = vec![Vec::new(); scc_count];
    for (j, st) in state.into_iter().enumerate() {
        if let NodeState::InSCC(c) = st {
            scc[c].push(j);
        }
    }
    scc
}

pub fn contraction(graph: &Vec<Vec<bool>>, partition: &Vec<Vec<usize>>) -> Vec<Vec<bool>> {
    let s1 = graph.len();
    let s2 = partition.len();
    let inverse: Vec<usize> = (0..s1)
        .map(|j| {
            (0..s2)
                .filter(|j1| partition[*j1].contains(&j))
                .next()
                .unwrap()
        })
        .collect(); // maps each node to its partition
    let mut contracted_graph = vec![vec![false; s2]; s2];
    let mut p1;
    let mut p2;
    for j1 in 0..s1 {
        for j2 in 0..s1 {
            p1 = inverse[j1];
            p2 = inverse[j2];
            if p1 != p2 && graph[j1][j2] {
                contracted_graph[p1][p2] = true;
            }
        }
    }
    contracted_graph
}

pub fn topological_sort(dag: &Vec<Vec<bool>>) -> Vec<usize> {
    let mut dag = dag.clone();
    let s = dag.len();
    #[allow(non_snake_case)]
    let mut S: Vec<usize> = (0..s).filter(|&j| !(dag[j].contains(&true))).collect();
    let mut order: Vec<usize> = Vec::new();
    while let Some(node) = S.pop() {
        order.push(node);
        for k in 0..s {
            dag[k][node] = false;
        }
        S = (0..dag.len())
            .filter(|&j| !order.contains(&j) && !(dag[j].contains(&true)))
            .collect()
    }
    order
}
