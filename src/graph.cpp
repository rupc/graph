#include "graph.h"
// if adj list is initialized,
// and want to represent graph as matrix 
// for algorithm like a floyd warshall(matrix more suitable)
void AdjGraph::init_matrix() {
    if (Mats.size() > 0) {
        return;
    }
    int V = Nodes.size();
    Mats.resize(V);
    for (auto &p : Mats) {
        p.resize(V, DIST_INF);

    }
    for (auto &e : Edges) {
        Mats[EDGE_SRC(e)-1][EDGE_DST(e)-1] = 
            EDGE_WEIGHT(e);
    }
}
void init_romania(std::istream &pin, AdjGraph &graph) {
    int V;
    int id;
    std::string name;
    pin >> V ;
    for (size_t i = 0; i < V; i++) {
        pin >> id >> name;
        graph.insert_node(id, name);
    }
    int src, dst, weight;
    while(pin >> src >> dst >> weight) {
        graph.insert_edge(src, dst, weight);
    }
}


void AdjGraph::reverse_edge(EdgeKey &e) {
    std::swap(Edges[e].src, Edges[e].dst);
}
void AdjGraph::reverse_edge_all() {
    // re
    /*Adjlist adj_copy(Adjs);
    for (auto &p : Nodes) {
        for (auto &q : Adjs[p.first]) {
            Adjs[p.first].erase(q);
            Adjs[q].insert(p.first);
        }
    }*/
}
inline void AdjGraph::init_single_source(int sid) {
    for (auto &p : Nodes) {
        p.second.init_single_node();
    }
    Nodes[sid].dist = 0;
    Nodes[sid].color = COLOR_GRAY;
    Nodes[sid].phi = nullptr;
}

inline int AdjGraph::get_weight(int src, int dst) {
    return edge_map[std::make_pair(src, dst)].weight;
}
// relax edge(u, v) used for finding shortest path 
inline void AdjGraph::relax_distance(int uid, int vid, int weight) {
    NodeType *v = get_node(vid);
    NodeType *u = get_node(uid);
    if(v->dist > u->dist + weight) {
        v->dist = u->dist + weight;
        v->phi = u;
    }
}
bool comp_pair(std::pair<int*, int> &lhs, std::pair<int*, int> &rhs) {
    if (*lhs.first > *rhs.first) {
        return true;
    } else {
        return false;
    }
} 
bool compare_by_dist(NodeType *lhs, NodeType *rhs) {
    if (lhs->dist > rhs->dist) return true;
    else return false;
}
inline void AdjGraph::relax_distance(int uid, int vid, int weight, FakeHeap &heap) {
    NodeType *v = get_node(vid);
    NodeType *u = get_node(uid);
    if(v->dist > u->dist + weight) {
        v->dist = u->dist + weight;
        v->phi = u;
        std::make_heap(heap.begin(), heap.end(), compare_by_dist);
    }
}
// Heap version, O(ElgV);
void AdjGraph::dijkstra_clrs(int src) {
    init_single_source(src);
    std::set<int> S; // actually, it's useless
    // v.dist represent assumed shortest path from source
    // (value type of vector also can be <NodeType *>)(maybe more elegant version)
    FakeHeap PQ;
    for (auto &p : Nodes) {
        PQ.push_back(&p.second);
    }
    std::make_heap(PQ.begin(), PQ.end(), compare_by_dist);
    while (!PQ.empty()) {
        // retrival root, which is shortest distance from src
        int u = (PQ.front())->id;
        // send it to back
        std::pop_heap(PQ.begin(), PQ.end());
        // and remove back
        PQ.pop_back();
        // reconstruct to recover heap property
        //std::make_heap(PQ.begin(), PQ.end(),compare_by_dist);
        // add u to visited set
        S.insert(u);
        for (auto &v : Adjs[u]) {
            relax_distance(u, v, get_weight(u, v), PQ);
        }
    }
}
Matrix AdjGraph::floyd_warshall(int src) {
    if (Mats.size() == 0) {
        init_matrix();
    }
    Matrix fw = Mats;
    int i, j, k;
    int V = Nodes.size();
    for (k = 0; k < V; k++) {
        for (i = 0; i < V; i++) {
            for (j = 0; j < V; j++) {
                if (fw[i][j] > fw[i][k] + fw[k][j]) {
                    fw[i][j] = fw[i][k] + fw[k][j];
                }
            }
        }
    }
    return fw;
}
bool AdjGraph::bellman_ford(int src) {
    init_single_source(src);
    
    // Iterate |V| - 1 through All edges
    for (size_t i = 0; i < Adjs.size() - 1; i++) {
        for (auto &e : Edges) {
            relax_distance(EDGE_SRC(e), EDGE_DST(e), EDGE_WEIGHT(e));
        }
        std::cout << "----" << std::endl;
        print_node_dist(std::cout);
    }
    for (auto &e : Edges) {
        // each edge(u, v) in E
        NodeType *u = get_node(EDGE_SRC(e));
        NodeType *v = get_node(EDGE_DST(e));
        int weight = EDGE_WEIGHT(e);
        if (v->dist > u->dist + weight) {
            return false;
        }
    }
    return true;
}
inline bool AdjGraph::check_path_to(int src, int dst) {
}
// initialize adj mat as adj list. 
void AdjGraph::init_mat_to_lst(Matrix &adj_mat) {
    const unsigned N = adj_mat.size();
    for (unsigned i = 0; i < N; ++i) {
        for (unsigned j = 0; j < N; ++j) {
            if(adj_mat[i][j] != 0) {
                insert_edge(i+1, j+1, adj_mat[i][j]);
            }
        }
    } 
}

inline void AdjGraph::init_all_nodes_property() {
    for (auto &p : Nodes) {
        p.second.init_single_node();
    }
}
nForest AdjGraph::dfs_forest(int sid) {
    nForest forest;
    init_single_source(sid);
}
// it doesn't work if graph is dis-connected
nodePath AdjGraph::dfs_stack(int sid, nVisitor visitor = nullptr) {
    //if(Adjs[sid].size() == 0) return nodePath();
    init_single_source(sid);
    nodePath vo; // save visit order
    std::stack<NodeType *> after_visits;
    int _time = 0;
    get_node(sid)->dist = _time;
    after_visits.push(get_node(sid));
    // Visiting smallest id as first is slightly poor performance
    // Because exploring node(i.e. root of forest) also can be pushed into a stack
    while (!after_visits.empty()) {
        NodeType *u = after_visits.top();
        after_visits.pop();
        // because it pushes every adjacent & grey node to visit the node with smallest id as first
        if(u->color == COLOR_BLACK) {
            continue;
        } 
        u->color = COLOR_BLACK;
        // not correct finish_time in stack based dfs 
        //u->finish_time = _time;
        vo.push_back(u);
        // visit smallest id as first. 
        // Adjs is sorted increasingly, and stack -> reverse visit
        auto rev_p = Adjs[u->id].rbegin();
        for(; rev_p != Adjs[u->id].rend(); ++rev_p) {
            NodeType *v = get_node(*rev_p);
            if (v->color != COLOR_BLACK) {
                v->phi = u;
                v->color = COLOR_GRAY;
                // dist is same meaning as depth
                v->dist = v->phi->dist + 1;
                after_visits.push(v);
            }
        }
        _time++;
    }
    return vo;
}
// topological_sort is possible
// only if given graph is DAG(direct acyclic graph)
nodePath AdjGraph::dfs_recursive(int sid, nVisitor *visitor = nullptr) {
    init_single_source(sid);
    nodePath path;
    int _time = 0;
    dfs_visit(get_node(sid), _time, path);
    for (auto &p : Nodes) {
        if (p.second.color == COLOR_WHITE) {
            dfs_visit(&p.second, _time, path);
        }
    }   
    return path;
}
// also visit small id first 
void AdjGraph::dfs_visit(NodeType *u, int &_time, nodePath &path) {
    _time++;
    u->dist = _time;
    u->color = COLOR_GRAY;
    for (auto &p : Adjs[u->id]) {
        NodeType *v = get_node(p);
        if (v->color == COLOR_WHITE) {
            v->phi = u;
            dfs_visit(v, _time, path);
        }
    }
    u->color = COLOR_BLACK;
    _time++;
    u->finish_time = _time;
    path.push_front(u);
}
nodePath AdjGraph::tsort_kahn() {
    std::queue<NodeType *> q;
    std::vector<int> in_degs(Nodes.size());
    int cnt = 0;

    for (auto &v : Nodes) {
        in_degs[cnt++] = v.second.in_deg;
        if (v.second.in_deg == 0)
            q.push(&v.second);
    }
    nodePath kahn_order;
    cnt = 0;
    while (!q.empty()) {
        NodeType* u = q.front();
        q.pop();
        kahn_order.push_back(u);
        for (auto &p : Adjs[u->id]) {
            NodeType *v = get_node(p);
            if (--(v->in_deg) == 0) {
               q.push(v); 
            }
        }
        ++cnt;
    }
    if (cnt != Nodes.size()) {
        std::cerr << "There exists a cycle in the graph\n";
    }
    cnt = 0;
    // recover in-deg
    for (auto &v : Nodes) {
        v.second.in_deg = in_degs[cnt++];
    }
    return kahn_order;
}
nodePath AdjGraph::topological_sort() {
    return dfs_recursive(1, nullptr);
}
nodePath AdjGraph::bfs(int sid, nVisitor visitor = nullptr) {
    init_all_nodes_property();
    NodeType *start_node = get_node(sid);
    start_node->color = COLOR_GRAY;
    start_node->dist = DIST_INF;
    start_node->phi = nullptr;
    queue<NodeType *> after_visits;
    after_visits.push(start_node);
    nodePath vo; // save visit order
    while(!after_visits.empty()) {
        NodeType *u = after_visits.front();
        after_visits.pop();
        vo.push_back(u);
        for (auto p : Adjs[u->id]) {
            NodeType *v = get_node(p);
            if (v->color == COLOR_WHITE) {
                v->color = COLOR_GRAY;
                v->dist = u->dist + 1;
                v->phi = get_node(u->id);
                after_visits.push(v); // copy by value
            }
        }
        u->color = COLOR_BLACK;
    }
    return vo;
}
inline bool AdjGraph::is_already(int id) {
    if(Nodes.find(id) != Nodes.end()) return true;
    else return false;
}

NodeType* AdjGraph::insert_node(int id) {
    if (!is_already(id)) {
        Nodes.insert(make_pair(id, NodeType(id)));
    }
    return &Nodes[id];
}
NodeType* AdjGraph::insert_node(int id, std::string name) {
    if (!is_already(id)) {
        Nodes.insert(make_pair(id, NodeType(id, name)));
    }
    return &Nodes[id];
}
NodeType* AdjGraph::get_node(int id) {
    return &Nodes[id];
}
// when grpah input is given as pair of number
EdgeType& AdjGraph::insert_edge(int src, int dst) {
    insert_node(src)->out_deg++;
    insert_node(dst)->in_deg++;

    Adjs[src].insert(dst);
    std::pair<int, int> end_point = std::make_pair(src, dst);
    EdgeType edge(src, dst);
    edge_map.insert(std::make_pair(end_point, edge));

    if (!directed) {
        Adjs[dst].insert(src);
        std::pair<int, int> end_point = std::make_pair(dst, src);
        EdgeType edge(dst, src);
        edge_map.insert(std::make_pair(end_point, edge));
    }
    return edge_map[end_point];
}
void AdjGraph::insert_edge(int src, int dst, int weight) {
    this->insert_edge(src, dst).weight = weight;
    if( !directed) {
        this->insert_edge(dst, src).weight = weight;
    }
}
void AdjGraph::delete_edge(EdgeKey &e) {
    get_node(EDGE_S(e))->out_deg--;
    get_node(EDGE_D(e))->in_deg--;
    Adjs[EDGE_S(e)].erase(EDGE_D(e));
    Edges.erase(e);
}
void AdjGraph::delete_node(int id) {
    Nodes.erase(id);
    Adjs.erase(id);
    int V = Adjs.size();
    for (int u = 0; u < V; u++) {
        for (auto &v : Adjs[u]) {
            if (v == id) {
                Edges.erase({u+1, v});
                get_node(u)->out_deg--;
            }
        }
        Adjs[u].erase(id);
    }
}

void AdjGraph::print_matrix(std::ostream &out) {
    int V = Nodes.size();
    out << "adjacent matrix\n----------\n";
    for (int i = 0; i < V; ++i) {
        for (int j = 0; j < V; ++j) {
            if (Mats[i][j] == DIST_INF) 
                out << "* ";
            else 
                out << Mats[i][j] << " ";
        }
        out << "\n";
    } 
}
void gNode::print(std::ostream &pout) {
    pout << left << setw(4) << id 
         << setw(4) << in_deg
         << setw(4) << out_deg;
    if (phi == nullptr) {
        pout << setw(5) << "null";
    } else {
        pout << setw(5) << phi->id;
    }
    if (dist == DIST_INF) {
        pout << setw(5) << "inf";
    } else {
        pout << setw(5) << dist;
    }
    if (color == COLOR_GRAY) {
        pout << setw(8) << "grey";
    } else if(color == COLOR_BLACK) {
        pout << setw(8) << "black";
    } else if(color == COLOR_WHITE) {
        pout << setw(8) << "white";
    } else {
        pout << setw(8) << "none";
    }
    pout << setw(5) << finish_time << "\n";
}
void AdjGraph::print_nodes(std::ostream &pout) {
    pout << left << setw(4) << "id"
         << setw(4) << "in"
         << setw(4) << "out"
         << setw(5) << "phi"
         << setw(5) << "dist"
         << setw(8) << "color"
         << setw(5) << "fin\n";
    pout << "\n";
    for (auto &v : Nodes) {
        v.second.print(pout);
    }
}
void AdjGraph::print(std::ostream &pout) {
    pout << "AdjGraph Info : " << name << "\n";
    if(directed) pout << "Direct graph";
    else pout << "Indirect graph with";
    pout << " with (|V| = " << Nodes.size() << ", |E| = ";
    pout << edge_map.size() << ")" << std::endl;

    std::cout << std::endl;
    pout << "id\t" << "adj list" << endl; 
    for (const auto &p : Nodes) {
        pout << right << p.first << "\t-> ";
        for (const auto &q : Adjs[p.first]) {
            pout << q << " ";
        }
        pout << std::endl;
    }
    
}

void AdjGraph::print_node_dist(std::ostream &out) {
    for (const auto p : Nodes) {
        out << p.first << "(";
        out << p.second.dist << ")\n";
    }
}
void AdjGraph::print_dfs(std::ostream &out) {
    const int format_width = 6;
    out << "\nDepth First Search\n";
    out << "------------------------\n";
    out << left << "id\td\tf\n";
    for (const auto &p : Nodes) {
        out << p.first << "\t";
        if (p.second.dist == DIST_INF) {
            out << "Unreachable\n";
            continue;
        }
        out << p.second.dist << "\t"
            << p.second.finish_time << "\n"; 
    }
}
void AdjGraph::print_edges(std::ostream &out) {
    int cnt = 1;
    out << "Edge info : ";
    out << "|E| = " << edge_map.size() << "\n";
    for (const auto &p : edge_map) {
        out << "(" << get_node(p.first.first)->id << ", "
            << get_node(p.first.second)->id << ", "
            << p.second.weight << "),";
        out << "\n";
    }
}
void AdjGraph::mprint(std::ostream &pout) {
    for (size_t i = 0; i < adj_mat.size(); ++i) {
        for (size_t j = 0; j < adj_mat.size(); ++j) {
            std::cout << adj_mat[i][j];
        }
        std::cout << std::endl;
    } 
}
AdjGraph::AdjGraph(Matrix mat) {
    const size_t sz = mat.size();
    adj_mat.resize(sz);
    for (size_t i = 0; i < sz; i++) {
        adj_mat[i].resize(sz);
        for (size_t j = 0; j < sz; j++) {
            adj_mat[i][j] = mat[i][j];
        }
    }
}
// Fill 1 if path v to u is exist. If not, fill 0 instead.
Matrix AdjGraph::check_all_path_to() {
    if (adj_mat.size() == 0) {
        init_matrix();
    }
    const size_t sz = adj_mat.size();
    Matrix res(sz, vector<int>(sz));
    std::stack<int> ds;
    std::vector<bool> visited(sz, false);
    for (size_t i = 0; i < sz; i++) {
        ds.push(i);
        while (!ds.empty()) {
            int u = ds.top(); ds.pop();
            if(!visited[u]) {
                visited[u] = true;
                for (size_t j = 0; j < sz; j++) {
                    if (adj_mat[u][j] != 0) res[i][j] = 1 ;
                    if (!visited[j] && adj_mat[u][j] != 0) {
                        ds.push(j);
                        //cout << j+1 << " ";
                    }
                }
            }
        }
        std::fill(visited.begin(), visited.end(), false);
    }
    return res;
}

inline bool AdjGraph::check_cycle(int src) {
    return check_path_to(src, src);
}

// initialize graph from file, which representation of edge is u v w
void AdjGraph::init_from_weighted_edge(istream &pin) {
    string line("");
    int u, v, w;
    while(getline(pin, line)) {
        stringstream iss(line);
        if(line[0] == '#') {
            std::cout << line << std::endl;
            continue;
        }
        // given only one node
        if (line.size() < 2) {
            iss >> u;
            insert_node(u);
        } else {
            iss >> u >> v >> w;
            insert_edge(u, v, w);
        }
    }
}
