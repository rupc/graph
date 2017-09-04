#include <iostream>
#include <iomanip>
#include <sstream>
#include "graph.h"
using namespace std;
const char input_romania[] = "data/romania.in";
const char input_mat2[] = "data/mat2.in";
const char input_bellmanford[] = "data/bellmanford2.in";
const char input_dijkstra[] = "data/dijkstra1.in";
const char input_dfs1[] = "data/dfs1.in";
const char input_topo1[] = "data/topo1.in";
const char input_topo2[] = "data/topo2.in";
void bellmanford_test(istream &pin);
void dijkstra_test(istream &pin);
void dfs_recur_test();
void topo_test();

ostream& operator << (ostream &out, nodePath path) {

    for (auto p : path) {
        out << p->id << "-";
    }
    out << "\n\n";
    return out;
}

int main(void) {
    AdjGraph graph1(false);
    std::ifstream pin(input_romania);
    if(!pin.is_open()) {
        std::cout << "Open Failed" << std::endl;
    }
    //init_romania(pin, graph1);
    topo_test();
    //dfs_recur_test();
    //graph1.print_edges(cout);
    //graph1.dijkstra(1, 0);
    //graph1.print_edges(cout);
    //bellmanford_test(pin);
    //dijkstra_test(pin);
    return 0;
}

void dfs_recur_test() {
    AdjGraph DfsRecursive(true);
    std::ifstream pin(input_dfs1);
    DfsRecursive.init_from_weighted_edge(pin);
    //DfsRecursive.print(cout);
    //DfsRecursive.dfs_recursive(3, nullptr);
    cout << DfsRecursive.topological_sort();

    //DfsRecursive.print_dfs(cout);
}

void topo_test() {
    AdjGraph topoGraph(true);
    ifstream pin(input_topo2);
    topoGraph.init_from_weighted_edge(pin);
    cout << topoGraph.tsort_kahn();
    topoGraph.print(cout);

    topoGraph.init_matrix();
    topoGraph.print_matrix(cout);
    //topoGraph.print(cout);
    //cout << topoGraph.dfs_recursive(5, nullptr);
    //topoGraph.print_nodes(cout);
    //topoGraph.print_dfs(cout);
}

void dijkstra_test(istream &pin) {
    AdjGraph DijkstraGraph(false);
    DijkstraGraph.init_from_weighted_edge(pin);
    DijkstraGraph.print(cout);
    DijkstraGraph.reverse_edge_all();
    /*DijkstraGraph.dijkstra_clrs(0);
    DijkstraGraph.print_node_dist(cout);*/
    DijkstraGraph.print(cout);
}

void bellmanford_test(istream &pin) {
    AdjGraph BellmanFordGraph(true);
    string line("");
    int u, v, w;
    while(getline(pin, line)) {
        stringstream iss(line);
        if(line[0] == '#') continue;
        iss >> u >> v >> w;
        BellmanFordGraph.insert_edge(u, v, w);
    }
    //BellmanFordGraph.print(cout);
    BellmanFordGraph.print_edges(cout);
    bool nCycle = BellmanFordGraph.bellman_ford(1);
    std::cout << "isCycle : " << boolalpha << nCycle << std::endl;
    BellmanFordGraph.print_node_dist(cout);
}
