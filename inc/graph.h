/*
 * Copyright (c) 2017 Yongrae, Jo <memex@postech.ac.kr>
 * Author: Yongrae, Jo <memex@postech.ac.kr>
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy of
 * this software and associated documentation files (the "Software"), to deal in
 * the Software without restriction, including without limitation the rights to
 * use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of
 * the Software, and to permit persons to whom the Software is furnished to do so,
 * subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS
 * FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR
 * COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER
 * IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
 * CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 */


#ifndef GRAPH_H
#define GRAPH_H
// input, output 
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
// container
#include <vector>
#include <list>
#include <map>
#include <set>
#include <string>
// container adaper
#include <stack>
#include <queue>
// STL algorithm related
#include <algorithm>
#include <functional>
#include <numeric>
// C standard library
#include <cmath>
#include <cstdlib>
#include <climits>

#define COLOR_WHITE 0x01 // for unexplored node
#define COLOR_GRAY  0x02 // for exploring node.
#define COLOR_BLACK 0x03 // for explored node
#define DIST_INF 999999999
#define EDGE_WEIGHT_DEFAULT 1

using namespace std;

struct gNode {
    // Followed property is dynamic
    // need to be initialze for graph traverse
    int color;
    int dist;
    gNode* phi;
    int finish_time = 0; 
    // Follwing property is static.
    // No need to initialze for search algorithms
    // except for insertion, deletion, modification ...
    int id; // Unique identifier by number or string name
    std::string name;
    int in_deg = 0;
    int out_deg = 0;
    gNode() = default;
    gNode(int node_id)
        : id(node_id), color(COLOR_WHITE),
          dist(DIST_INF), phi(nullptr), name("") {}
    gNode(int i, std::string s)
        : id(i), color(COLOR_WHITE), 
          dist(DIST_INF), phi(nullptr), name(s) {}
    void init_single_node() {
        color = COLOR_WHITE;
        dist = DIST_INF;
        phi = nullptr;
        finish_time = 0;
    }
    void print(ostream &out);
};

struct gEdge {
    int src;
    int dst;
    int weight;
    gEdge() = default;
    gEdge(int s, int d, int w) 
        : src(s), dst(d), weight(w) {}
    gEdge(int s, int d)
        : src(s), dst(d), weight(EDGE_WEIGHT_DEFAULT) {}
};

bool act_on_node(gNode * v);
typedef gNode NodeType;
typedef gEdge EdgeType;
// id values of visit order
typedef std::vector<int> iPath; 
// node pointers of visit order
typedef std::list<NodeType *> nodePath;    
typedef std::vector<int> iPath; 
typedef std::list<nodePath> nForest;
// action of id based traverse
using iActor = std::function<bool(int)>; 
// visit and work on node while traverse
typedef std::function<bool(NodeType *)> nVisitor; 

typedef std::vector<std::vector<int>> Matrix;
typedef std::vector<std::vector<int>> AdjMat;
typedef std::map<int, std::set<int>> AdjList;

typedef std::map<int, NodeType> IdMap;
typedef std::set<int> NodeSet;
// Collection of eges
typedef std::pair<int, int> EdgeKey;
typedef std::map<EdgeKey, EdgeType> EdgeMap;
// shortest path
typedef std::vector<int> shortDists; 
typedef std::vector<std::pair<int*, int>> Heap;
typedef std::vector<NodeType *> FakeHeap;
/*
 * Fundamental Graph Algorithms 
 *     1. Depth First Search
 *     2. Breadth First Search
 *     3. Dijkstra 
 *     4. Bellman-Ford Algorithm
 *     5. Floyd Washall Algorithm
 *     5. Topological Sorting 
 *     6. Find the exsistence of path from source to other vertexes
 *     7. Cycle Detection Algorithm
 */
// get edge property of EdgeMap
#define EDGE_SRC(e) (e.first.first)
#define EDGE_DST(e) (e.first.second)
#define EDGE_WEIGHT(e) (e.second.weight)
// get id from EdgeKey
#define EDGE_S(e) (e.first)
#define EDGE_D(e) (e.second)
/*class _Heap {

};*/
class AdjGraph {
    // ailias of internal data structure for neat name
    IdMap &Nodes = id_node_map;
    AdjList &Adjs = id_adj_map;
    EdgeMap &Edges = edge_map;
    AdjMat &Mats = adj_mat;
    public:
        AdjGraph() = default; // default constructor
        // copy constructor
        AdjGraph(const AdjGraph&);
        AdjGraph(bool dir, std::string s = "") 
            : directed(dir), name(s) {
        }
        AdjGraph(int V, int E, bool dir, bool allow_minus) {}
        AdjGraph(Matrix mat);
        ~AdjGraph() {}
        // initialize internal data structures
        inline void init_all_nodes_property();
        void init_mat_to_lst(Matrix &mat);
        void init_matrix();
        
        void init_from_weighted_edge(istream &pin);
        inline void init_single_source(int sid);

        inline bool is_already(int id);
        // insertion function
        NodeType* insert_node(int id);
        NodeType* insert_node(int id, std::string name);
        EdgeType& insert_edge(int src, int dst);
        void insert_edge(int src, int dst, int weight);
        void delete_edge(EdgeKey &e);
        void delete_node(int id);
        // get internal data structure
        NodeType* get_node(int id);
        NodeSet* get_adjacent_node(int id);
        inline int get_weight(int src, int dst);
        size_t node_size() {return id_adj_map.size();}
        size_t edge_size() {return edge_map.size();}
        size_t out_degree(NodeType *v) {return Adjs[v->id].size();}
        size_t out_degree(int id) {return Adjs[id].size();}
        size_t in_degree(NodeType *v) {}
        size_t in_degree(int id) {}
        // fundamental graph algorithm
        // graph traverse
        nodePath bfs(int sid, nVisitor visitor);
        nodePath dfs_stack(int sid, nVisitor visitor);
        nForest dfs_forest(int sid);

        nodePath dfs_recursive(int sid, nVisitor *visitor);
        void dfs_visit(NodeType *u, int &_time, nodePath &path);
        // find shortest path
        void dijkstra_clrs(int src);
        bool bellman_ford(int src);
        Matrix floyd_warshall(int src);
        nodePath topological_sort();
        nodePath tsort_kahn(); 
        // test if path from src to dst exsit
        Matrix check_all_path_to();
        bool check_path_to(int src, int dst);
        bool is_acyclic();
        bool check_cycle(int src); // check cycle from src to src
        bool check_cycle();
        void reverse_edge(EdgeKey &e);
        void reverse_edge_all();
        bool has_path(int src, int dst);
        bool has_cycle();
        bool has_hamilton_cycle();
        bool has_euler_circuit();
        bool has_universal_sink();
        // aux function
        inline void relax_distance(int u, int v, int w);
        inline void relax_distance(int u, int v, int w, FakeHeap &heap);

        //inline void relax_distance(int u, int v, int w, Heap &heap);
        // pretty print function
        void print(std::ostream &pout);
        void print_edges(std::ostream &out);
        void print_shortest_path();
        void print_node_dist(std::ostream &out);
        void mprint(std::ostream &pout);
        void print_nodes(std::ostream &out);
        void print_dfs(std::ostream &out);
        void print_matrix(std::ostream &out);
        std::ostream& operator << (ostream &);
    private:
        std::string name;
        IdMap id_node_map; // id-node mapping table
        AdjList id_adj_map; // adjacent list table
        AdjMat adj_mat; // adjacent matrix
        // edge can be replicated ex) multi graph
        std::vector<EdgeType> edge_set;
        EdgeMap edge_map;
        bool directed;
};
void init_romania(std::istream &pin, AdjGraph &graph);
#endif
