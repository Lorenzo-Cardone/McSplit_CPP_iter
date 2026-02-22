#include "graph.h"

#include <stdio.h>
#include <stdlib.h>

#include <iostream>
#include <string>
#include <unordered_set>
#include <cmath>
#include <cstdint>

constexpr int BITS_PER_UNSIGNED_INT (CHAR_BIT * sizeof(unsigned int));

static void fail(std::string msg) {
    std::cerr << msg << std::endl;
    exit(1);
}

Graph::Graph(unsigned int n) {
    this->n = n;
    label = std::vector<unsigned int>(n, 0u);
    adjset = {n, std::unordered_map<size_t, unsigned int>()};
    label_count_per_node_fan_in = {n, std::unordered_map<unsigned int, float>()};
    label_count_per_node_fan_out = {n, std::unordered_map<unsigned int, float>()};
}

Graph induced_subgraph(struct Graph& g, std::vector<int> vv) {
    Graph subg(vv.size());
    #pragma omp parallel for
    for (size_t i=0; i<subg.n; i++) {
        //for (int j=0; j<subg.n; j++) {
        //    subg.adjmat[i][j] = g.adjmat[vv[i]][vv[j]];
        //}
        for (size_t j=0; j<subg.n; j++) {
            if (g.adjset[vv[i]].contains(vv[j])) {
                subg.adjset[i][j] = g.adjset[vv[i]][vv[j]];
            }
        }
    }

    for (size_t i=0; i<subg.n; i++)
        subg.label[i] = g.label[vv[i]];
    return subg;
}

void add_edge(Graph& g, int v, int w, bool directed=false, unsigned int val=1) {
    if (v != w) {
        if (directed) {
            g.adjset[v][w] |= val;
            g.adjset[w][v] |= (val<<16);
        } else {
            g.adjset[v][w] = val;
            g.adjset[w][v] = val;
        }
    } else {
        // To indicate that a vertex has a loop, we set the most
        // significant bit of its label to 1
        g.label[v] |= (1u << (BITS_PER_UNSIGNED_INT-1));
    }
}

struct Graph readDimacsGraph(char* filename, bool directed, bool vertex_labelled) {
    struct Graph g(0);

    FILE* f;
    
    if ((f=fopen(filename, "r"))==NULL)
        fail("Cannot open file");

    char* line = NULL;
    size_t nchar = 0;

    int nvertices = 0;
    int medges = 0;
    int v, w;
    int edges_read = 0;
    int label;

    while (getline(&line, &nchar, f) != -1) {
        if (nchar > 0) {
            switch (line[0]) {
            case 'p':
            {
                if ( sscanf(line, "p edge %d %d", &nvertices, &medges) != 2 )
                    fail("Error reading a line beginning with p.\n");
                g = Graph(nvertices);
                break;
            }
            case 'e':
            {
                unsigned int val = 1;
                int ret_val = sscanf(line, "e %d %d %d", &v, &w, &val);
                if ( ret_val !=2 && ret_val != 3) {
                    fail("Error reading a line beginning with e.\n");
                }
                add_edge(g, v-1, w-1, directed, val);
                edges_read++;
                break;
            }
            case 'n':
                if (sscanf(line, "n %d %d", &v, &label)!=2)
                    fail("Error reading a line beginning with n.\n");
                if (vertex_labelled)
                    g.label[v-1] |= label;
                break;
            }
        }
    }

    if (medges>0 && edges_read != medges) fail("Unexpected number of edges.");

    fclose(f);
    return g;
}

struct Graph readLadGraph(char* filename, bool directed) {
    struct Graph g(0);
    FILE* f;
    
    if ((f=fopen(filename, "r"))==NULL)
        fail("Cannot open file");

    int nvertices = 0;
    int w;

    if (fscanf(f, "%d", &nvertices) != 1)
        fail("Number of vertices not read correctly.\n");
    g = Graph(nvertices);

    for (int i=0; i<nvertices; i++) {
        int edge_count;
        if (fscanf(f, "%d", &edge_count) != 1)
            fail("Number of edges not read correctly.\n");
        for (int j=0; j<edge_count; j++) {
            if (fscanf(f, "%d", &w) != 1)
                fail("An edge was not read correctly.\n");
            add_edge(g, i, w, directed);
        }
    }

    fclose(f);
    return g;
}

int read_word(FILE *fp) {
    unsigned char a[2];
    if (fread(a, 1, 2, fp) != 2)
        fail("Error reading file.\n");
    return (int)a[0] | (((int)a[1]) << 8);
}

struct Graph readBinaryGraph(char* filename, bool directed, bool edge_labelled,
        bool vertex_labelled){
    struct Graph g(0);
    FILE* f;
    
    if ((f=fopen(filename, "rb"))==NULL)
        fail("Cannot open file");

    int nvertices = read_word(f);
    g = Graph(nvertices);

    // Labelling scheme: see
    // https://github.com/ciaranm/cp2016-max-common-connected-subgraph-paper/blob/master/code/solve_max_common_subgraph.cc
    int m = g.n * 33 / 100;
    int p = 1;
    int k1 = 0;
    int k2 = 0;
    while (p < m && k1 < 16) {
        p *= 2;
        k1 = k2;
        k2++;
    }
    
    for (int i=0; i<nvertices; i++) {
        int label = (read_word(f) >> (16-k1));
        if (vertex_labelled)
            g.label[i] |= label;
    }

    for (int i=0; i<nvertices; i++) {
        int len = read_word(f);
        for (int j=0; j<len; j++) {
            int target = read_word(f);
            int label = (read_word(f) >> (16-k1)) + 1;
            add_edge(g, i, target, directed, edge_labelled ? label : 1);
        }
    }
    fclose(f);
    return g;
}

struct Graph readGraph(char* filename, char format, bool directed, bool edge_labelled, bool vertex_labelled) {
    struct Graph g(0);
    if (format=='D') g = readDimacsGraph(filename, directed, vertex_labelled);
    else if (format=='L') g = readLadGraph(filename, directed);
    else if (format=='B') g = readBinaryGraph(filename, directed, edge_labelled, vertex_labelled);
    else fail("Unknown graph format\n");
    return g;
}

void Graph::printGraphMtx () {
    std::cout << "{";
    for (size_t x = 0; x < n; x++) {
        std::cout << "{";
        for (size_t y = 0; y < n; y++) {
            std::cout << (adjset[x].contains(y) ? adjset[x].at(y) : 0) << ((y != n-1) ? "," : "");
        }
        std::cout << "}" << ((x != n-1) ? "," : "");
    }
    //for (const auto &x : adjmat) {
    //    std::cout << "{";
    //    for (const auto &y : x) {
    //        std::cout << y << ((&y != &x.back()) ? "," : "");
    //    }
    //    std::cout << "}" << ((&x != &adjmat.back()) ? "," : "");
    //}
    std::cout << "}" << std::endl;
}

struct Graph graphFromMtx(std::vector<std::unordered_map<size_t, unsigned int>> mat) {
    struct Graph g(0);
    g.adjset = mat;
    g.n = mat.size();
    g.label.resize(g.n, 0u);
    return g;
}

enum Direction { FORWARD, BACKWARD, BOTH };

void computeNodeDesctriptors(Graph &g, int neighbourhood_radius, bool limit_fan_in_fan_out, bool fanout_only, float distance_effect_dampening)
{

    // for each node run BFS up to neighbourhood_radius
    #pragma omp parallel for
    for (size_t start_node = 0; start_node < g.n; start_node++) {
        std::vector<std::pair<uint64_t, std::pair<int, Direction>>> frontier; // node_idx, distance
        std::unordered_map<uint64_t, int> visited_backward; // node_idx, distance
        std::unordered_map<uint64_t, int> visited_forward; // node_idx, distance

        frontier.reserve(g.n);
        frontier.emplace_back(std::pair<uint64_t, std::pair<int, Direction>>{start_node, {0, BOTH}});

        size_t index = 0;

        while (index < frontier.size()) {
            std::pair<uint64_t, std::pair<int, Direction>> current = frontier[index++];

            int distance = current.second.first + 1;
            if (distance > neighbourhood_radius)
                continue;
            // add all neighbours of nodes at distance d to distance d+1
            for (auto neigh : g.adjset[current.first]) {
                size_t neigh_idx = neigh.first;
                unsigned int edge_label = neigh.second;
                // if limited fan in/fan out, check if current node is in forward or backward direction and add neighbour accordingly
                if (limit_fan_in_fan_out) {
                    if (edge_label & 0xFFFFu) { // forward edge
                        if (current.second.second == BACKWARD) {
                            continue; // cannot go forward from a node we reached backwards
                        }
                    } 
                    if (edge_label & 0xFFFF0000u) { // backward edge
                        if (current.second.second == FORWARD) {
                            continue; // cannot go backward from a node we reached forwards
                        }
                    }
                }
                // exit edge is val, enter edge is val << 16
                if ((edge_label & 0xFFFFu) || !limit_fan_in_fan_out) {
                    if (visited_forward.contains(neigh_idx)) {
                        continue;
                    }
                    visited_forward[neigh_idx] = distance;
                    frontier.emplace_back(std::pair<uint64_t, std::pair<int, Direction>>{neigh_idx, {distance, FORWARD}});
                }
                if ((edge_label & 0xFFFF0000u) && !fanout_only && limit_fan_in_fan_out) {
                    if (visited_backward.contains(neigh_idx)) {
                        continue;
                    }
                    visited_backward[neigh_idx] = distance;
                    frontier.emplace_back(std::pair<uint64_t, std::pair<int, Direction>>{neigh_idx, {distance, BACKWARD}});
                }
            }
        }

        // std::cout << "Index: " << index << std::endl;
        // std::cout << visited_forward.size() << std::endl;
        // std::cout << visited_backward.size() << std::endl;

        // visited all nodes, add them to label_count_per_node
        if (!fanout_only) {
            for (const auto & [node_idx, distance] : visited_backward) {
                unsigned int node_label = g.label[node_idx];
                float effect = std::pow(distance_effect_dampening, static_cast<float>(distance - 1));
                g.label_count_per_node_fan_in[start_node][node_label] += effect;
            }
        }
        for (const auto & [node_idx, distance] : visited_forward) {
            unsigned int node_label = g.label[node_idx];
            float effect = std::pow(distance_effect_dampening, static_cast<float>(distance - 1));
            g.label_count_per_node_fan_out[start_node][node_label] += effect;
        }
    }
    return;
}
