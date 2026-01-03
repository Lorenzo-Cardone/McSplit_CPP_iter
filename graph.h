#include <limits.h>
#include <stdbool.h>

#include <vector>
#include <unordered_map>

struct Graph {
    size_t n;
    std::vector<std::unordered_map<size_t, unsigned int>> adjset; // map of node_idx and edge_label
    std::vector<unsigned int> label;
    std::vector<std::unordered_map<unsigned int, float>> label_count_per_node;
    Graph(unsigned int n);
    void printGraphMtx();
};

Graph induced_subgraph(struct Graph& g, std::vector<int> vv);

Graph readGraph(char* filename, char format, bool directed, bool edge_labelled, bool vertex_labelled);

Graph graphFromMtx(std::vector<std::unordered_map<size_t, unsigned int>> mat);

void computeNodeDesctriptors(struct Graph& g, int neighbourhood_radius, float distance_effect_dampening, bool limit_fan_in_fan_out);