#include <limits.h>
#include <stdbool.h>

#include <vector>
#include <unordered_map>

struct Graph {
    size_t n;
    std::vector<std::unordered_map<size_t, unsigned int>> adjset;
    std::vector<unsigned int> label;
    Graph(unsigned int n);
    void printGraphMtx();
};

Graph induced_subgraph(struct Graph& g, std::vector<int> vv);

Graph readGraph(char* filename, char format, bool directed, bool edge_labelled, bool vertex_labelled);

Graph graphFromMtx(std::vector<std::unordered_map<size_t, unsigned int>> mat);