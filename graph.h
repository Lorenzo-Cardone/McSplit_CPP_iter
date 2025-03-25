#include <limits.h>
#include <stdbool.h>

#include <vector>

struct Graph {
    int n;
    std::vector<std::vector<unsigned int>> adjmat;
    std::vector<unsigned int> label;
    Graph(unsigned int n);
    void printGraphMtx();
};

Graph induced_subgraph(struct Graph& g, std::vector<int> vv);

Graph readGraph(char* filename, char format, bool directed, bool edge_labelled, bool vertex_labelled);

Graph graphFromMtx(std::vector<std::vector<unsigned int>> mat);