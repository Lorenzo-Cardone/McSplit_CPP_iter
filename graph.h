#include <limits.h>
#include <stdbool.h>

#include <vector>
#include <unordered_map>
#include <unordered_set>

struct Graph {
    size_t n;
    std::vector<std::unordered_map<size_t, unsigned int>> adjset; // the map is a map<other_node_id, edge_label>
    std::vector<unsigned int> label;
    Graph(unsigned int n);
    void printGraphMtx();

    std::vector<std::unordered_map<long int, size_t>> neighboring_labels; // map<label, amount>
    void get_neighboring_nodes (std::unordered_set<size_t> &node_set, size_t node_id, size_t distance); // returns a set<node_id>
    std::unordered_map<long int, size_t> get_neighboring_labels (size_t node_id, size_t distance); // returns a map<label, amount>
    void initialize_neighboring_labels (size_t distance); // initializes the neighboring_labels vector with all the maps
};

Graph induced_subgraph(struct Graph& g, std::vector<int> vv);

Graph readGraph(char* filename, char format, bool directed, bool edge_labelled, bool vertex_labelled);

Graph graphFromMtx(std::vector<std::unordered_map<size_t, unsigned int>> mat);