#include "graph.h"
#include "k_core.h"
#include "pmc_heu.h"

#include <iostream>
#include <fstream>
#include <numeric>


std::pair<kcm::Graph, std::vector<kcm::Node>> read_graph(std::ifstream &file)
{
    int n;
    file >> n;

    kcm::Graph g(n);

    float w;
    for(int i = 0; i < n; ++i) {
        for(int j = 0; j < n; ++j) {
            file >> w;
            if(w > 0) {
                g.add_edge(i, j, w);
            }
        }
    }

    file >> n;
    std::vector<kcm::Node> true_association;
    for(int i = 0; i < n; ++i) {
        file >> w;
        true_association.push_back(w);
    }

    return {g, true_association};
}


int main(int argc, char **argv)
{
    if(argc < 2) {
        std::cerr << "Usage: " << argv[0] << "input_file" << std::endl;
        return EXIT_FAILURE;
    }

    std::ifstream data_file;
    data_file.open(argv[1]);
    if(!data_file.is_open()) {
        std::cerr << "File " << argv[1] << "doesn't exist" << std::endl;
        return EXIT_FAILURE;
    }

    auto &&[g, true_association] = read_graph(data_file);

    kcm::KCore k(g);

    std::vector<kcm::Node> ordered_nodes(g.size(), 0);
    std::iota(ordered_nodes.begin(), ordered_nodes.end(), 0);
    std::stable_sort(ordered_nodes.begin(), ordered_nodes.end(), [&](auto u, auto v){
        return k.get(u) > k.get(v);
    });

    for(auto v : ordered_nodes) {
        std::cout << v << ": " << k.get(v) << std::endl;
    }

    auto H = kcm::find_heuristic_clique(g, k);
    std::sort(H.begin(), H.end());

    std::sort(true_association.begin(), true_association.end());


    std::cout << "RESULT: ";
    for(auto v : H) {
        std::cout << v << " ";
    }
    std::cout << std::endl;

    std::cout << "GROUND: ";
    for(auto v : true_association) {
        std::cout << v << " ";
    }
    std::cout << std::endl;

    return EXIT_SUCCESS;
}