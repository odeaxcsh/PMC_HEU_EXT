#include "pmc_heu.h"

#include <numeric>
#include <iostream>


namespace kcm 
{

std::vector<Node> find_heuristic_clique(const Graph &g, const KCore &k)
{
    std::vector<Node> ordered_nodes(g.size(), 0);
    std::iota(ordered_nodes.begin(), ordered_nodes.end(), 0);

    std::stable_sort(ordered_nodes.begin(), ordered_nodes.end(), [&](auto u, auto v){
        return k.get(u) > k.get(v);
    });


    std::vector<Node> H;
    Weight max_clique_size = 0;   
    for(auto u : ordered_nodes) {
        if(k.get(u) < max_clique_size) {
            continue;
        }

        std::vector<Node> C = { u };

        Weight C_k_core = 0; // k-core of G[C]
        std::vector<Weight> cd = { 0 }; // core numbers of G[C]

        auto S = g.neighbor_nodes(u);
        std::erase_if(S, [=](auto node){ return k.get(node) < max_clique_size; });
        std::sort(S.begin(), S.end(), [&](Node u, Node v){ 
            return k.get(u) > k.get(v); 
        });

        for(auto v : S) {
            Weight induced_degree = 0;
            bool extends_clique = true;

            for(auto u : C) {
                if(!g.get_edge(u, v).has_value()) {
                    extends_clique = false;
                    break;
                } else {
                    induced_degree += g.get_edge(u, v).value();
                }
            }

            if(!extends_clique || induced_degree < C_k_core) {
                continue;
            }

            for(int i = 0; i < C.size(); ++i) {
                cd[i] += g.get_edge(C[i], v).value();
            }

            C.push_back(v);
            cd.push_back(induced_degree);
            C_k_core = *std::min_element(cd.cbegin(), cd.cend());
        }

        Node least_core_number_vertex = *std::min_element(C.cbegin(), C.cend(), [&](const auto u, const auto v) {
            return k.get(u) < k.get(v);
        });

        if(C_k_core > max_clique_size) {
            std::cout << "Maximum Clique incraesed " << C_k_core << std::endl;
            max_clique_size = C_k_core;
            H = C;
        }
    }

    return H;
}

}
