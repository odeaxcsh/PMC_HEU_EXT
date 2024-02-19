#pragma once

#include <memory>

#include "graph.h"



namespace kcm
{

class KCore
{
public:    
    KCore(Graph &g) : graph(g), k_core_number(graph.size()) 
    {
        calculate_initial_k_cores();
    };


    KCore &remove_edge(Node u, Node v);

    KCore &add_edge(Node u, Node v, Weight w);
    

    Weight get(Node u) const
    {
        return k_core_number[u];
    }


private:
    
    std::tuple<Graph, std::vector<bool>, std::vector<Weight>> remove_subcore(Node r, Weight w);
    
    std::tuple<Graph, std::vector<bool>, std::vector<Weight>> add_subcore(Node r, Weight w);

    void calculate_initial_k_cores();
    

private:
    Graph &graph;
    std::vector<Weight> k_core_number;
};

}

