#include "k_core.h"

#include <queue>
#include <ranges>
#include <numeric>
#include <iostream>


namespace kcm
{

KCore &KCore::remove_edge(Node u, Node v)
{
    auto w_opt = graph->get_edge(u, v);
    if(!w_opt.has_value()) {
        return *this;
    }
    auto w = w_opt.value();

    if(k_core_number[u] > k_core_number[v]) {
        std::swap(u, v);
    }

    Graph subcore(0);
    std::vector<Weight> cd;
    std::vector<bool> visisted;
    if(k_core_number[v] - k_core_number[u] >= w) {
        graph->remove_edge(u, v);
        std::tie(subcore, visisted, cd) = remove_subcore(u, w);
    } else {
        std::tie(subcore, visisted, cd) = remove_subcore(u, w);
        graph->remove_edge(u, v);
        cd[u] = cd[u] - w;
        cd[v] = cd[v] - w;
        subcore.remove_edge(u, v);
    }

    // update with subcore:
    std::vector<Node> ordered_nodes(graph->size(), 0);
    std::iota(ordered_nodes.begin(), ordered_nodes.end(), 0);

    std::stable_sort(ordered_nodes.begin(), ordered_nodes.end(), [=](Node i, Node j){
        return cd[i] < cd[j];
    });

    Weight a = k_core_number[u], l = -w;

    for(auto u : ordered_nodes) {
        if(!visisted[u]) {
            continue;
        }
        visisted[u] = false;

        for(auto [v, w] : subcore.neighbors(u)) {
            if(cd[v] > cd[u]) {
                cd[v] = cd[v] - w;
            }
        }

        if(cd[u] <= a + l) {
            k_core_number[u] = a + l;
        } else {
            l = k_core_number[u] - a;
            k_core_number[u] = cd[u];            
        }

        std::stable_sort(ordered_nodes.begin(), ordered_nodes.end(), [=](Node i, Node j){
            return cd[i] < cd[j];
        });
    }

    return *this;
}


KCore &KCore::add_edge(Node u, Node v, Weight w)
{
    graph->add_edge(u, v, w);

    if(k_core_number[u] > k_core_number[v]) {
        std::swap(u, v);
    }
    
    auto [subcore, visisted, cd] = add_subcore(u, w);

    std::vector<Node> ordered_nodes(graph->size(), 0);
    std::iota(ordered_nodes.begin(), ordered_nodes.end(), 0);

    // increasing order
    std::stable_sort(ordered_nodes.begin(), ordered_nodes.end(), [=](Node i, Node j){
        return cd[i] < cd[j];
    });


    Weight a = k_core_number[u], l = 0;

    // update with subcore:
    for(auto u : ordered_nodes) {
        if(!visisted[u]) {
            continue;
        }
        visisted[u] = false;
        for(auto [v, w] : subcore.neighbors(u)) {
            if(cd[v] > cd[u]) {
                cd[v] = cd[v] - w;
            }
        }

        if(cd[u] <= a + l) {
            k_core_number[u] = a + l;
        } else {
            l = k_core_number[u] - a;
            k_core_number[u] = cd[u];
        }

        std::stable_sort(ordered_nodes.begin(), ordered_nodes.end(), [=](Node i, Node j){
            return cd[i] < cd[j];
        });
    }
    
    return *this;
}


std::tuple<Graph, std::vector<bool>, std::vector<Weight>> KCore::remove_subcore(Node r, Weight w)
{
    std::vector<Weight> cd(graph->size(), 0);

    auto output = Graph(graph->size());

    auto visisted = graph->DFS(r, [&](Node node, Node neighbor, Weight weight) {
        if(k_core_number[neighbor] > k_core_number[r] - w) {
            cd[node] += weight;
            if(k_core_number[neighbor] <= k_core_number[r]) {
                output.add_edge(node, neighbor, weight);
                return true;
            }
        }
        return false; 
    });

    return {output, visisted, cd};
}


std::tuple<Graph, std::vector<bool>, std::vector<Weight>> KCore::add_subcore(Node r, Weight w)
{
    std::vector<Weight> cd(graph->size(), 0);
    
    auto output = Graph(graph->size());

    auto visisted = graph->DFS(r, [&](Node node, Node neighbor, Weight weight) {
        if(k_core_number[neighbor] >= k_core_number[r]) {
            cd[node] += weight;
            if(k_core_number[neighbor] < k_core_number[r] + w) {
                output.add_edge(node, neighbor, weight);
                return true;
            }
        }
        return false;
    });

    return {output, visisted, cd};
}


void KCore::calculate_initial_k_cores()
{
    // only cause std::priority queue donsn't support decrese-key
    std::vector<bool> removed(graph->size(), false);

    auto degrees = graph->degrees();

    // value: Degree / key: Node
    using HeapElement = std::tuple<Degree, Node>; 
    
    std::priority_queue<
        HeapElement, std::vector<HeapElement>, std::greater<HeapElement>
    > queue;

    for(auto node : graph->nodes()) {
        queue.push({degrees[node], node});
    }

    Degree k = 0;
    while(!queue.empty()) {
        auto [_, node] = queue.top(); 
        queue.pop();

        if(removed[node]) { // decrece key instead of this
            continue;
        }

        removed[node] = true;

        k = std::max(k, degrees[node]);        
        k_core_number[node] = k;

        for(auto &[neighbor, weight] : graph->neighbors(node)) {
            degrees[neighbor] -= weight;
            queue.push({degrees[neighbor], neighbor}); // decrece key instead
        }
    }
}

}
