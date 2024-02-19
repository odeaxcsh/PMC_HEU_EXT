#include "graph.h"
#include <ranges>
#include <iostream>
#include <stack>
#include <ranges>


namespace kcm
{

size_t Graph::size() const
{
    return neighborhoods.size();
}


Graph &Graph::add_edge(Node u, Node v, Weight w)
{
    if(neighborhoods[u].find(v) != neighborhoods[u].end()) {
        remove_edge(u, v);
    }
    
    neighborhoods[u][v] = w;
    neighborhoods[v][u] = w;

    _degrees[u] += w;
    _degrees[v] += w;

    return *this;
}


Graph &Graph::remove_edge(Node u, Node v)
{
    auto it = neighborhoods[u].find(v);
    if(it != neighborhoods[u].end()) {
        _degrees[u] -= (*it).second;
        _degrees[v] -= (*it).second;

        neighborhoods[u].erase(it);
        
        // if one edge exists it's reverse must exist too
        neighborhoods[v].erase(neighborhoods[v].find(u));
    }
    return *this;
}


std::optional<Weight> Graph::get_edge(Node u, Node v) const
{
    auto it = neighborhoods[u].find(v);
    if(it == neighborhoods[u].end()) {
        return std::nullopt;
    }

    return (*it).second;
}


int Graph::degree(Node u) const
{
    return neighborhoods[u].size();
}


Degree Graph::wdegree(Node u) const
{
    return _degrees[u];
}


std::vector<Degree> Graph::degrees() const
{
    return _degrees;
}

std::vector<bool> Graph::DFS(Node root, std::function<bool(Node, Node, Weight)> on_edge)
{

    std::stack<Node> queue({root});
    std::vector<bool> visisted(size(), false);

    while(!queue.empty()) {
        Node u = queue.top();
        queue.pop();
        if(visisted[u]) {
            continue;
        }
        visisted[u] = true;

        for(auto [neighbor, weight] : neighbors(u)) {
            if(on_edge(u, neighbor, weight)) {
                queue.push(neighbor);
            }
        }
    }

    return visisted;
}

std::vector<Node> Graph::neighbor_nodes(Node u) const
{
    auto neighbors = std::views::keys(neighborhoods[u]);


    return std::vector<Node>(neighbors.begin(), neighbors.end());
}


}
