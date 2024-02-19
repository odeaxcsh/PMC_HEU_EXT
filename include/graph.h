#pragma once

#include <map>
#include <tuple>
#include <vector>
#include <optional>
#include <ranges>
#include <queue>
#include <functional>



namespace kcm
{

using Weight = float;
using Degree = Weight;

using Node = unsigned int;
using NeighborsMap = std::map<Node, Weight>;


class Graph
{
public:
    Graph(size_t n) : neighborhoods(n), _degrees(n) {}


    size_t size() const;

    Graph &add_edge(Node u, Node v, Weight w);

    Graph &remove_edge(Node u, Node v);

    std::optional<Weight> get_edge(Node u, Node v) const;

    int degree(Node u) const;

    Degree wdegree(Node u) const;

    std::vector<Degree> degrees() const;

    std::vector<bool> DFS(Node root, std::function<bool(Node, Node, Weight)> on_edge);
    

    auto nodes() const
    {
        return std::views::iota(0, (int)size());
    }

    auto neighbors(Node u) const
    {
        return neighborhoods[u];
    }

    std::vector<Node> neighbor_nodes(Node u) const;


private:

    std::vector<NeighborsMap> neighborhoods;
    std::vector<Weight> _degrees;
};

}
