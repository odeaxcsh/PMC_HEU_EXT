#include <gtest/gtest.h>

#include "graph.h"
#include <string>


using namespace kcm;

TEST(Graph, Edges) {
    kcm::Graph g(3);
    g.add_edge(0, 1, 3);
    g.add_edge(0, 2, 5);

    ASSERT_EQ(g.get_edge(0, 1).value(), 3);
    ASSERT_EQ(g.get_edge(0, 2).value(), 5);

    ASSERT_FALSE(g.get_edge(1, 2).has_value()) << "Edge shouldn't be there";
}


TEST(Graph, WDegree) {
    kcm::Graph g(3);
    g.add_edge(0, 1, 3);
    g.add_edge(0, 2, 5);

    ASSERT_EQ(g.wdegree(0), 8);
    ASSERT_EQ(g.wdegree(1), 3);
    ASSERT_EQ(g.wdegree(2), 5);
}


TEST(Graph, WDegreeComplete) {
    kcm::Graph g(5);
    for(Node i = 0; i < 5; ++i) {
        for(Node j = 0; j < i; ++j) {
            g.add_edge(i, j, Weight(0.1));
        }
    }

    for(Node i = 0; i < 5; ++i) {
        ASSERT_EQ(g.wdegree(i), Weight(0.4)) << "Failed at " << std::to_string(i);
    }
}


TEST(Graph, RemoveEdge)
{
    kcm::Graph g(3);
    g.add_edge(0, 1, 3);
    g.add_edge(0, 2, 5);

    ASSERT_EQ(g.get_edge(0, 2).value(), 5);
    g.remove_edge(0, 2);
    ASSERT_FALSE(g.get_edge(0, 2).has_value());


    ASSERT_EQ(g.get_edge(0, 1).value(), 3);
    g.remove_edge(0, 1);
    ASSERT_FALSE(g.get_edge(0, 1).has_value());
}
