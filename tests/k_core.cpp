#include <gtest/gtest.h>


#include "k_core.h"
#include "graph.h"


using namespace kcm;

TEST(KCore, Simple)
{
    auto g = std::make_shared<Graph>(5);
    for(Node i = 0; i < 4; ++i) {
        g->add_edge(i, i + 1, Weight(0.5));
    }

    KCore k(g);
    for(Node i = 0; i < 5; ++i) {
        ASSERT_EQ(k.get(i), Weight(0.5));
    }
}


TEST(KCore, CompleteGraph)
{
    auto g = std::make_shared<Graph>(5);
    for(Node i = 0; i < 5; ++i) {
        for(Node j = 0; j < i; ++j) {
            g->add_edge(i, j, Weight(0.1));
        }
    }

    kcm::KCore k(g);
    for(Node i = 0; i < 5; ++i) {
        ASSERT_EQ(k.get(i), Weight(0.4));
    }
}


TEST(KCore, AddEdge)
{
    auto g = std::make_shared<kcm::Graph>(5);
    for(Node i = 0; i < 4; ++i) {
        g->add_edge(i, i + 1, Weight(0.5));
    }

    kcm::KCore k(g);
    k.add_edge(0, 4, 0.5);
    for(Node i = 0; i < 5; ++i) {
        ASSERT_EQ(k.get(i), Weight(1.0)) << "failed at " << i;
    }
}


TEST(KCore, AddEdgeComplex)
{
    
}


TEST(KCore, RemoveEdge)
{

}
