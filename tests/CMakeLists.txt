add_executable(
    graph
    graph.cpp
    k_core.cpp
)

target_link_libraries(
    graph
    GTest::gtest_main
    heulib
)

include(GoogleTest)
gtest_discover_tests(graph)
