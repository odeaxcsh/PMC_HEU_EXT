#include "graph.h"
#include "k_core.h"
#include "pmc_heu.h"

#include <iostream>
#include <iomanip>
#include <fstream>

#include <chrono>
#include <numeric>
#include <cmath>
#include <set>


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


std::pair<kcm::Graph, std::vector<kcm::Node>> 
read_points(
    std::ifstream &reference, std::ifstream &sample, std::ifstream &gt, 
    float threshold = 0.05, float sigma = 1, bool binary = false
) {
    int n, m;
    reference >> n;
    sample >> m;

    int n_assoc = n * m;

    auto distance = [](std::pair<float, float> a, std::pair<float, float> b) {
        float dx = a.first - b.first;
        float dy = a.second - b.second;
        return std::sqrt(dx * dx + dy * dy);
    };

    std::vector<std::pair<float, float>> reference_points, sample_points;


    std::cout << "SAMPLE POINTS: " << m << std::endl;
    std::cout << "REFERENCE POINTS: " << n << std::endl;
    std::cout << "NUMBER OF ASSOCATIONS: " << n_assoc << std::endl;
    
    float x, y;
    for(int i = 0; i < n; ++i) {
        reference >> x >> y;
        reference_points.push_back({x, y});
    }

    for(int i = 0; i < m; ++i) {
        sample >> x >> y;
        sample_points.push_back({x, y});
    }

    auto g = kcm::Graph(n_assoc);

    auto start = std::chrono::high_resolution_clock::now();
    int count = 0;
    double sum = 0;
    for(int i = 0; i < n_assoc; ++i) {
        for(int j = 0; j < n_assoc; ++j) {
            float d = distance(sample_points[i % m], sample_points[j % m]);
            float d_p = distance(reference_points[i / m], reference_points[j / m]);
            float conssitency = std::exp(-sigma * std::abs(d - d_p));
            if(conssitency > threshold) {
                g.add_edge(i, j, binary ? 1 : conssitency);
                ++count;
                sum += binary ? 1 : conssitency;
            }
        }
    }
    auto end = std::chrono::high_resolution_clock::now();


    std::cout << "NUMBER OF EDGES: " << count << std::endl;
    std::cout << "GRAPH DENSITY: " << std::fixed << std::setprecision(4) << float(count) / (n_assoc * (n_assoc - 1)) << std::endl;
    std::cout << "WEIGHTED DENSITY: " << std::fixed << std::setprecision(4) << sum / (n_assoc * (n_assoc - 1)) << std::endl;

    std::cout << std::endl;

    std::cout << "Finding Maximum Clique took: " <<
        std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count() << " ms" << std::endl;
    std::cout << std::endl;

    int gt_n;
    gt >> gt_n;

    std::vector<kcm::Node> true_associations;

    for(int i = 0; i < gt_n; ++i) {
        gt >> x;
        true_associations.push_back(x);
    }
    return {g, true_associations};
    
}


std::pair<float, float> precition_recall(std::vector<kcm::Node> retrieved, std::vector<kcm::Node> relevant)
{
    std::set<kcm::Node> retrieved_set(retrieved.begin(), retrieved.end());
    std::set<kcm::Node> relevant_set(relevant.begin(), relevant.end());
    std::set<kcm::Node> intersect;

    set_intersection(
        retrieved_set.begin(), retrieved_set.end(), 
        relevant_set.begin(), relevant_set.end(),
        std::inserter(intersect, intersect.begin())
    );

    return {
        float(intersect.size()) / retrieved_set.size(),
        float(intersect.size()) / relevant_set.size()
    };
}


int main(int argc, char **argv)
{
    if(argc < 2) {
        std::cerr << "Usage: " << argv[0] << "input_file" << std::endl;
        return EXIT_FAILURE;
    }


    kcm::Graph g(0);
    std::vector<kcm::Node> true_association;

    if(argc == 2) {
        std::ifstream data_file;
        data_file.open(argv[1]);
        if(!data_file.is_open()) {
            std::cerr << "File " << argv[1] << "doesn't exist" << std::endl;
            return EXIT_FAILURE;
        }

        std::tie(g, true_association) = read_graph(data_file);
    }

    if(argc == 4) {
        std::ifstream reference, sample, gt;
        reference.open(argv[1]);
        if(!reference.is_open()) {
            std::cerr << "File " << argv[1] << "doesn't exist" << std::endl;
            return EXIT_FAILURE;
        }

        sample.open(argv[2]);
        if(!sample.is_open()) {
            std::cerr << "File " << argv[2] << "doesn't exist" << std::endl;
            return EXIT_FAILURE;
        }

        gt.open(argv[3]);
        if(!gt.is_open()) {
            std::cerr << "File " << argv[3] << "doesn't exist" << std::endl;
            return EXIT_FAILURE;
        }

        std::tie(g, true_association) = read_points(reference, sample, gt);
    }

    auto start = std::chrono::high_resolution_clock::now();
    kcm::KCore k(g);
    auto end = std::chrono::high_resolution_clock::now();

    std::cout << "Calculating K core numbers took: " <<
        std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count() << " ms" << std::endl;


    start = std::chrono::high_resolution_clock::now();
    auto H = kcm::find_heuristic_clique(g, k);
    end = std::chrono::high_resolution_clock::now();

    std::cout << "Finding Maximum Clique took: " <<
        std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count() << " ms" << std::endl;
    std::cout << std::endl;

    std::sort(H.begin(), H.end());

    std::sort(true_association.begin(), true_association.end());
    
    std::cout << "RESULT: ";
    for(int i = 0; i < H.size(); ++i) {
        if(i > 25) {
            break;
        }
        std::cout << H[i] << " ";

    }
    if(H.size() > 25) {
        std::cout << "... (" << H.size() - 25 << " more)" << std::endl;
    }
    std::cout << std::endl;

    std::cout << "GROUND: ";
    for(int i = 0; i < true_association.size(); ++i) {
        if(i > 25) {
            break;
        }
        std::cout << true_association[i] << " ";

    }
    if(true_association.size() > 25) {
        std::cout << "... (" << true_association.size() - 25 << " more)" << std::endl;
    }
    std::cout << std::endl;

    auto [precition, recall] = precition_recall(H, true_association);
    std::cout << "PRECISION: " << (int) (precition * 100) << " %" << std::endl; ;
    std::cout << "RECALL: " << (int) (recall * 100) << " %" << std::endl; ;
    std::cout << "RHO: " << std::fixed << std::setprecision(3) <<  1 - ((float) true_association.size() / g.size()) << std::endl;
    return EXIT_SUCCESS;
}
