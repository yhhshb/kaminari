#include <iostream>
#include <fstream>
#include <vector>
#include <unordered_set>
#include <cmath>
#include <cassert>
#include <algorithm> // For std::min
#include "../../bundled/unordered_dense/include/ankerl/unordered_dense.h"

#include "../include/rbo.hpp"

namespace kaminari::rbo {

struct options_t {
    std::string file1;
    std::string file2;
    std::string output_filename;
    bool raw_lists;
};

options_t check_args(const argparse::ArgumentParser& parser);
double RBO(const std::vector<int>& A, const std::vector<int>& B, double p);
double RBO_MIN(const std::vector<int>& A, const std::vector<int>& B, double p);
double RBO_weight(double p, int d);

int main(const argparse::ArgumentParser& parser) 
{
    auto opts = check_args(parser);
    if (opts.raw_lists) {
        std::cerr << "Raw lists mode is not implemented yet" << std::endl;
        return 1;
    }

    std::cerr << "Reading lists from " << opts.file1 << " and " << opts.file2 << std::endl;
    std::cerr << "Output file: " << opts.output_filename << std::endl;

    std::ifstream file1(opts.file1);
    std::ifstream file2(opts.file2);
    std::ofstream output(opts.output_filename);

    if (!file1.is_open() || !file2.is_open() || !output.is_open()) {
        std::cerr << "Error opening files" << std::endl;
        return 1;
    }

    double sum_rbo = 0.0;

    std::string line1, line2;
    std::vector<int> list1, list2;
    std::vector<std::tuple<int, int>> values;
    
    //lambda func
    auto parse_line = [&](const std::string& line) {
        std::istringstream iss(line);
        std::string name;
        int list_size;
        iss >> name >> list_size;  // Read sequence name and list size

        std::string pair;
        while (iss >> pair) {
            int first, second;
            sscanf(pair.c_str(), "(%d,%d)", &first, &second);  // Parse "(x,y)" pairs
            values.emplace_back(first, second);
        }
    };

    size_t i = 0;
    for (i = 0; std::getline(file1, line1) && std::getline(file2, line2); ++i) {
        values.clear();
        parse_line(line1); //fill values

        int j = 0;
        while (j < values.size()) {
            std::vector<int> group;
            int score = std::get<1>(values[j]);
            group.push_back(std::get<0>(values[j]));
            j++;

            while (j < values.size() && std::get<1>(values[j]) == score) {
                group.push_back(std::get<0>(values[j]));
                j++;
            }

            std::sort(group.begin(), group.end());
            list1.insert(list1.end(), group.begin(), group.end());
        }

        values.clear();
        parse_line(line2); //fill values

        j = 0;
        while (j < values.size()) {
            std::vector<int> group;
            int score = std::get<1>(values[j]);
            group.push_back(std::get<0>(values[j]));
            j++;

            while (j < values.size() && std::get<1>(values[j]) == score) {
                group.push_back(std::get<0>(values[j]));
                j++;
            }

            std::sort(group.begin(), group.end());
            list2.insert(list2.end(), group.begin(), group.end());
        }


        double rbo = RBO(list1, list2, 0.9);
        list1.clear();
        list2.clear();
        sum_rbo += rbo;
        output << rbo << std::endl;
        
    }

    output << sum_rbo / i << std::endl;

    return 0;
}

argparse::ArgumentParser get_parser()
{
    argparse::ArgumentParser parser("rbo", "1.0.0", argparse::default_arguments::help);
    parser.add_description("Process RBO measure between lists (rbo -h for more information)");

    parser.add_argument("-f1", "--file1")
        .help("First file containing lists (1 list / line)")
        .required();
    parser.add_argument("-f2", "--file2")
        .help("Second file containing lists (1 list / line), this tool will process RBO between list in line n of file1 and line n of file2")
        .required();
    parser.add_argument("-o", "--output-filename")
        .help("RBOs output file, last line is average RBO values of all lists, others are RBO values of each list")
        .default_value("rbos.txt");
    parser.add_argument("--raw_lists") 
        .help("flag - lists are not considered to be output of index queries anymore, but raw lists of numbers separated by ',' or ';' or '\\t' [default: raw_lists is off]")
        .default_value(false)   
        .implicit_value(true);  

    return parser;
}

options_t check_args(const argparse::ArgumentParser& parser)
{
    options_t opts;
    opts.file1 = parser.get<std::string>("-f1");
    opts.file2 = parser.get<std::string>("-f2");
    opts.output_filename = parser.get<std::string>("-o");
    opts.raw_lists = parser.get<bool>("--raw_lists");

    return opts;
}


// Function to calculate RBO (Rank-Biased Overlap)
double RBO(const std::vector<int>& A, const std::vector<int>& B, double p) {
    
    assert(0.0 < p && p < 1.0);
    assert(ankerl::unordered_dense::set<int>(A.begin(), A.end()).size() == A.size()); // Check no duplicates
    assert(ankerl::unordered_dense::set<int>(B.begin(), B.end()).size() == B.size());

    size_t D = std::min(A.size(), B.size());
    ankerl::unordered_dense::set<int> A_prefix, B_prefix;
    double rbo = 0.0, P = 1.0;  // Current power of p
    int X_d = 0;  // Overlap at depth d

    for (size_t i = 0; i < D; ++i) {
        int d = i + 1;  // depth
        int overlap = (A_prefix.count(B[i]) + B_prefix.count(A[i]) + (A[i] == B[i]));
        X_d += overlap;
        assert(X_d == ankerl::unordered_dense::set<int>(A.begin(), A.begin() + d).size() 
                      + ankerl::unordered_dense::set<int>(B.begin(), B.begin() + d).size() 
                      - ankerl::unordered_dense::set<int>(A.begin(), A.begin() + d).count(B[i]));

        rbo += (X_d / (double)d) * P;  
        P *= p;  // Update power of p
        A_prefix.insert(A[i]);
        B_prefix.insert(B[i]);
    }
    return (1 - p) * rbo;
}

// Function to calculate RBO_MIN
double RBO_MIN(const std::vector<int>& A, const std::vector<int>& B, double p) {
    assert(0.0 < p && p < 1.0);
    assert(ankerl::unordered_dense::set<int>(A.begin(), A.end()).size() == A.size());
    assert(ankerl::unordered_dense::set<int>(B.begin(), B.end()).size() == B.size());

    size_t D = std::min(A.size(), B.size());
    ankerl::unordered_dense::set<int> A_prefix, B_prefix;
    double rbo_min = 0.0, P = 1.0;
    int X_d = 0;
    int X_D = 0;
    ankerl::unordered_dense::set<int> setA(A.begin(), A.begin() + D);
    ankerl::unordered_dense::set<int> setB(B.begin(), B.begin() + D);

    for (const auto& elem : setA) {
        if (setB.find(elem) != setB.end()) {
            ++X_D;  // Count common elements in the intersection
        }
    }

    for (size_t i = 0; i < D; ++i) {
        int d = i + 1;  // depth
        int overlap = (A_prefix.count(B[i]) + B_prefix.count(A[i]) + (A[i] == B[i]));
        X_d += overlap;
        assert(X_d == ankerl::unordered_dense::set<int>(A.begin(), A.begin() + d).size() 
                      + ankerl::unordered_dense::set<int>(B.begin(), B.begin() + d).size() 
                      - ankerl::unordered_dense::set<int>(A.begin(), A.begin() + d).count(B[i]));

        rbo_min += ((X_d - X_D) / (double)d) * P;  
        std::cerr << "i = " << i << ", X_d = " << X_d << ", X_D = " << X_D << ", rbo_min = " << rbo_min << std::endl;
        P *= p;  
        A_prefix.insert(A[i]);
        B_prefix.insert(B[i]);
    }
    return (1 - p) * (rbo_min - X_D * std::log(1 - p) / p);
}

// Function to calculate RBO_weight
double RBO_weight(double p, int d) {
    assert(0.0 < p && p < 1.0);
    double sum = 0.0, P = 1.0;
    for (int i = 0; i < d; ++i) {
        P *= p;
        sum += P / (i + 1);
    }
    return 1 - P + (1 - p) / p * d * (-std::log(1 - p) - sum);
}

} // namespace kaminari