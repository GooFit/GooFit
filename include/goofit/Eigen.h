#pragma once

#include <Eigen/Dense>
#include <fstream>
#include <string>

namespace GooFit {

template <typename M>
M read_csv(std::string name, bool comma = true) {
    using namespace Eigen;

    std::ifstream input(name);

    std::string line;
    std::vector<typename M::Scalar> values;
    size_t rows = 0;
    while(std::getline(input, line)) {
        std::stringstream lineStream(line);
        std::string cell;
        typename M::Scalar value;

        if(comma)
            while(std::getline<char>(lineStream, cell, ','))
                values.push_back(std::stod(cell));
        else
            while(lineStream >> value)
                values.push_back(value);

        ++rows;
    }
    return Map<const Matrix<typename M::Scalar, M::RowsAtCompileTime, M::ColsAtCompileTime, ColMajor>>(
        values.data(), values.size() / rows, rows);
};

} // namespace GooFit
