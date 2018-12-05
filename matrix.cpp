#include <iostream>
#include <vector>
#include <fstream>
#include <sstream>
#include <utility>
#include <functional>


class invalid_input : public std::runtime_error {

public:

    invalid_input(const std::string & input) 
        : std::runtime_error("Invalid input: " + input) { }

};


class Matrix {

    std::vector<std::vector<double>> matrix;

public:

    Matrix() { }
    Matrix(const std::vector<std::vector<double>> & matrix) : matrix(matrix) { }

    void appendLine(const std::vector<double> & line) {

        if (matrix.size() and matrix[0].size() != line.size()) {
            throw std::runtime_error(std::string("Error: Row too ") + 
                                    (line.size() > matrix[0].size() ? "long" : "short"));
        }

        matrix.emplace_back(line);

    }

    const std::vector<std::vector<double>> & getMatrix() const {
        return matrix;
    }

    double determinant() const {
        return 0;
    }

};

void swapTwo(std::vector<std::vector<double>> & matrix, const size_t r1, const size_t r2) {

    const std::vector<double> v1 = matrix[r1];
    const std::vector<double> v2 = matrix[r2];

    matrix[r1] = v2;
    matrix[r2] = v1;

}

void swap(std::vector<std::vector<double>> & matrix, const size_t pivot) {

    if (matrix[pivot][pivot]) {
        return;
    }

    size_t nonzero = 0;

    for (size_t i = pivot + 1; i < matrix.size(); ++i) {
        if (matrix[i][pivot]) {
            nonzero = i;
            break;
        }
    }

    if (nonzero) {
        swapTwo(matrix, pivot, nonzero);
    }

}

void reduce(std::vector<std::vector<double>> & matrix, const size_t pivot) {

    if (not matrix[pivot][pivot]) {
        return;
    }

    for (size_t i = pivot + 1; i < matrix.size(); ++i) {
        if (not matrix[i][pivot]) {
            continue;
        }

        const double ratio = matrix[i][pivot] / matrix[pivot][pivot];
        for (size_t column = pivot; column < matrix[0].size(); ++column) {
            matrix[i][column] -= matrix[pivot][column] * ratio;
        }
    }

}

Matrix gauss(const Matrix & matrix) {

    std::vector<std::vector<double>> m = matrix.getMatrix();

    const size_t rows = m.size();
    const size_t columns = m[0].size();

    std::cout << "rows: " << rows << " columns: " << columns << std::endl;

    if (rows + 1 != columns) {
        throw std::runtime_error("Cannot perform gaussian method on current matrix");
    }

    for (size_t pivot = 0; pivot < rows - 1; ++pivot) {

        swap(m, pivot);
        reduce(m, pivot);

    }

    return Matrix(m);

}

std::ostream & operator<<(std::ostream & os, const Matrix & matrix) {

    for (const auto vector : matrix.getMatrix()) {
        for (const double i : vector) {
            os << i << " ";
        }
        std::endl(os);
    }

    return os;

}

void checkLineIsValid(const std::string & line) {

    bool dashAllowed = true;

    for (const char c : line) {

        if (c == '-' and dashAllowed) {
            dashAllowed = false;
            continue;
        }

        if (isspace(c)) {
            dashAllowed = true;
            continue;
        }

        if (not isdigit(c)) {
            throw invalid_input(line);
        }

    }

}

std::vector<double> parseLine(const std::string & line) {

    checkLineIsValid(line);

    std::vector<double> vector;
    std::string number;

    for (const char c : line) {

        if (isspace(c) and number.size()) {

            vector.emplace_back(std::stod(number));
            number = "";

        } else {
            number.append(1, c);
        }

    }

    if (number.size()) {
        vector.emplace_back(std::stod(number));
    }

    return vector;

}

bool isWSOnly(const std::string & str) {

    for (const char c : str) {
        if (not isspace(c)) {
            return false;
        }
    }

    return true;

}

Matrix readMatrix(std::istream & is) {

    Matrix matrix;

    while (not is.eof()) {

        std::string input;
        std::getline(is, input);

        if (isWSOnly(input)) {
            break;
        }

        matrix.appendLine(parseLine(input));

    }

    return matrix;

}

void matrixFromStdin() {

    std::cout << "Input matrix (empty line denotes end of matrix): " << std::endl;

    Matrix matrix = readMatrix(std::cin);

    std::cout << matrix;

}

void matrixFromFile(const std::string & filename) {

    std::ifstream input(filename);

    if (not input) {
        throw std::runtime_error("Error opening file: " + filename);
    }

    Matrix matrix = readMatrix(input);
    std::cout << matrix;
    Matrix transformed = gauss(matrix);
    std::cout << transformed;

}


int main(int argc, const char * argv[]) {

    std::vector<std::string> args;

    matrixFromFile("matrix.txt");
    return 0;

    for (int i = 1; i < argc; ++i) {
        args.emplace_back(argv[i]);
    }

    if (not args.size()) {

        try {
            matrixFromStdin();
        } catch (const std::exception & e) {
            std::cout << e.what() << std::endl;
        }
        return -1;

    }

    for (const std::string & arg : args) {

        try {
            matrixFromFile(arg);
        } catch (const std::exception & e) {
            std::cout << e.what() << std::endl;
        }

    }

}
