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

class infinite_solutions : public std::runtime_error {

public:

    infinite_solutions()
        : std::runtime_error("Error: Matrix could not be solved because it has infinite number of solutions") { }

};

class no_solutions : public std::runtime_error {

public:

    no_solutions()
        : std::runtime_error("Error: Matrix has no solutions.") { }
};

typedef std::vector<std::vector<double>> matrix_t;

class Matrix {

    matrix_t matrix;

public:

    Matrix() { }
    Matrix(const matrix_t & matrix) : matrix(matrix) { }

    void appendLine(const std::vector<double> & line) {

        if (matrix.size() and matrix[0].size() != line.size()) {
            throw std::runtime_error(std::string("Error: Row too ") + 
                                    (line.size() > matrix[0].size() ? "long" : "short"));
        }

        matrix.emplace_back(line);

    }

    const matrix_t & getMatrix() const {
        return matrix;
    }

    size_t rows() const {
        return matrix.size();
    }

    size_t columns() const {
        
        if (not rows()) {
            return 0;
        }

        return matrix.front().size();

    }

    std::vector<double> & getLine(const size_t index) {
        return matrix[index];
    }

    double & get(const size_t line, const size_t index) {
        return matrix[line][index];
    }

    const std::vector<double> & getLine(const size_t index) const {
        return matrix[index];
    }

    double get(const size_t line, const size_t index) const {
        return matrix[line][index];
    }

    bool isSquare() const {

        if (not rows()) {
            return false;
        }

        return rows() == columns();

    }

    double determinant() const {
        
        if (not isSquare()) {
            throw std::runtime_error("Error: attempting to compute determinant of a non-square matrix");
        }

        double result = 1;

        for (size_t i = 0; i < rows(); ++i) {
            result *= get(i, i);
        }

        return result;

    }

    bool isSolvable() const {
        return rows() + 1 == columns();
    }

};

void swapTwo(matrix_t & matrix, const size_t r1, const size_t r2) {

    const std::vector<double> v1 = matrix[r1];
    const std::vector<double> v2 = matrix[r2];

    matrix[r1] = v2;
    matrix[r2] = v1;

}

void swap(matrix_t & matrix, const size_t pivot) {

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

void reduce(matrix_t & matrix, const size_t pivot) {

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

bool isZeroOnly(const std::vector<double> & line) {

    for (double d : line) {

        if (d) {
            return false;
        }

    }

    return true;

}

matrix_t removeZeroOnly(const matrix_t & matrix) {

    matrix_t res;

    for (const std::vector<double> & line : matrix) {

        if (not isZeroOnly(line)) {
            res.emplace_back(line);
        }

    }

    return res;

}

Matrix gauss(const Matrix & matrix) {

    matrix_t m = matrix.getMatrix();

    const size_t rowCount = matrix.rows();
    const size_t columnCount = m.front().size();

    for (size_t pivot = 0; pivot < rowCount - 1 and pivot < columnCount; ++pivot) {

        swap(m, pivot);
        reduce(m, pivot);

    }

    m = removeZeroOnly(m);

    return Matrix(m);

}

double * reallocate(double * arr, const size_t oldSize, const size_t newSize) {

    double * newArr = new double[newSize];

    for (size_t i = 0; i < oldSize; ++i) {
        newArr[i] = arr[i];
    }

    delete[] arr;

    return newArr;

}

std::vector<double> reverseAndConvert(const double * solved, const size_t sSize) {

    std::vector<double> retval;

    for (size_t i = sSize - 1;; --i) {

        retval.emplace_back(solved[i]);

        if (not i) {
            break;
        }

    }

    return retval;

}

double solveLine(const std::vector<double> line, const double * solved, const size_t sSize) {

    size_t lineIter = line.size() - 1;
    size_t solutionsIter = 0;

    double solution = line[lineIter--];

    while (solutionsIter < sSize) {
        solution -= solved[solutionsIter] * line[lineIter];

        --lineIter;
        ++solutionsIter;
    }

    if (not line[lineIter] and solution) {
        throw no_solutions();
    }
    if (not (line[lineIter] or solution)) {
        throw infinite_solutions();
    }

    return solution / line[lineIter];

}

std::vector<double> solve(const Matrix & matrix) {

    if (not matrix.isSolvable()) {
        throw std::runtime_error("");
    }

    size_t solutionsSize = 1;
    size_t solutionsCount = 0;
    double * solutions = new double[solutionsSize];

    for (size_t i = matrix.rows() - 1; i != (size_t)-1; --i) {

        if (solutionsCount == solutionsSize) {
            solutions = reallocate(solutions, solutionsSize, solutionsSize * 2);
            solutionsSize <<= 1;
        }

        solutions[solutionsCount] = solveLine(matrix.getLine(i), solutions, solutionsCount);

        ++solutionsCount;

    }

    std::vector<double> retval = reverseAndConvert(solutions, solutionsCount);
    delete[] solutions;
    return retval;

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

Matrix matrixFromStdin() {

    std::cout << "Input matrix (empty line denotes end of matrix): " << std::endl;

    Matrix matrix = readMatrix(std::cin);

    return matrix;

}

Matrix matrixFromFile(const std::string & filename) {

    std::ifstream input(filename);

    if (not input) {
        throw std::runtime_error("Error opening file: " + filename);
    }

    Matrix matrix = readMatrix(input);

    return matrix;
    

}

void writeMatrix(std::ofstream & os, const Matrix & matrix) {

    std::cout << matrix;

    os << "<table>\n";

    const matrix_t & innerMatrix = matrix.getMatrix();

    for (const auto & line : innerMatrix) {

        os << "<tr>";

        for (const double d : line) {
            os << "<td>" << d << "</td>";
        }

        os << "</tr>";
    }

    os << "</table>" << std::endl;

}

void writeDeterminant(std::ofstream & os, const Matrix & matrix) {

    if (matrix.isSquare()) {

        std::cout << "\n" << "Determinant of current matrix is " << matrix.determinant() << std::endl;
        os << "Determinant: " << matrix.determinant() << std::endl;

    } else {

        std::cout << "\n" << "Determinant could not be determined for current matrix." << std::endl;
        os << "Determinant could not be determined for current matrix." << std::endl;
    }

}

void writeSolution(std::ofstream & os, const Matrix & matrix) {

    if (not matrix.isSolvable()) {
        std::cout << "\n" << "Matrix could not be solved" << std::endl;
        os << "Matrix could not be solved." << std::endl;
    } else {

        try {
            const std::vector<double> solutions = solve(matrix);
            std::endl(std::cout);

            for (size_t i = 0; i < solutions.size(); ++i) {
                std::cout << "x" << (i + 1) << ": " << solutions[i] << std::endl;
                os << "x" << (i + 1) << ": " << solutions[i] << "<br>\n";
            }

        } catch (const no_solutions & e) {
            std::cout << "\n" << "Matrix has no solutions." << std::endl;
            os << "Matrix has no solutions." << std::endl;
        } catch (const infinite_solutions & e) {
            std::cout << "\n" << "Matrix has infinite number of solutions." << std::endl;
            os << "Matrix has infinite number of solutions." << std::endl;
        }

    }

}

void writeSection(std::ofstream & os, const std::string & sectionName, 
                  const std::function<void(std::ofstream & os, const Matrix&)> function, 
                  const Matrix & matrix) {

    std::cout << sectionName << ":" << std::endl;

    os << "<p>\n";
    os << "<h2>" << sectionName << "</h2>\n";

    function(os, matrix);

    os << "</p>\n";
    os << "<hr>" << std::endl;
    

}

void writeHeader(std::ofstream & os) {

    os << "<head>\n";

    os << "<title>Matrix</title>\n";

    os << "</head>\n";

}

void handleMatrix(const Matrix & matrix, const std::string & file) {

    std::ofstream os(file);

    if (not os) {
        throw std::runtime_error("Error: Output file \"" + file + "\" could not be opened.");
    }

    os << "<html>\n";

    writeHeader(os);

    Matrix transformed = gauss(matrix);

    os << "<body>\n";
    writeSection(os, "Input", writeMatrix, matrix);    
    writeSection(os, "Gaussian elimination", writeMatrix, transformed);
    writeSection(os, "Determinant", writeDeterminant, transformed);
    writeSection(os, "Solution", writeSolution, transformed);
    os << "</body>\n";

    os << "</html>";

    std::endl(std::cout);
    std::endl(os);

}

Matrix menuFile() {

    std::cout << "File: ";
    std::string filename;
    std::getline(std::cin, filename);

    return matrixFromFile(filename);

}

bool menuInputIsValid(const std::string & input) {

    if (input.size() != 1) {
        return false;
    }

    if (input[0] != '1' and input[0] != '2' and input[0] != '3') {
        return false;
    }

    return true;

}

std::string getOFile() {

    std::cout << "Output file: ";
    std::string file;
    std::cin >> file;

    return file;

}

void menu() {

    while (true) {
        
        std::cout << "  [1] Matrix from keyboard input [2] Matrix from file [3] Exit" << std::endl;

        std::string input;
        std::getline(std::cin, input);

        if (not menuInputIsValid(input)) {
            continue;
        }

        Matrix matrix;

        try {

            switch (input[0]) {
                case '1':
                    matrix = matrixFromStdin();
                    break;
                case '2':
                    matrix = menuFile();
                    break;
                case '3':
                    return;
                default:
                    break;
            }

        } catch (const std::exception & e) {
            std::cout << e.what() << std::endl;
            continue;
        }

        handleMatrix(matrix, getOFile());
    }
}

int fromArgs(const std::string & in, const std::string & out) {

    try {
        
        Matrix input = matrixFromFile(in);
        handleMatrix(input, out);

    } catch (const std::exception & e) {

        std::cout << e.what() << std::endl;
        return -1;

    }

    return 0;

}

int main(int argc, const char * argv[]) {

    if (argc == 3) {
        return fromArgs(argv[1], argv[2]);
    }

    menu();

}
