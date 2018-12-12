/**
* @file matrix.cpp
* @brief This file includes all classes and functions required to 
*  perform gaussian elimination on matrices, calculate their determinant
*  and solve them
*
* @author Filip Peterek
* @date 12/12/2018
*/

#include <iostream>
#include <vector>
#include <fstream>
#include <sstream>
#include <utility>
#include <functional>


/** 
* @brief Exception thrown when invalid input is received
*/
class invalid_input : public std::runtime_error {

public:

    invalid_input(const std::string & input) 
        : std::runtime_error("Invalid input: " + input) { }

};


/** 
* @brief Exception thrown if matrix cannot be solved because it has
*  an infinite amount of solutions
*/
class infinite_solutions : public std::runtime_error {

public:

    infinite_solutions()
        : std::runtime_error("Error: Matrix could not be solved because it has infinite number of solutions") { }

};

/** 
* @brief Exception thrown if matrix cannot be solved because it has
*  no solutions
*/
class no_solutions : public std::runtime_error {

public:

    no_solutions()
        : std::runtime_error("Error: Matrix has no solutions.") { }
};

/** 
* @brief Raw matrix collection
*/
typedef std::vector<std::vector<double>> matrix_t;

/** 
* @brief Matrix class. Holds raw matrix using matrix_t. Provides additional
*  information about matrix.
*/
class Matrix {

    matrix_t matrix; /**< Raw matrix */

public:

    /** @brief Default constructor */
    Matrix() { }
    /** @brief Copy constructor */
    Matrix(const matrix_t & matrix) : matrix(matrix) { }

    /** 
    * @brief Appends a line to the end of matrix if line size matches the rest of the matrix 
    * @param line Line to be appended
    */
    void appendLine(const std::vector<double> & line) {

        if (matrix.size() and matrix[0].size() != line.size()) {
            throw std::runtime_error(std::string("Error: Row too ") + 
                                    (line.size() > matrix[0].size() ? "long" : "short"));
        }

        matrix.emplace_back(line);

    }

    /** @brief Get underlying raw matrix */
    const matrix_t & getMatrix() const {
        return matrix;
    }

    /** @brief Get number of rows in matrix */
    size_t rows() const {
        return matrix.size();
    }

    /** @brief Get number of columns in matrix */
    size_t columns() const {
        
        if (not rows()) {
            return 0;
        }

        return matrix.front().size();

    }

    /** 
    * @brief Get a single line from matrix 
    * @param index Line index 
    * @return Reference to a single line
    */
    std::vector<double> & getLine(const size_t index) {
        return matrix[index];
    }

    /** 
    * @brief Get a single element from matrix 
    * @param line Line index
    * @param index Element index 
    * @return Reference to a single element
    */
    double & get(const size_t line, const size_t index) {
        return matrix[line][index];
    }

    /** 
    * @brief Get a single line from matrix 
    * @param index Line index 
    * @return Const reference to a single line
    */
    const std::vector<double> & getLine(const size_t index) const {
        return matrix[index];
    }

    /** 
    * @brief Get a single element from matrix 
    * @param line Line index
    * @param index Element index 
    * @return Copy of a single element
    */
    double get(const size_t line, const size_t index) const {
        return matrix[line][index];
    }

    /** @brief Checks whether matrix is a square (number of rows == number of columns) */
    bool isSquare() const {

        if (not rows()) {
            return false;
        }

        return rows() == columns();

    }

    /** @brief Calculates the determinant of a square matrix after 
    *    gaussian elimination has been performed */
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

    /** @brief Checks whether matrix can be solved */
    bool isSolvable() const {
        return rows() + 1 == columns();
    }

};

/**
* @brief Swaps two rows of a matrix
* @param matrix Raw matrix
* @param r1 Index of first row to be swapped
* @param r2 Index of second row to be swapped
*/
void swapTwo(matrix_t & matrix, const size_t r1, const size_t r2) {

    const std::vector<double> v1 = matrix[r1];
    const std::vector<double> v2 = matrix[r2];

    matrix[r1] = v2;
    matrix[r2] = v1;

}

/**
* @brief Attempts to swap rows in matrix if pivot is zero
* @param matrix Raw matrix
* @param pivot Pivot index
*/
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

/**
* @brief Performs subtraction between lines to reduce all elements underneath pivot to zero
* @param matrix Raw matrix
* @param pivot Pivot index
*/
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

/**
* @brief Checks whether line consists of only zeroes
* @param line Line to be checked
*/
bool isZeroOnly(const std::vector<double> & line) {

    for (double d : line) {

        if (d) {
            return false;
        }

    }

    return true;

}

/**
* @brief Removes zero only lines from matrix
* @param matrix Matrix to have zero only lines removed
* @return Matrix of type matrix_t without any zero only lines
*/
matrix_t removeZeroOnly(const matrix_t & matrix) {

    matrix_t res;

    for (const std::vector<double> & line : matrix) {

        if (not isZeroOnly(line)) {
            res.emplace_back(line);
        }

    }

    return res;

}

/**
* @brief Performs gaussian elimination on a matrix
* @param matrix Matrix to be transformed
* @return New matrix, which has been transformed using gaussian elimination
*/
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

/**
* @brief Reallocates an array of doubles
* @param arr Old array
* @param oldSize Current size of array
* @param newSize Desired size
* @return Reallocated array
*/
double * reallocate(double * arr, const size_t oldSize, const size_t newSize) {

    double * newArr = new double[newSize];

    for (size_t i = 0; i < oldSize; ++i) {
        newArr[i] = arr[i];
    }

    delete[] arr;

    return newArr;

}

/**
* @brief Solutions are stored in a reverse order. This function orders solutions
*  properly and stores them in an std::vector rather than a C array
* @param solved C array of solutions
* @param sSize Size of solved
* @return std::vector with it's elements in proper order
*/
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

/**
* @brief Solves a single line and finds a solution to an equation using solutions
*  to previous equations
* @param line Line to be solved
* @param solved C array of previous solutions
* @param sSize Size of solved
* @return Solution to current line
*/
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

/**
* @brief Solves a matrix
* @param matrix Matrix to be solved
* @return Vector of solutions
*/
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

/**
* @brief Overload for std::ostream::operator<< which allows writing matrices
*  right into ostreams
*/
std::ostream & operator<<(std::ostream & os, const Matrix & matrix) {

    for (const auto vector : matrix.getMatrix()) {
        for (const double i : vector) {
            os << i << " ";
        }
        std::endl(os);
    }

    return os;

}

/**
* @brief Checks if line is a valid matrix definition. Throws invalid_input
*  if not.
* @param line Line to be checked
*/
void checkLineIsValid(const std::string & line) {

    bool dashAllowed = true;

    for (const char c : line) {

        if ((c == '-' and dashAllowed) or isdigit(c))  {
            dashAllowed = false;
            continue;
        }

        if (isspace(c)) {
            dashAllowed = true;
            continue;
        }

        throw invalid_input(line);
    }

}

/**
* @brief Parses a line containing a matrix row definition
* @param line Line to be checked
* @return An std::vector of doubles which represents a single line
*/
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


/**
* @brief Checks whether string consists only of whitespace characters
* @param str String to be checked
*/
bool isWSOnly(const std::string & str) {

    for (const char c : str) {
        if (not isspace(c)) {
            return false;
        }
    }

    return true;

}

/**
* @brief Reads a matrix object from selected stream
* @param is Input stream from which the matrix is read
*/
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

/**
* @brief Reads a matrix from stdin
*/
Matrix matrixFromStdin() {

    std::cout << "Input matrix (empty line denotes end of matrix): " << std::endl;

    Matrix matrix = readMatrix(std::cin);

    return matrix;

}

/**
* @brief Reads a matrix from a file
* @param filename File containing a matrix 
*/
Matrix matrixFromFile(const std::string & filename) {

    std::ifstream input(filename);

    if (not input) {
        throw std::runtime_error("Error opening file: " + filename);
    }

    Matrix matrix = readMatrix(input);

    return matrix;
    

}

/**
* @brief Writes a matrix both to stdout and a specified ostream as an html table
* @param os Output stream to an html file
* @param matrix Matrix to be written 
*/
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

/**
* @brief Writes determinant of a matrix to both stdout and an ostream as an html element
* @param os Output stream to an html file
* @param matrix Matrix after gaussian elimination application
*/
void writeDeterminant(std::ofstream & os, const Matrix & matrix) {

    if (matrix.isSquare()) {

        std::cout << "\n" << "Determinant of current matrix is " << matrix.determinant() << std::endl;
        os << "Determinant: " << matrix.determinant() << std::endl;

    } else {

        std::cout << "\n" << "Determinant could not be determined for current matrix." << std::endl;
        os << "Determinant could not be determined for current matrix." << std::endl;
    }

}

/**
* @brief Writes the solution to a matrix to both stdout and an ostream as an html element
* @param os Output stream to an html file
* @param matrix Matrix after gaussian elimination application
*/
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

/**
* @brief Writes a section to both stdout and an ostream as an html element
* @param os Output stream to an html file
* @param sectionName Name of section, used as a heading in html file
* @param function Function which handles desired functionality (eg. determinant, solutions, etc...)
* @param matrix Matrix object passed to function
*/
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

/**
* @brief Writes an html header to an output stream
* @param os Output stream to an html file
*/
void writeHeader(std::ofstream & os) {

    os << "<head>\n";
    os << "<title>Matrix</title>\n";
    os << "</head>\n";

}

/**
* @brief Performs calculations on selected matrix and writes the results to both a file and stdout
* @param matrix Matrix on which calculations should be performed
* @param file Output file
*/
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

/**
* @brief Asks for a file containing a matrix definition and returns the matrix as a Matrix object
* @return Matrix object parsed from file
*/
Matrix menuFile() {

    std::cout << "File: ";
    std::string filename;
    std::getline(std::cin, filename);

    return matrixFromFile(filename);

}

/**
* @brief Checks if input is a valid menu option
* @param input Input from user
*/
bool menuInputIsValid(const std::string & input) {

    if (input.size() != 1) {
        return false;
    }

    if (input[0] != '1' and input[0] != '2' and input[0] != '3') {
        return false;
    }

    return true;

}

/**
* @brief Asks for name of output file
* @return Returs name of output file specified by user
*/
std::string getOFile() {

    std::cout << "Output file: ";
    std::string file;
    std::getline(std::cin, file);

    return file;

}

/**
* @brief Displays a menu and handles user input
*/
void menu() {

    while (true) {
        
        std::cout << "[1] Matrix from keyboard input [2] Matrix from file [3] Exit" << std::endl;
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

/**
* @brief Reads matrix from a file and writes it to html file specified using command line arguments
* @param in Input file passed as a command line argument
* @param out Output file passed as a command line argument
* @return Returns -1 on failure and 0 on success
*/
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

/**
* @brief Calls fromArgs() if exactly two command line arguments were passed to the programs,
*  otherwise calls menu()
*/
int main(int argc, const char * argv[]) {

    if (argc == 3) {
        return fromArgs(argv[1], argv[2]);
    }

    menu();

}


