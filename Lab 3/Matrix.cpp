#include "Matrix.h"

double Matrix::CalculateDeterminant()
{
    size_t n = this->rank;
    double det = 1.0;

    for (int i = 0; i < n; ++i) {
        int pivot_row = i;
        while (pivot_row < n && matrix[pivot_row][i] == 0) {
            pivot_row++;
        }

        if (pivot_row == n) {
            return 0.0; 
        }

        if (pivot_row != i) {
            swap(matrix[pivot_row], matrix[i]); 
            det *= -1; 
        }

        double pivot = matrix[i][i]; 

        det *= pivot; 

        for (int j = i + 1; j < n; ++j) {
            double factor = matrix[j][i] / pivot; 
            for (int k = i; k < n; ++k) {
                matrix[j][k] -= factor * matrix[i][k];
            }
        }
    }
}

double Matrix::Determinant(double** data, size_t n)
{
    double det = 1;

    for (size_t i = 0; i < n; i++) {
        size_t k = i;
        for (size_t j = i + 1; j < n; j++) {
            if (abs(data[j][i]) > abs(data[k][i])) {
                k = j;
            }
        }
        if (abs(data[k][i]) < 1e-6) {
            return 0;
        }
        std::swap(data[i], data[k]);
        if (i != k) {
            det = -det;
        }

        det *= data[i][i];
        for (size_t j = i + 1; j < n; j++) {
            data[i][j] /= data[i][i];
        }
        for (size_t j = 0; j < n; j++) {
            if (j != i && abs(data[j][i]) > 1e-6) {
                for (size_t k = i + 1; k < n; k++) {
                    data[j][k] -= data[i][k] * data[j][i];
                }
            }
        }
    }

    return det;
}

void Matrix::Print()
{
    for (size_t i = 0; i < rank; ++i)
    {
        for (size_t j = 0; j < rank; ++j) {
            cout << matrix[i][j] << " ";
    }
    cout << endl;
}
}

double* Matrix::CramerSolution(double* constants)
{
    double det_A = this->CalculateDeterminant();
    if (det_A == 0) {
        std::cout << "Determinant of the matrix is zero, cannot solve the system of equations by Cramer's rule." << std::endl;
        exit(EXIT_FAILURE);
    }

    double* solutions = new double[this->rank];
    for (int i = 0; i < this->rank; i++) {
        double** matrix_copy = new double* [this->rank];
        for (int j = 0; j < this->rank; j++) {
            matrix_copy[j] = new double[this->rank];
            for (int k = 0; k < this->rank; k++) {
                matrix_copy[j][k] = this->matrix[j][k];
            }
            matrix_copy[j][i] = constants[j];
        }
        double det_Ai = Determinant(matrix_copy, this->rank);
        solutions[i] = det_Ai / det_A;
        for (int j = 0; j < this->rank; j++) {
            delete[] matrix_copy[j];
        }
        delete[] matrix_copy;
    }

    return solutions;
}
