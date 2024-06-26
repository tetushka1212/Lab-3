#pragma once
#include<iostream>
using namespace std;
class Matrix
{
private:
	double** matrix;
	size_t rank;
public:
	Matrix(size_t n)
	{
		this->rank = n;
		this->matrix = new double* [n];
		for (size_t i = 0; i < n; i++) {
			matrix[i] = new double[n];
			for (size_t j = 0; j < n; j++)
			{
				matrix[i][j] = rand() % 11;
			}
		}
	}
	double CalculateDeterminant();
	double Determinant(double** data, size_t n);
	void Print();
	double* CramerSolution(double* constants);
	~Matrix() {
		for (int i = 0; i < this->rank; i++) {
			delete[] this->matrix[i];
		}
		delete[] this->matrix;
	}
};

