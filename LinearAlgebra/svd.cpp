#include <iostream>
#include <vector>
#include <cmath>
#include <Eigen/Dense>

using namespace std;
using Matrix = vector<vector<double>>;

Matrix transposeMatrix(const Matrix &A) {
    int m = A.size(), n = A[0].size();
    Matrix AT(n, vector<double>(m));
    for (int i = 0; i < m; ++i)
        for (int j = 0; j < n; ++j)
            AT[j][i] = A[i][j];
    return AT;
}

Matrix multiplyMatrix(const Matrix &A, const Matrix &B) {
    int m = A.size(), p = B.size(), n = B[0].size();
    Matrix result(m, vector<double>(n, 0));
    for (int i = 0; i < m; ++i)
        for (int j = 0; j < n; ++j)
            for (int k = 0; k < p; ++k)
                result[i][j] += A[i][k] * B[k][j];
    return result;
}

vector<double> eigenValues(const Matrix &A) {
    int n = A.size();
    Eigen::MatrixXd mat(n, n);
    for (int i = 0; i < n; ++i)
        for (int j = 0; j < n; ++j)
            mat(i, j) = A[i][j];

    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> solver(mat);
    Eigen::VectorXd eigenvalues = solver.eigenvalues();

    vector<double> result(n);
    for (int i = 0; i < n; ++i)
        result[i] = eigenvalues(i);
    return result;
}

Matrix eigenVectors(const Matrix &A, const vector<double> &eigenValues) {
    int n = A.size();
    Eigen::MatrixXd mat(n, n);
    for (int i = 0; i < n; ++i)
        for (int j = 0; j < n; ++j)
            mat(i, j) = A[i][j];

    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> solver(mat);
    Eigen::MatrixXd eigenvectors = solver.eigenvectors();

    Matrix result(n, vector<double>(n));
    for (int i = 0; i < n; ++i)
        for (int j = 0; j < n; ++j)
            result[i][j] = eigenvectors(i, j);
    return result;
}

Matrix orthonormalize(const Matrix &V) {
    int m = V.size(), n = V[0].size();
    Matrix result(m, vector<double>(n, 0));
    for (int i = 0; i < n; ++i) {
        vector<double> vi = V[i];
        double norm = 0;
        for (double x : vi)
            norm += x * x;
        norm = sqrt(norm);
        for (int j = 0; j < m; ++j)
            result[j][i] = vi[j] / norm;
    }
    return result;
}

Matrix calculateSigma(const vector<double> &eigenValues, int m, int n) {
    int r = min(m, n);
    Matrix Sigma(m, vector<double>(n, 0));
    for (int i = 0; i < r; ++i)
        Sigma[i][i] = sqrt(eigenValues[i]);
    return Sigma;
}

Matrix calculateU(const Matrix &A, const Matrix &V, const Matrix &Sigma) {
    int m = A.size(), n = V[0].size();
    Matrix U(m, vector<double>(n, 0));

    for (int i = 0; i < n; ++i) {
        vector<double> vi(n);
        for (int j = 0; j < n; ++j)
            vi[j] = V[j][i];

        vector<double> ui(m, 0);
        for (int j = 0; j < m; ++j)
            for (int k = 0; k < n; ++k)
                ui[j] += A[j][k] * vi[k];
        
        double norm = 0;
        for (double x : ui)
            norm += x * x;
        norm = sqrt(norm);
        for (int j = 0; j < m; ++j)
            U[j][i] = ui[j] / norm;
    }
    return U;
}

void roundMatrix(Matrix &A) {
    for (auto &row : A)
        for (auto &x : row)
            x = round(x * 1000) / 1000.0;
}

int main() {
    // Ma trận A đầu vào
    Matrix A = {
        {1, 2},
        {3, 4},
        {5, 6}
    };

    int m = A.size();
    int n = A[0].size();

    // Tính A^T và A^T A
    Matrix AT = transposeMatrix(A);
    Matrix ATA = multiplyMatrix(AT, A);

    // Tính trị riêng và vector riêng
    vector<double> EigenValues = eigenValues(ATA);
    Matrix V = eigenVectors(ATA, EigenValues);
    V = transposeMatrix(V);
    V = orthonormalize(V);

    // Tính ma trận Sigma
    Matrix Sigma = calculateSigma(EigenValues, m, n);

    // Tính ma trận U
    Matrix U = calculateU(A, V, Sigma);

    // Làm tròn ma trận
    roundMatrix(U);
    roundMatrix(Sigma);
    roundMatrix(V);

    // Xuất kết quả
    cout << "Ma trận U:\n";
    for (const auto &row : U) {
        for (double x : row)
            cout << x << " ";
        cout << "\n";
    }

    cout << "Ma trận Sigma:\n";
    for (const auto &row : Sigma) {
        for (double x : row)
            cout << x << " ";
        cout << "\n";
    }

    cout << "Ma trận V:\n";
    for (const auto &row : V) {
        for (double x : row)
            cout << x << " ";
        cout << "\n";
    }

    return 0;
}