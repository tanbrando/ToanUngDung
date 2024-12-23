#include <iostream>
#include <vector>
#include <cmath>
#include <Eigen/Dense>

using namespace std;
using Matrix = vector<vector<double>>;

int epsilon = 1e-5;

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
    for (int i = 0; i < r; ++i) {
        if (eigenValues[i] > epsilon) { // Bỏ qua trị riêng nhỏ hơn ngưỡng
            Sigma[i][i] = sqrt(eigenValues[i]);
        } else {
            Sigma[i][i] = 0; // Đặt giá trị bằng 0 nếu trị riêng nhỏ
        }
    }
    return Sigma;
}

Matrix calculateU(const Matrix &A, const Matrix &V, const Matrix &Sigma) {
    int m = A.size();        // Số dòng của A
    int n = A[0].size();     // Số cột của A
    int r = Sigma.size();    // Số giá trị kỳ dị (rank của A)

    Matrix U(m, vector<double>(m, 0));  // Kích thước U: m x m

    int validColumns = 0; // Đếm số cột hợp lệ trong U
    for (int i = 0; i < r; ++i) {
        // Lấy giá trị kỳ dị sigma_i
        double sigma = Sigma[i][i];
        if (sigma < 1e-6) { // Bỏ qua trường hợp sigma bằng 0
            cerr << "Cảnh báo: Giá trị kỳ dị sigma[" << i << "] = 0, bỏ qua cột tương ứng." << endl;
            continue;
        }

        // Lấy cột thứ i của V (V[:, i])
        vector<double> vi(n);
        for (int j = 0; j < n; ++j)
            vi[j] = V[j][i];

        // Tính toán A * V[:, i]
        vector<double> ui(m, 0);
        for (int j = 0; j < m; ++j)
            for (int k = 0; k < n; ++k)
                ui[j] += A[j][k] * vi[k];

        // Chuẩn hóa bằng cách chia cho sigma_i
        for (int j = 0; j < m; ++j)
            U[j][validColumns] = ui[j] / sigma;

        validColumns++;
    }

    // Hoàn thiện các cột còn thiếu bằng vector trực giao
    for (int i = validColumns; i < m; ++i) {
        vector<double> orthogonal(m, 0);
        orthogonal[i] = 1.0; // Vector đơn vị
        for (int j = 0; j < validColumns; ++j) {
            vector<double> uj(m);
            for (int k = 0; k < m; ++k)
                uj[k] = U[k][j];

            // Loại bỏ thành phần trùng tuyến
            double dot = 0;
            for (int k = 0; k < m; ++k)
                dot += orthogonal[k] * uj[k];

            for (int k = 0; k < m; ++k)
                orthogonal[k] -= dot * uj[k];
        }

        // Chuẩn hóa vector trực giao
        double norm = 0;
        for (double x : orthogonal)
            norm += x * x;
        norm = sqrt(norm);

        for (int j = 0; j < m; ++j)
            U[j][i] = orthogonal[j] / norm;
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
        {2, 4, 8},
        {4, 8, 16},
        {1, 2, 3}
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
    cout << "Ma tran U:\n";
    for (const auto &row : U) {
        for (double x : row)
            cout << x << " ";
        cout << "\n";
    }

    cout << "Ma tran Sigma:\n";
    for (const auto &row : Sigma) {
        for (double x : row)
            cout << x << " ";
        cout << "\n";
    }

    cout << "Ma tran V:\n";
    for (const auto &row : V) {
        for (double x : row)
            cout << x << " ";
        cout << "\n";
    }

    return 0;
}