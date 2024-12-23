#include <iostream>
#include <vector>
#include <cmath>
#include <Eigen/Dense>

using namespace std;
using Matrix = vector<vector<double>>;

double epsilon = 1e-5;

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

void NormalizeVector(vector<double>& A) { 
    double sum = 0; 
    int size = A.size(); 
    for (int i = 0; i < size; i++) { 
        if (abs(A[i]) < epsilon) { 
            A[i] = 0; 
        } 
        sum += A[i] * A[i]; 
    } 
    double length = sqrt(sum); 
    for (int i = 0; i < size; i++) { 
        A[i] = A[i] / length; 
    } 
} 

void FindUi_WithEigenValueIsZero(vector<vector<double>>& U, int i, int m) { 
    // Tạo vector ngẫu nhiên
    vector<double> ei(m, 0);
    for (int k = 0; k < m; ++k) {
        ei[k] = static_cast<double>(rand()) / RAND_MAX; // Random values
    }

    vector<double> res(m, 0); 
    NormalizeVector(ei); 

    // Thực hiện trực giao hóa với các cột đã có
    for (int k = 0; k < i; k++) { 
        // Tính dot product
        double dot = 0;
        for (int j = 0; j < m; ++j) {
            dot += U[j][k] * ei[j];  // U.col(k) dot ei
        }

        // Trừ projection từ ei
        for (int j = 0; j < m; ++j) {
            res[j] += dot * U[j][k];  // res = res + dot * U.col(k)
        }
    }

    // Subtract the result from ei to make it orthogonal
    for (int j = 0; j < m; ++j) {
        ei[j] -= res[j];
    }

    // Chuẩn hóa lại vector
    NormalizeVector(ei); 

    // Lưu cột vào ma trận U
    for (int j = 0; j < m; ++j) {
        U[j][i] = ei[j];  // U.col(i) = ei
    }
}

Matrix calculateU(const Matrix &A, const Matrix &V, const Matrix &Sigma) {
    int m = A.size();
    int n = A[0].size();
    int r = min(Sigma.size(), Sigma[0].size());
    Matrix U(m, vector<double>(m, 0));

    int validColumns = 0;  // Đếm số lượng các cột hợp lệ đã được chuẩn hóa
    vector<int> validColumnOrders;

    for (int i = 0; i < r; ++i) {
        double sigma = Sigma[i][i];

        // Kiểm tra sigma gần bằng 0 và bỏ qua cột nếu cần
        if (fabs(sigma) < epsilon) continue; 

        // Lấy cột thứ i của V (V[:, i])
        vector<double> vi(n);
        for (int j = 0; j < n; ++j)
            vi[j] = V[j][i];

        // Tính toán A * V[:, i]
        vector<double> ui(m, 0);
        for (int j = 0; j < m; ++j)
            for (int k = 0; k < n; ++k)
                ui[j] += A[j][k] * vi[k];

        // Chuẩn hóa
        double norm = 0;
        for (int j = 0; j < m; ++j) {
            ui[j] /= sigma;
            norm += ui[j] * ui[j];
        }
        norm = sqrt(norm);

        // Kiểm tra xem norm có đủ lớn không
        if (norm < epsilon) {
            fill(ui.begin(), ui.end(), 1.0); // Hoặc chọn giá trị khác thay thế
            norm = 1.0;
        }

        // Lưu cột đã chuẩn hóa vào ma trận U
        for (int j = 0; j < m; ++j)
            U[j][validColumns] = ui[j] / norm;

        // Tăng validColumns sau mỗi lần chuẩn hóa thành công
        validColumnOrders.push_back(i);
        validColumns++;
    }

    // Hoàn thiện các cột còn thiếu bằng trực giao hóa
    for (int i = validColumns; i < m; ++i) {
        FindUi_WithEigenValueIsZero(U,i,m);
    }

    Matrix newU(m, vector<double>(m, 0));
    int count = 0;
    for (int i = 0; i < m; ++i) {
        if (count < validColumns && i == validColumnOrders[count]) {
            for (int j = 0; j < m; ++j)
                newU[j][i] = U[j][count];
            count++;
        } else {
            for (int j = 0; j < m; ++j)
                newU[j][i] = U[j][validColumns + count];
        }
    }


    return newU;
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