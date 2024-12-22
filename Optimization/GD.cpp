#include <iostream>
#include <vector>
#include <cmath>
#include <random>
#include <Eigen/Dense>

using namespace Eigen;
using namespace std;

// Hàm gradient tính đạo hàm của hàm mất mát
VectorXd grad(const MatrixXd& Xbar, const VectorXd& y, const VectorXd& w) {
    int N = Xbar.rows();
    return (1.0 / N) * Xbar.transpose() * (Xbar * w - y);
}

double l(const MatrixXd& Xbar, const VectorXd& y, const VectorXd& w) {
    int N = Xbar.rows();
    return 0.5 * (1.0 / N) * (Xbar * w - y).squaredNorm();
}

// Thuật toán Gradient Descent
pair<vector<VectorXd>, int> myGradientDescent(const MatrixXd& Xbar, const VectorXd& y, VectorXd w_init, double alpha, int loop = 1000, double epsilon = 1e-4) {
    vector<VectorXd> w_history;
    w_history.push_back(w_init);

    for (int i = 0; i < loop; ++i) {
        VectorXd w_new = w_history.back() - alpha * grad(Xbar, y, w_history.back());

        if (grad(Xbar, y, w_new).norm() / w_new.size() < epsilon) {
            return {w_history, i};
        }

        w_history.push_back(w_new);
    }

    return {w_history, loop};
}

int main() {
    // Khởi tạo dữ liệu
    int N = 1000;
    double noise_std = 0.2;

    // Dữ liệu X (1000 điểm)
    VectorXd X = VectorXd(N);
    VectorXd y = VectorXd(N);
    random_device rd;
    mt19937 gen(rd());
    normal_distribution<> dist(0, noise_std);

    for (int i = 0; i < N; ++i) {
        X(i) = static_cast<double>(i) / N;
        y(i) = 4 + 3 * X(i) + dist(gen);  // Mối quan hệ tuyến tính với nhiễu
    }

    // Xbar = [1, X] (chèn thêm cột 1 cho bias)
    MatrixXd Xbar(N, 2);
    Xbar.col(0) = VectorXd::Ones(N);  // Cột 1 chứa toàn bộ giá trị 1
    Xbar.col(1) = X;                  // Cột 2 là X

    // Khởi tạo w
    VectorXd w_init(2);
    w_init << 2, 1;

    // Chạy Gradient Descent
    auto [w1, it1] = myGradientDescent(Xbar, y, w_init, 0.1);
    cout << "Phương pháp Gradient Descent: w = " << w1.back().transpose() << ", after " << it1 + 1 << " iterations, l = " << l(Xbar, y, w1.back()) << endl;

    // Phương pháp hồi quy tuyến tính chuẩn (nghịch đảo giả)
    MatrixXd A = Xbar.transpose() * Xbar;
    VectorXd b = Xbar.transpose() * y;
    VectorXd w_lr = A.colPivHouseholderQr().solve(b);
    cout << "Phương pháp Nghịch đảo: w = " << w_lr.transpose() << endl;

    // Vẽ kết quả
    vector<double> x_vals(N);
    vector<double> y_vals(N);
    vector<double> x_fit = {0, 1};
    vector<double> y_fit = {w_lr(0) + w_lr(1) * x_fit[0], w_lr(0) + w_lr(1) * x_fit[1]};

    for (int i = 0; i < N; ++i) {
        x_vals[i] = X(i);
        y_vals[i] = y(i);
    }

    return 0;
}