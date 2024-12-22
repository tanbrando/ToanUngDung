#include <iostream>
#include <cmath>

using namespace std;

// Hàm mục tiêu: f(x) = x^2 - 4x + 4
double f(double x) {
    return (exp(2 * x) + 8 * pow(x, 2) + 4 * x) / (5 - x) - 10 * x; 
}

// Đạo hàm của hàm mục tiêu: f'(x) = 2x - 4
double f_prime(double x) {
    return ((2*exp(2*x)+16*x+4)*(5-x)+(exp(2*x)+8*pow(x,2)+4*x))/pow(5-x,2)-10; 
}

// Thuật toán Gradient Descent với Momentum
double gradientDescentWithMomentum(double initial_x, double learning_rate, double momentum, int max_iterations, double epsilon) {
    double x = initial_x;
    double v = 0;  // Khởi tạo momentum (v)

    // Lặp qua các vòng lặp để cập nhật x
    for (int i = 0; i < max_iterations; ++i) {
        double gradient = f_prime(x);  // Tính gradient tại x
        v = momentum * v + learning_rate * gradient;  // Cập nhật momentum
        x = x - v;  // Cập nhật giá trị x theo momentum

        if (abs(f(x)) < epsilon) {
            break;
        }
    }

    return x;
}

int main() {
    double initial_x = 0.0;  // Giá trị khởi tạo của x
    double learning_rate = 0.001;  // Tốc độ học (learning rate)
    double momentum = 0.1;  // Hệ số momentum
    int max_iterations = 1000;  // Số vòng lặp tối đa
    double epsilon = 1e-5;  // Mức đ�� chấp nhận của sai số (ε)

    // Thực hiện Gradient Descent với Momentum
    double optimized_x = gradientDescentWithMomentum(initial_x, learning_rate, momentum, max_iterations,epsilon);

    cout << "Optimized x = " << optimized_x << endl;
    cout << "f(optimized x) = " << f(optimized_x) << endl;

    return 0;
}