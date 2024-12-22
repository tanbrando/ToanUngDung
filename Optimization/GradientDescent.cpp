#include <iostream>
#include <cmath>

using namespace std;

// Hàm mục tiêu: f(x) = x^2 - 4x + 4
double f(double x) {
    return 3 * exp(pow(x, 5) - pow(x, 4)) + pow(x, 2) - 20 * x + log(x + 25) - 10;
}

// Đạo hàm của hàm mục tiêu: f'(x) = 2x - 4
double f_prime(double x) {
    return 3 * (5 * pow(x, 4) - 4 * pow(x, 3)) * exp(pow(x, 5) - pow(x, 4)) + 2 * x - 20 + 1 / (x + 25); 
}

// Thuật toán Gradient Descent
double gradientDescent(double initial_x, double learning_rate, int max_iterations, double esp) {
    double x = initial_x;

    // Lặp qua các vòng lặp để cập nhật x
    for (int i = 0; i < max_iterations; ++i) {
        double gradient = f_prime(x);  // Tính gradient tại x
        x = x - learning_rate * gradient;  // Cập nhật giá trị x theo gradient
        if (fabs(initial_x - x) < esp) 
            break; 
        else 
            initial_x = x; 
    }

    return x;
}

int main() {
    double initial_x = 0.0;  // Giá trị khởi tạo của x
    double learning_rate = 0.01;  // Tốc độ học (learning rate)
    int max_iterations = 10000;  // Số vòng lặp tối đa
    double esp = 0.00001;  // Mức epsilon để dừng cập nhật x (thông số tối thiểu giá trị đ�� sai số)

    // Thực hiện Gradient Descent
    double optimized_x = gradientDescent(initial_x, learning_rate, max_iterations,esp);

    cout << "Optimized x = " << optimized_x << endl;
    cout << "f(optimized x) = " << f(optimized_x) << endl;

    return 0;
}