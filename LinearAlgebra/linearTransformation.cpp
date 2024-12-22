#include <iostream>
#include <vector>

using namespace std;

// Hàm đổi chỗ dòng của ma trận
void swapRows(vector<vector<double>>& A, vector<double>& b, int row1, int row2) {
    swap(A[row1], A[row2]);
    swap(b[row1], b[row2]);
}

// Hàm giải hệ phương trình tuyến tính Ax = b bằng phương pháp Gauss
vector<double> gaussianElimination(vector<vector<double>> A, vector<double> b) {
    int n = A.size();

    // Biến để lưu nghiệm
    vector<double> x(n);

    // Bước 1: Biến đổi ma trận A thành dạng tam giác trên
    for (int i = 0; i < n; ++i) {
        // Tìm dòng có phần tử lớn nhất ở cột i để tránh lỗi số học
        int maxRow = i;
        for (int j = i + 1; j < n; ++j) {
            if (abs(A[j][i]) > abs(A[maxRow][i])) {
                maxRow = j;
            }
        }

        // Đổi chỗ dòng i và dòng maxRow
        if (i != maxRow) {
            swapRows(A, b, i, maxRow);
        }

        // Biến đổi ma trận A thành dạng tam giác trên
        for (int j = i + 1; j < n; ++j) {
            double factor = A[j][i] / A[i][i];
            for (int k = i; k < n; ++k) {
                A[j][k] -= A[i][k] * factor;
            }
            b[j] -= b[i] * factor;
        }
    }

    // Bước 2: Giải hệ phương trình theo chiều ngược lại
    for (int i = n - 1; i >= 0; --i) {
        x[i] = b[i] / A[i][i];
        for (int j = i - 1; j >= 0; --j) {
            b[j] -= A[j][i] * x[i];
        }
    }

    return x;
}

int main() {
    // Ma trận A (3x3) và vector b
    vector<vector<double>> A = {
        {2, -1, 1},
        {3, 3, 9},
        {3, 3, 5}
    };
    vector<double> b = {8, 0, -2};

    // Giải hệ phương trình Ax = b
    vector<double> x = gaussianElimination(A, b);

    // In nghiệm
    cout << "Nghiem x = ";
    for (double xi : x) {
        cout << xi << " ";
    }
    cout << endl;

    return 0;
}