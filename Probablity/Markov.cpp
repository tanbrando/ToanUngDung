#include <iostream>
#include <vector>
#include <iomanip>

using namespace std;

// Hàm nhân ma trận (mat1 * mat2)
vector<vector<double>> multiplyMatrix(const vector<vector<double>>& mat1, const vector<vector<double>>& mat2) {
    int n = mat1.size();
    vector<vector<double>> result(n, vector<double>(n, 0.0));

    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            for (int k = 0; k < n; k++) {
                result[i][j] += mat1[i][k] * mat2[k][j];
            }
        }
    }
    return result;
}

// Hàm nhân vector và ma trận (vec * mat)
vector<double> multiplyVectorMatrix(const vector<double>& vec, const vector<vector<double>>& mat) {
    int n = vec.size();
    vector<double> result(n, 0.0);

    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            result[i] += vec[j] * mat[j][i];
        }
    }
    return result;
}

// Hàm in vector
void printVector(const vector<double>& vec) {
    for (double val : vec) {
        cout << fixed << setprecision(4) << val << " ";
    }
    cout << endl;
}

int main() {
    // Ma trận chuyển đổi trạng thái (P)
    vector<vector<double>> transitionMatrix = {
        {0.09, 0.31, 0.38, 0.22},
        {0.18, 0.16, 0.44, 0.22},
        {0.17, 0.26, 0.23, 0.34},
        {0.23, 0.37, 0.24, 0.16}
    };

    // Trạng thái ban đầu (π0)
    vector<double> initialState = {0.0, 1.0, 0.0, 0.0};

    // Số bước chuyển đổi (n)
    int steps = 4;

    cout << "Trạng thái ban đầu (π0): ";
    printVector(initialState);

    vector<double> currentState = initialState;

    // Lặp qua n bước
    for (int step = 1; step <= steps; step++) {
        currentState = multiplyVectorMatrix(currentState, transitionMatrix);
        cout << "Buoc" << step << ": ";
        printVector(currentState);
    }

    return 0;
}