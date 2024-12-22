#include <iostream>
#include <vector>
#include <iomanip>
#include <algorithm>
#include <limits>

using namespace std;

// In một bảng simplex
void printTable(const vector<vector<double>>& table) {
    for (const auto& row : table) {
        for (double val : row) {
            cout << setw(10) << val << " ";
        }
        cout << endl;
    }
    cout << endl;
}

// Tính chỉ số lãi biên (relative profit)
vector<double> calculateRelativeProfit(const vector<vector<double>>& table, const vector<double>& c, int m, int n) {
    vector<double> relProf(n, 0);
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < m; ++j) {
            relProf[i] += table[j][1] * table[j][3 + i];
        }
        relProf[i] -= c[i];
    }
    return relProf;
}

// Phương pháp Simplex
vector<vector<double>> simplex(vector<double> c, vector<vector<double>> A, vector<double> b, vector<int> B = {}, vector<double> cb = {}) {
    int m = A.size(), n = A[0].size();

    // Nếu không có B và cb, khởi tạo
    if (B.empty()) {
        B.resize(m);
        for (int i = 0; i < m; ++i) B[i] = n - m + i;
    }
    if (cb.empty()) {
        cb.resize(m);
        for (int i = 0; i < m; ++i){
            cb[i] = c[c.size()-m + i];
        }  
    }

    // Khởi tạo bảng simplex
    vector<vector<double>> table(m, vector<double>(3 + n));
    for (int i = 0; i < m; ++i) {
        table[i][0] = B[i];
        table[i][1] = cb[i];
        table[i][2] = b[i];
        for (int j = 0; j < n; ++j) {
            table[i][3 + j] = A[i][j];
        }
    }


    bool reached = false;
    while (!reached) {
        vector<double> relProf = calculateRelativeProfit(table, c, m, n);
        // Kiểm tra nếu đã tối ưu
        bool optimal = true; 

        for (double p : relProf) {
            if (p > 0) {
                optimal = false; 
                break; 
            }
        }
        if (optimal) {
            reached = true;
            break;
        }

        // Tìm cột vào (chỉ số biến vào)
        int k = max_element(relProf.begin(), relProf.end()) - relProf.begin();

        // Kiểm tra không bị chặn (unbounded)
        double minRatio = 99999;
        int r = -1;
        for (int i = 0; i < m; ++i) {
            if (table[i][2] > 0 && table[i][3 + k] > 0) {
                double ratio = table[i][2] / table[i][3 + k];
                if (ratio < minRatio) {
                    minRatio = ratio;
                    r = i;
                }
            }
        }

        if (r == -1) {
            throw runtime_error("Unbounded solution.");
        }

        // Cập nhật bảng simplex
        double pivot = table[r][3 + k];
        for (int j = 2; j < 3 + n; ++j) {
            table[r][j] /= pivot;
        }

        for (int i = 0; i < m; ++i) {
            if (i != r) {
                double factor = table[i][3 + k];
                for (int j = 2; j < 3 + n; ++j) {
                    table[i][j] -= factor * table[r][j];
                }
            }
        }

        // Cập nhật biến cơ sở
        table[r][0] = k;
        table[r][1] = c[k];
    }

    return table;
}

// Phương pháp Two-Phase Simplex
pair<vector<double>, double> bigM(vector<double> c, vector<vector<double>> A, vector<double> b) {
    int m = A.size(), n = A[0].size();

    // Khởi tạo bài toán giai đoạn 1
    vector<double> c2(n + m, 0);
    for (int i = 0; i < m; ++i) c2[n + i] = 1;

    vector<vector<double>> table = simplex(c2, A, b);

    if (table.empty()) {
        return {};
    }

    vector<double> sol(n, 0);
    for (int i = 0; i < m; ++i) {
        sol[table[i][0]] = table[i][2];
    }

    double rs = 0;
    for (int i = 0; i < m; ++i) {
        rs += table[i][1] * table[i][2];
    }

    return {sol, rs};
}

int main() {
    double M = 1e6;
    vector<double> c = {-1, -2, 1, 0, 0, M, M};
    vector<vector<double>> A = {{-1, 4, -2, 1, 0, 0, 0},
                  {1, 1, 2, 0, -1, 1, 0},
                  {2, -1, 2, 0, 0, 0, 1}};

    vector<double> b = {6, 6,4};

    try {
        auto [sol, rs] = bigM(c, A, b);
        cout << "The optimal solution is: " << rs << endl;
        cout << "The objective function value is: ";
        for (double x : sol) cout << x << " ";
        cout << endl;
    } catch (const exception& e) {
        cout << e.what() << endl;
    }

    return 0;
}