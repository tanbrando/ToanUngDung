#include <iostream>
#include <vector>
#include <cmath>
#include <iomanip>

using namespace std;

typedef vector<vector<double>> Matrix;

void printMatrix(const Matrix& mat) {
    for (const auto& row : mat) {
        for (double val : row) {
            cout << setw(10) << val << " ";
        }
        cout << endl;
    }
    cout << endl;
}

Matrix inverseMatrix(vector<vector<double>>& matrix) {
    int n = matrix.size();
    Matrix augmentedMatrix(n, vector<double>(2 * n, 0));

    // Tạo ma trận kết hợp [matrix | I]
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            augmentedMatrix[i][j] = matrix[i][j];  // Sao chép ma trận ban đầu
        }
        augmentedMatrix[i][n + i] = 1;  // Tạo ma trận đơn vị
    }

    // Áp dụng khử Gauss-Jordan
    for (int i = 0; i < n; i++) {
        // Tìm phần tử chính và hoán đổi hàng nếu cần
        if (augmentedMatrix[i][i] == 0) {
            for (int j = i + 1; j < n; j++) {
                if (augmentedMatrix[j][i] != 0) {
                    swap(augmentedMatrix[i], augmentedMatrix[j]);
                    break;
                }
            }
        }

        // Chuẩn hóa hàng hiện tại
        double divisor = augmentedMatrix[i][i];
        if (divisor == 0) {
            cout << "Ma trận không có nghịch đảo!" << endl;
            return Matrix();  // Trả về ma trận rỗng
        }
        for (int j = 0; j < 2 * n; j++) {
            augmentedMatrix[i][j] /= divisor;
        }

        // Khử Gauss-Jordan cho các hàng khác
        for (int j = 0; j < n; j++) {
            if (j != i) {
                double factor = augmentedMatrix[j][i];
                for (int k = 0; k < 2 * n; k++) {
                    augmentedMatrix[j][k] -= factor * augmentedMatrix[i][k];
                }
            }
        }
    }

    // Tách ma trận nghịch đảo ra khỏi ma trận kết hợp
    Matrix inverse(n, vector<double>(n, 0));
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            inverse[i][j] = augmentedMatrix[i][n + j];
        }
    }

    return inverse;
}

Matrix transposeMatrix(vector<vector<double>>& matrix) {
    int rows = matrix.size();
    int cols = matrix[0].size();
    Matrix transpose(cols, vector<double>(rows, 0));

    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            transpose[j][i] = matrix[i][j];
        }
    }

    return transpose;
}
double myRound(double x) {
    return round(x * 1000000) / 1000000.0;
}

void roundMatrix(Matrix& mat) {
    for (auto& row : mat) {
        for (auto& val : row) {
            val = myRound(val);
        }
    }
}

Matrix identityMatrix(int n) {
    Matrix I(n, vector<double>(n, 0.0));
    for (int i = 0; i < n; ++i) {
        I[i][i] = 1.0;  // Thiết lập các phần tử đường chéo chính bằng 1
    }
    return I;
}

Matrix subtractMatrices(const Matrix& A, const Matrix& B) {
    int n = A.size();
    Matrix result(n, vector<double>(n, 0.0));
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            result[i][j] = A[i][j] - B[i][j];
        }
    }
    return result;
}
vector<double> subtractVector(const vector<double> &A, const vector<double> &B) {
    int n = A.size();
    vector<double> result(n, 0.0);
    for (int i = 0; i < n; ++i) {
        result[i] = A[i] - B[i];
    }
    return result;
}

double dotProduct(const vector<double>& v1, const vector<double>& v2) {
    int n = v1.size();
    double sum = 0;
    for (int i = 0; i < n; ++i) {
        sum += v1[i] * v2[i];
    }
    return sum;
}

vector<double> getVectorFromMatrix(const Matrix& mat,int n) {
    int dem = mat.size();
    vector<double> result;
    for (int i = 0; i < dem; ++i) {
        result.push_back(mat[i][n]);
    }
    return result;
}

Matrix multiplyMatrix(const vector<vector<double>>& A, const vector<vector<double>>& B) {
    int m = A.size();
    int n = B[0].size();
    int common = B.size();
    Matrix C(m, vector<double>(n, 0.0));
    for (int i = 0; i < m; ++i) {
        for (int j = 0; j < n; ++j) {
            for (int k = 0; k < common; ++k) {
                C[i][j] += A[i][k] * B[k][j];
            }
        }
    }
    return C;
}

// Nhân ma trận với vector
vector<double> multiplyMatrixVector(const Matrix& A, const vector<double>& v) {
    int n = A.size();
    vector<double> result(n, 0.0);
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            result[i] += A[i][j] * v[j];
        }
    }
    return result;
}

vector<double> multiplyScalarMatrixVector(const double scalar,const Matrix& A, const vector<double>& v) {
    int n = A.size();
    vector<double> result(n, 0.0);
    result = multiplyMatrixVector(A,v);
    for (int i = 0; i < n; ++i) result[i] *= scalar;
    return result;
}

vector<double> multiplyScalarVector(const double scalar, const vector<double>& v) {
    vector<double> result = v;
    for (int i = 0; i < result.size(); ++i) result[i] *= scalar;
    return result;
}

Matrix convertToDiagonalMatrix(const vector<double>& v) {
    Matrix result(v.size(), vector<double>(v.size(), 0.0));
    for (int i = 0; i < v.size(); ++i) {
        result[i][i] = v[i];
    }
    return result;
}
double vectorNorm(const vector<double>& v) {
    double sum = 0.0;
    for (double val : v) sum += val * val;
    return sqrt(sum);
}

pair<Matrix,Matrix> qrDecomposition(const Matrix& A) {
    int n = A.size();
    Matrix Q(n, vector<double>(n, 0.0));
    Matrix R(n, vector<double>(n, 0.0));

    Matrix a = A;

    for (int k = 0; k < n; k++) {
        vector<double> u = getVectorFromMatrix(a,k);
        for (int j = 0; j < k; j++) {
            vector<double> v = getVectorFromMatrix(Q, j);
            R[j][k] = dotProduct(u, v);
            for (int i = 0; i < n; i++) {
                u[i] -= R[j][k] * v[i];
            }
        }

        R[k][k] = vectorNorm(u);
        for (int i = 0; i < n; i++) {
            Q[i][k] = u[i] / R[k][k];
        }
    }

    return {Q, R};    
}

vector<double> qrAlgorithm(Matrix A, int maxIterations = 1000, double tolerance = 1e-9) {
    int n = A.size();
    vector<double> eigenvalues(n, 0.0);

    for (int iter = 0; iter < maxIterations; ++iter) {
        // Phân rã QR
        auto [Q, R] = qrDecomposition(A);

        // Cập nhật ma trận A mới
        A = multiplyMatrix(R, Q);

        // Kiểm tra điều kiện hội tụ
        bool converged = true;
        for (int i = 0; i < n - 1; ++i) {
            if (fabs(A[i + 1][i]) > tolerance) {
                converged = false;
                break;
            }
        }
        if (converged) break;
    }

    // Trị riêng nằm trên đường chéo chính
    for (int i = 0; i < n; ++i) {
        eigenvalues[i] = A[i][i];
    }

    return eigenvalues;
}

bool gaussElimination(Matrix& A, vector<double>& b) {
    int n = A.size();
    for (int i = 0; i < n; ++i) {
        // Tìm phần tử lớn nhất trong cột hiện tại
        double maxEl = abs(A[i][i]);
        int maxRow = i;
        for (int k = i + 1; k < n; k++) {
            if (abs(A[k][i]) > maxEl) {
                maxEl = abs(A[k][i]);
                maxRow = k;
            }
        }
        // Hoán đổi hàng nếu cần
        for (int k = i; k < n; k++) {
            swap(A[maxRow][k], A[i][k]);
        }
        swap(b[maxRow], b[i]);

        // Giảm bậc ma trận
        for (int k = i + 1; k < n; k++) {
            double coeff = -A[k][i] / A[i][i];
            for (int j = i; j < n; j++) {
                A[k][j] += coeff * A[i][j];
            }
            b[k] += coeff * b[i];
        }
    }

    roundMatrix(A);

    // Kiểm tra xem có giải pháp hay không
    for (int i = n - 1; i >= 0; i--) {
        if (A[i][i] == 0) {
            if (b[i] != 0) return false; // Không có giải pháp
        }
    }

    return true; // Có giải pháp
}


// Hàm tính vector riêng từ trị riêng
vector<double> eigenVector(const Matrix& A, double eigenvalue) {
    int n = A.size();
    Matrix modifiedMatrix = subtractMatrices(A, identityMatrix(n));
    for (int i = 0; i < n; ++i) {
        modifiedMatrix[i][i] -= eigenvalue - 1;
    }
    vector<double> b(n, 0.0);
    if (!gaussElimination(modifiedMatrix, b)) {
        throw runtime_error("Không có vector riêng tương ứng với trị riêng này.");
    }


    vector<double> eigenvector(n, 0.0);
    for (int i = 0; i < n; ++i) {
        if (modifiedMatrix[i][i] == 0) {
            eigenvector[i] = 1; 
        }
    }
    for (int i = n - 2; i >= 0; --i) {
        for (int j = i + 1; j < n; ++j) {
            eigenvector[i] -= modifiedMatrix[i][j] * eigenvector[j];
        }
        eigenvector[i] /= modifiedMatrix[i][i];
    }

    return eigenvector;
}
Matrix eigenVectors(const Matrix& A, const vector<double>& eigenvalues) {
    Matrix eigenvectors;
    for (double eigenvalue : eigenvalues) {
        try {
            vector<double> ev = eigenVector(A, eigenvalue);
            eigenvectors.push_back(ev);
        } catch (const runtime_error& e) {
            cout << e.what() << endl;
        }
    }
    return eigenvectors;
}

Matrix orthonormalize(const vector<vector<double>>& vectors) {
    int n = vectors.size();
    int m = vectors[0].size();
    Matrix orthonormalVectors(n, vector<double>(m, 0.0));

    for (int i = 0; i < m; ++i) {
        // Tính độ dài của vector
        double norm = 0.0;
        for (int j = 0; j < m; ++j) {
            norm += vectors[j][i] * vectors[j][i];
        }
        norm = sqrt(norm);

        // Chuẩn hóa vector
        for (int j = 0; j < n; ++j) {
            orthonormalVectors[j][i] = vectors[j][i] / norm;
        }

        // Đảm bảo các vector khác là trực giao (nếu cần thiết)
        for (int j = 0; j < i; ++j) {
            double dotProduct = 0.0;
            for (int k = 0; k < n; ++k) {
                dotProduct += orthonormalVectors[k][i] * orthonormalVectors[k][j];
            }

            // Trừ phần tử chồng chéo
            for (int k = 0; k < n; ++k) {
                orthonormalVectors[k][i] -= dotProduct * orthonormalVectors[k][j];
            }

            // Chuẩn hóa lại vector sau khi loại bỏ chồng chéo
            double newNorm = 0.0;
            for (int k = 0; k < n; ++k) {
                newNorm += orthonormalVectors[k][i] * orthonormalVectors[k][i];
            }
            newNorm = sqrt(newNorm);

            for (int k = 0; k < n; ++k) {
                orthonormalVectors[k][i] /= newNorm;
            }
        }
    }

    return orthonormalVectors;
}

Matrix calculateSigma(const vector<double>& eigenValues, int m, int n) {
    Matrix Sigma(m, vector<double>(n, 0.0));

    for (int i = 0; i < eigenValues.size() && i < min(m, n); ++i) {
        if (eigenValues[i] > 0) { 
            Sigma[i][i] = sqrt(eigenValues[i]);
        }
    }

    return Sigma;
}

Matrix calculateU(const Matrix& A, const Matrix& V, const Matrix& Sigma) {
    int m = A.size();
    int pos = min(Sigma.size(),Sigma[0].size());
    while (Sigma[pos-1][pos-1] < 0.000001) {
        pos--;
    }
    Matrix U = Matrix(m,vector<double>(m,0.0));

    for (int i = 0; i < min(Sigma.size(),Sigma[0].size()); i++) {
        if (Sigma[i][i] > 0.000001) {
            vector<double> Ui = multiplyScalarMatrixVector((1.0/Sigma[i][i]),A,getVectorFromMatrix(V,i));
            for(int j = 0; j < m; j++) {
                U[j][i] = Ui[j];
            }
        }
    }
    
    while (pos < m) {
        vector<double> e3(m,0);
        e3[pos] = 1.0;
        vector<double> Ui = e3;
        for(int i = 0; i < pos; i++) {
            vector<double> ui = getVectorFromMatrix(U,i);
            double tmp = dotProduct(e3,ui) / (vectorNorm(ui)*vectorNorm(ui));
            ui = multiplyScalarVector(tmp, ui);
            Ui = subtractVector(Ui,ui);
        }


        long double chuanhoa = vectorNorm(Ui);

        for(int j = 0 ; j < m ; j++) {
            U[j][pos] = Ui[j] / chuanhoa;
        }
        pos++;
    }


    return U;
}

void svdDecomposition(Matrix& A) {
    int m = A.size();
    int n = A[0].size();

    Matrix AT = transposeMatrix(A);

    Matrix ATA = multiplyMatrix(AT,A);

    vector<double> EigenValues = qrAlgorithm(ATA); 

    Matrix V = eigenVectors(ATA,EigenValues);
    V = transposeMatrix(V);
    V = orthonormalize(V);

    Matrix Sigma = calculateSigma(EigenValues,m,n);

    Matrix U = calculateU(A, V, Sigma);

    roundMatrix(U);



    // Matrix U = calculateU(A, V, Sigma);
    // cout << "U : " << endl;
    // printMatrix(U);
    // U = orthonormalize(U);

    cout << "U : " << endl;
    printMatrix(U);
    cout << "Sigma:\n";
    printMatrix(Sigma);
    cout << "V:\n";
    printMatrix(V);
    cout << "V^T:\n";
    printMatrix(transposeMatrix(V));
    cout << "Product of U, Sigma, and V^T:\n";
    Matrix RS = multiplyMatrix(multiplyMatrix(U, Sigma), transposeMatrix(V));
    roundMatrix(RS);
    printMatrix(RS);

}

int main() {
    int m, n;
    cout << "Nhap so dong : "; cin >> m;
    cout << "Nhap so cot : "; cin >> n;
    Matrix A(m, vector<double>(n));
    cout << "Nhap ma tran:\n";
    for(int i = 0; i < m; i++) {
        for(int j = 0; j < n; j++) {
            cout << "A[" << i << "][" << j << "] = ";  
            cin >> A[i][j];
        }
    }
    cout << "Ma tran A:\n";
    printMatrix(A);
    svdDecomposition(A);
}