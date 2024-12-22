#include <bits/stdc++.h> 

using namespace std;
const int MAX =100;
bool isSquare(int matrix[][MAX], int n) {
    return n > 0;
}

// Kiểm tra ma trận đối xứng
bool isSymmetric(int matrix[][MAX], int n) {
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            if (matrix[i][j] != matrix[j][i]) {
                return false;
            }
        }
    }
    return true;
}

// Kiểm tra ma trận xác định dương
bool isPositiveDefinite(int matrix[][MAX], int n) {
    for (int i = 0; i < n; i++) {
        if (matrix[i][i] <= 0) {
            return false;
        }
    }
    return true;
}

void Cholesky_Decomposition(int matrix[][MAX],int n)
{
    if (!isSquare(matrix, n)) {
        cout << "Matrix is not square!" << endl;
        return;
    }

    if (!isSymmetric(matrix, n)) {
        cout << "Matrix is not symmetric!" << endl;
        return;
    }

    if (!isPositiveDefinite(matrix, n)) {
        cout << "Matrix is not positive definite!" << endl;
        return;
    }
    int lower[n][n];
    memset(lower,0,sizeof(lower));
    for(int i =0;i < n;i++) {
        for(int j= 0; j <= i ; j++) {
            int sum=0;
            if (j== i){
                for(int k=0;k < j; k++) sum+= pow(lower[j][k],2);
                lower[j][j]=sqrt(matrix[j][j]- sum); 
            }else {
                for(int k=0;k < j;k++) sum+= (lower[i][k]* lower[j][k]);
                (lower[i][j]=matrix[i][j]- sum)/ lower[j][j]; 
            }
        }
    }
    cout<<setw(6)<< " Lower Triangular" <<setw(30)<< "Transpose"<<endl;
    for(int i =0;i < n;i++) {
        for(int j= 0; j <n; j++)cout <<setw(6)<< lower[i][j]<< "\t";cout <<"\t";
        for(int j= 0; j <n; j++)cout <<setw(6)<< lower[j][i]<< "\t";cout <<endl; 
    }
}

int main() {
    int n = 3;
    int matrix[][MAX] = {
        {15, 12, -5},
        {12, 35, -13},
        {-5, -13, 28}
    };
    Cholesky_Decomposition(matrix, n);
    return 0;
}