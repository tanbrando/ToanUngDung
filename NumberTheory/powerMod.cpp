#include <iostream>
using namespace std;

// Hàm tính a^m mod n bằng thuật toán lũy thừa nhanh
long long powerMod(long long a, long long m, long long n) {
    long long result = 1;   // Kết quả ban đầu
    a = a % n;              // Đảm bảo a < n

    while (m > 0) {
        // Nếu m lẻ, nhân a vào kết quả
        if (m % 2 == 1) {
            result = (result * a) % n;
        }

        // Bình phương a và giảm m
        a = (a * a) % n;
        m /= 2;
    }

    return result;
}

int main() {
    long long a, m, n;
    cout << "Nhap a, m, n: ";
    cin >> a >> m >> n;

    long long result = powerMod(a, m, n);
    cout << a << "^" << m << " mod " << n << " = " << result << endl;

    return 0;
}