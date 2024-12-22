#include <iostream>
#include <tuple>
#include <cmath>
using namespace std;

// Hàm tính gcd và tìm nghiệm cơ bản (x, y) của phương trình a * x + b * y = gcd(a, b)
tuple<long long, long long, long long> extendedGCD(long long a, long long b) {
    if (b == 0) return {a, 1, 0};
    auto [g, x1, y1] = extendedGCD(b, a % b);
    return {g, y1, x1 - (a / b) * y1};
}

bool solveDiophantine(long long a, long long b, long long c, long long& x, long long& y) {
    auto [g, x0, y0] = extendedGCD(abs(a), abs(b));

    // Kiểm tra điều kiện có nghiệm
    if (c % g != 0) return false;

    // Tìm nghiệm cơ bản
    x = x0 * (c / g);
    y = y0 * (c / g);

    // Điều chỉnh dấu nếu cần
    if (a < 0) x = -x;
    if (b < 0) y = -y;

    return true;
}
int gcd(int a, int b) {
    while (b != 0) {
        int r = a % b;
        a = b;
        b = r;
    }
    return a;
}

int main() {
    long long a, b, c;
    cout << "Nhap a, b, c: ";
    cin >> a >> b >> c;

    long long x, y;
    if (solveDiophantine(a, b, c, x, y)) {
        cout << "Phuong trinh co nghiem: x = " << x << ", y = " << y << endl;
        cout << "Dang nghiem tong quat:" << endl;
        cout << "x = " << x << " + k * " << b / gcd(a, b) << endl;
        cout << "y = " << y << " - k * " << a / gcd(a, b) << endl;
    } else {
        cout << "Phuong trinh khong co nghiem nguyen." << endl;
    }

    return 0;
}