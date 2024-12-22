#include <iostream>
#include <vector>
#include <numeric> // Thư viện gcd
#include <tuple>   // Thư viện tuple để trả về nhiều giá trị
using namespace std;

tuple<long long, long long, long long> extendedGCD(long long a, long long b) {
    if (b == 0) return {a, 1, 0}; // gcd(a, b) = a; x = 1; y = 0
    auto [g, x1, y1] = extendedGCD(b, a % b);
    return {g, y1, x1 - (a / b) * y1};
}

// Hàm kiểm tra đôi một nguyên tố cùng nhau
bool arePairwiseCoprime(const vector<long long>& mods) {
    int n = mods.size();
    for (int i = 0; i < n; ++i) {
        for (int j = i + 1; j < n; ++j) {
            if (gcd(mods[i], mods[j]) != 1) {
                return false;
            }
        }
    }
    return true;
}
long long chineseRemainderTheorem(const vector<long long>& remainders, const vector<long long>& mods) {
    long long M = 1; // Tính tích các m_i
    for (auto mod : mods) M *= mod;

    long long result = 0;
    for (size_t i = 0; i < mods.size(); ++i) {
        long long Mi = M / mods[i];
        auto [g, MiInverse, _] = extendedGCD(Mi, mods[i]);
        if (g != 1) throw runtime_error("He phuong trinh khong thoa man đieu kien CRT");
        result += remainders[i] * Mi * MiInverse;
    }
    return (result % M + M) % M; // Đảm bảo kết quả là số dương
}

int main() {
    int k;
    cout << "Nhap so luong phuong trinh k: ";
    cin >> k;

    vector<long long> a(k), m(k);
    cout << "Nhap cac he so a1, a2, ..., ak: ";
    for (int i = 0; i < k; ++i) cin >> a[i];
    cout << "Nhap cac modulo m1, m2, ..., mk: ";
    for (int i = 0; i < k; ++i) cin >> m[i];

    // Kiểm tra điều kiện CRT
    if (!arePairwiseCoprime(m)) {
        cerr << "Cac modulo khong đoi mot nguyen to cung nhau. Khong the ap dung CRT.\n";
        return 1;
    }

    try {
        long long solution = chineseRemainderTheorem(a, m);
        cout << "Nghiem cua he phuong trinh thang du Trung Hoa la: " << solution << endl;
    } catch (const exception& e) {
        cerr << "Loi: " << e.what() << endl;
        return 1;
    }

    return 0;
}