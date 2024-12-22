#include <iostream>
#include <vector>
#include <cmath>
#include <iomanip>

using namespace std;

void printFibonacci(int n) {
    vector<int> fibonacci;
    int a = 0, b = 1;
    while (a < n) {
        fibonacci.push_back(a);
        int next = a + b;
        a = b;
        b = next;
    }

    cout << "Day Fibonacci nho hon " << n << ": ";
    for (int num : fibonacci) {
        cout << num << " ";
    }
    cout << endl;
}

vector<bool> sieveOfEratosthenes(int limit) {
    vector<bool> isPrime(limit + 1, true);
    isPrime[0] = isPrime[1] = false; // 0 và 1 không phải số nguyên tố
    for (int i = 2; i * i <= limit; ++i) {
        if (isPrime[i]) {
            for (int j = i * i; j <= limit; j += i) {
                isPrime[j] = false;
            }
        }
    }
    return isPrime;
}

// Hàm tìm số nguyên tố gần nhất với N
int findNearestPrime(int n, const vector<bool>& isPrime) {
    int lower = n, upper = n;
    while (true) {
        if (lower > 1 && isPrime[lower]) return lower; // Số nguyên tố phía dưới
        if (upper < isPrime.size() && isPrime[upper]) return upper; // Số nguyên tố phía trên
        --lower;
        ++upper;
    }
}

int main() {
    int n; 
    cout << "Nhap so nguyen n: ";  
    cin >> n;
    printFibonacci(n);
    vector<bool> isPrime = sieveOfEratosthenes(n * 10); // Tạo sieve gần đúng với n * 10
    int nearestPrime = findNearestPrime(n, isPrime);
    cout << "So nguyen to gan nhat voi " << n << " la: " << nearestPrime << endl;
    return 0;
}
