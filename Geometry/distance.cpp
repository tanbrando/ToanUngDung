#include <iostream>
#include <vector>
#include <cmath>
#include <numeric> // Để sử dụng hàm accumulate
using namespace std;

// Hàm tính khoảng cách Euclid
double euclideanDistance(const vector<double>& x, const vector<double>& y) {
    if (x.size() != y.size()) throw invalid_argument("Kích thước vector không khớp!");
    double sum = 0.0;
    for (size_t i = 0; i < x.size(); i++) {
        sum += pow(x[i] - y[i], 2);
    }
    return sqrt(sum);
}

// Hàm tính độ đo cosine similarity
double cosineSimilarity(const vector<double>& x, const vector<double>& y) {
    if (x.size() != y.size()) throw invalid_argument("Kích thước vector không khớp!");
    double dotProduct = 0.0, normX = 0.0, normY = 0.0;
    for (size_t i = 0; i < x.size(); i++) {
        dotProduct += x[i] * y[i];
        normX += x[i] * x[i];
        normY += y[i] * y[i];
    }
    return dotProduct / (sqrt(normX) * sqrt(normY));
}

// Hàm tính hệ số tương quan Pearson
double pearsonCorrelation(const vector<double>& x, const vector<double>& y) {
    if (x.size() != y.size()) throw invalid_argument("Kích thước vector không khớp!");
    double meanX = accumulate(x.begin(), x.end(), 0.0) / x.size();
    double meanY = accumulate(y.begin(), y.end(), 0.0) / y.size();
    double numerator = 0.0, denomX = 0.0, denomY = 0.0;
    for (size_t i = 0; i < x.size(); i++) {
        double diffX = x[i] - meanX;
        double diffY = y[i] - meanY;
        numerator += diffX * diffY;
        denomX += diffX * diffX;
        denomY += diffY * diffY;
    }
    return numerator / (sqrt(denomX) * sqrt(denomY));
}

int main() {
    vector<double> x = {1.0, 2.0, 3.0};
    vector<double> y = {4.0, 5.0, 6.0};

    try {
        cout << "Khoang cach Euclid: " << euclideanDistance(x, y) << endl;
        cout << "Do do Cosine Similarity: " << cosineSimilarity(x, y) << endl;
        cout << "Khoang cach Cosine: " << 1 - cosineSimilarity(x, y) << endl;
        cout << "He so tuong quan Pearson: " << pearsonCorrelation(x, y) << endl;
    } catch (const invalid_argument& e) {
        cerr << "Loi: " << e.what() << endl;
    }

    return 0;
}