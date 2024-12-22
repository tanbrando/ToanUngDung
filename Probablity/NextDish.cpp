#include <iostream>
#include <vector>
#include <map>
#include <string>
#include <iomanip>

using namespace std;

// Hàm tính ma trận chuyển đổi trạng thái
vector<vector<double>> calculateTransitionMatrix(const vector<int>& data, int numStates) {
    vector<vector<int>> transitionCount(numStates, vector<int>(numStates, 0));
    vector<int> stateCount(numStates, 0);

    int n = data.size();
    for (int i = 0; i < n - 1; i++) {
        stateCount[data[i]]++;
        transitionCount[data[i]][data[i + 1]]++;
    }
    // Thêm self-loop cho trạng thái cuối
    stateCount[data[n - 1]]++;
    transitionCount[data[n - 1]][data[n - 1]]++;

    // Tính xác suất chuyển đổi
    vector<vector<double>> transitionMatrix(numStates, vector<double>(numStates, 0.0));
    for (int i = 0; i < numStates; i++) {
        if (stateCount[i] > 0) {
            for (int j = 0; j < numStates; j++) {
                transitionMatrix[i][j] = static_cast<double>(transitionCount[i][j]) / stateCount[i];
            }
        }
    }

    return transitionMatrix;
}

// Hàm in ma trận
void printMatrix(const vector<vector<double>>& matrix, const vector<string>& stateNames) {
    cout << "Ma tran chuyen doi trang thai:\n";
    for (size_t i = 0; i < matrix.size(); i++) {
        cout << stateNames[i] << ": ";
        for (size_t j = 0; j < matrix[i].size(); j++) {
            cout << stateNames[j] << " (" << fixed << setprecision(2) << matrix[i][j] << ") ";
        }
        cout << endl;
    }
}

// Hàm dự đoán xác suất món ăn tiếp theo
void predictNextState(int currentState, const vector<vector<double>>& transitionMatrix, const vector<string>& stateNames) {
    cout << "\nXac suat mon an tiep theo sau '" << stateNames[currentState] << "':\n";
    for (size_t i = 0; i < transitionMatrix[currentState].size(); i++) {
        cout << stateNames[i] << ": " << fixed << setprecision(2) << transitionMatrix[currentState][i] * 100 << "%\n";
    }
}

int main() {
    // Ánh xạ từ tên món ăn sang chỉ số
    map<string, int> foodToIndex = {{"Banh mi", 0}, {"Com tam", 1}, {"Pho", 2}, {"Pizza", 3}};
    vector<string> indexToFood = {"Banh mi", "Com tam", "Pho", "Pizza"};

    // Dữ liệu đầu vào dạng chuỗi
    vector<string> data = {"Banh mi", "Pho", "Com tam", "Pizza", "Banh mi", "Pho", 
                           "Com tam", "Pizza", "Banh mi", "Banh mi", "Pho", "Banh mi"};

    // Chuyển dữ liệu sang dạng số nguyên
    vector<int> numericData;
    for (const string& food : data) {
        numericData.push_back(foodToIndex[food]);
    }

    // Tính ma trận chuyển đổi trạng thái
    int numStates = foodToIndex.size();
    vector<vector<double>> transitionMatrix = calculateTransitionMatrix(numericData, numStates);

    // In ma trận chuyển đổi trạng thái
    printMatrix(transitionMatrix, indexToFood);

    // Dự đoán xác suất món ăn tiếp theo sau "Phở"
    int currentState = foodToIndex["Pho"];
    predictNextState(currentState, transitionMatrix, indexToFood);

    return 0;
}