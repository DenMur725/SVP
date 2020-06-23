#include "svp.h"
#include <ctime>
#include <chrono>
#include <iostream>

vector<typeT> DecomposingNumber(typeT _delta) {
    vector<typeT> result;
    while (_delta % 2 == 0) {
        result.push_back(2);
        _delta /= 2;
    }
    for (typeT i = 3; i <= sqrt(_delta); i += 2) {
        while (_delta % i == 0) {
            result.push_back(i);
            _delta /= i;
        }
    }
    if (_delta != 1) {
        result.push_back(_delta);
    }
    return result;
}

vector<vector<typeT>> GenHermiteForm(int _n, typeT _delta) {
    vector<vector<typeT>> result(_n, vector<typeT>(_n, 0));
    vector<typeT> delta = DecomposingNumber(_delta);

    for (int i = 0; i < _n; i++) {
        result[i][i] = 1;
    }
    int k = 0, j = 0;

    while (k < delta.size()) {
        result[_n - 1][_n - 1] *= delta[j];
        j++;
        k++;
        for (int i = _n - 2; i >= 0 && k < delta.size(); i--) {
            if (delta[j] != delta[j - 1]) {
                break;
            }
            result[i][i] *= delta[j];
            j++;
            k++;
        }
    }
    for (int i = 1; i < _n; i++) {
        if (result[i][i] == 1) {
            continue;
        }
        if (i >= result[i][i]) {
            typeT tmp = 0;
            int c = rand() % i;
            for (int l = 0; l < c; l++) {
                result[i][l] = tmp + (rand() % ((result[i][i] / 2) - tmp));
                tmp = result[i][l];
            }
            for (int l = c; l < i; l++) {
                result[i][l] = tmp + (rand() % (result[i][i] - tmp));
                tmp = result[i][l];
            }
            continue;
        }
        typeT sizeStep = result[i][i] / i;
        for (int l = 0; l < i; l++) {
            result[i][l] = sizeStep * l + (rand() % sizeStep);
        }
    }
    return result;
}

vector<vector<typeT>> DataMatrix(int i) {
    if (i == 0) {
        int n = 5;
        typeT delta = 2*5*5*5*11;
        vector<vector<typeT>> m = GenHermiteForm(n, delta);
        return m;
    }
    if (i == 1) {
        int n = 4;
        vector<vector<typeT>> m(n, vector<typeT>(n, 0));
        m[0][0] = 1;    m[0][1] = 0;    m[0][2] = 0;    m[0][3] = 0;
        m[1][0] = 0;    m[1][1] = 1;    m[1][2] = 0;    m[1][3] = 0;
        m[2][0] = 0;    m[2][1] = 0;    m[2][2] = 1;    m[2][3] = 0;
        m[3][0] = 4;    m[3][1] = 9;    m[3][2] = 11;   m[3][3] = 17;
        return m;
    }
    if (i == 2) {
        int n = 4;
        vector<vector<typeT>> m(n, vector<typeT>(n, 0));
        m[0][0] = 1;    m[0][1] = 0;    m[0][2] = 0;    m[0][3] = 0;
        m[1][0] = 1;    m[1][1] = 2;    m[1][2] = 0;    m[1][3] = 0;
        m[2][0] = 3;    m[2][1] = 5;    m[2][2] = 6;    m[2][3] = 0;
        m[3][0] = 4;    m[3][1] = 9;    m[3][2] = 23;   m[3][3] = 30;
        return m;
    }
    if (i == 3) {
        int n = 4;
        vector<vector<typeT>> m(n, vector<typeT>(n, 0));
        m[0][0] = 1;    m[0][1] = 0;    m[0][2] = 0;    m[0][3] = 0;
        m[1][0] = 1;    m[1][1] = 2;    m[1][2] = 0;    m[1][3] = 0;
        m[2][0] = 3;    m[2][1] = 5;    m[2][2] = 6;    m[2][3] = 0;
        m[3][0] = 4;    m[3][1] = 9;    m[3][2] = 11;   m[3][3] = 12;
        return m;
    }
    if (i == 4) {
        int n = 3;
        vector<vector<typeT>> m(n, vector<typeT>(n, 0));
        m[0][0] = 1;    m[0][1] = 0;    m[0][2] = 0;
        m[1][0] = 0;    m[1][1] = 1;    m[1][2] = 0;
        m[2][0] = 1909;    m[2][1] = 1929;    m[2][2] = 1954;
        return m;
    }
    if (i == 5) {
        int n = 5;
        vector<vector<typeT>> m(n, vector<typeT>(n, 0));
        m[0][0] = 1;    m[0][1] = 0;    m[0][2] = 0;    m[0][3] = 0;    m[0][4] = 0;
        m[1][0] = 0;    m[1][1] = 1;    m[1][2] = 0;    m[1][3] = 0;    m[1][4] = 0;
        m[2][0] = 0;    m[2][1] = 0;    m[2][2] = 1;    m[2][3] = 0;    m[2][4] = 0;
        m[3][0] = 0;    m[3][1] = 0;    m[3][2] = 0;    m[3][3] = 1;    m[3][4] = 0;
        m[4][0] = 43;    m[4][1] = 932;    m[4][2] = 938;   m[4][3] = 1010;   m[4][4] = 1042;
        return m;
    }
}

void TestRun() {
    vector<vector<typeT>> matA = DataMatrix(5);
    int n = matA.size();
    cout << "matA: " << endl;
    if (n < 30) {
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                cout << matA[i][j] << "   ";
            }
            cout << endl;
        }
    }

    svp svp1(matA);

    auto start_time = std::chrono::steady_clock::now();

    vector<typeT> resX = svp1.StartSearch();

    auto finish_time = std::chrono::steady_clock::now();
    auto elapsed_ms = std::chrono::duration_cast<std::chrono::milliseconds>(finish_time - start_time);
    std::cout << "The time: " << elapsed_ms.count() << " ms\n";

    cout << "resX: " << endl;
    for (int j = 0; j < n; j++) {
        cout << resX[j] << "   ";
    }
    cout << endl;
}

int main()
{
    srand(time(0));
    TestRun();
    return 0;
}
