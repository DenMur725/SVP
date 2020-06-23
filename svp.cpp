#include "svp.h"

vector<typeT> svp::getMinX(){ return this->minX; }
typeT svp::getMin()         { return this->tabl[0][n - 1].minNormX; }
typeT svp::getDelta()       { return this->delta; }
int svp::getN()             { return this->n; }
int svp::getSizeE()         { return this->sizeE; }
int svp::getPowP()          { return this->powP; }
void svp::setPowP(int _newP){ this->powP = _newP; }

svp::svp()
{
    this->tabl = vector<vector<elemTabl>>(1, vector<elemTabl>(1));
    this->transpP = vector<vector<typeT>>(1, vector<typeT>(1));
    this->S = vector<typeT>(1);
    this->delta = 0;
    this->powP = 2;
    this->n = 1;
    this->sizeE = 1;
    this->minX = vector<typeT>(0);
}

svp::svp(vector<vector<typeT>> _matrix, int _powP)
{
    n = _matrix.size();
    powP = _powP;
    SmithStruct smithForm = SmithNormalForm(_matrix);
    delta = smithForm.determinant;
    if (delta < 0)
        delta *= -1;
    tabl = vector<vector<elemTabl>>(delta, vector<elemTabl>(n));
    transpP = vector<vector<typeT>>(n, vector<typeT>(n));
    tabl[0][n - 1].minNormX = 0;
    for (int i = 0; i < n; i++) {
        S.push_back(smithForm.S[i][i]);
        minX.push_back(_matrix[i][0]);
        tabl[0][n - 1].minNormX += pow(abs(_matrix[i][0]), _powP);
        for (int j = 0; j < n; j++) {
            transpP[j][i] = modS(smithForm.L[i][j], S[i]);
        }
    }

    for (int i = 0; i < delta; i++) {
        for (typeT j = 0; j < n; j++) {
            tabl[i][j].upMin = -1;
        }
    }
    sizeE = n - 1;
    for (int i = 0; i < n; i++) {
        if (S[i] != 1) {
            sizeE = i;
            break;
        }
    }
}

svp::~svp()
{
}

vector<typeT> svp::StartSearch() {
    typeT resultMin = RecursiveSearchMin(0, n - 1);
    typeT idUp = tabl[0][n - 1].upMin;
    if (idUp == -1) {
        return minX;
    }
    vector<typeT> resultMinX(n, 0);
    resultMinX[n - 1] = tabl[0][n - 1].Xm;
    for (int i = n - 2; i >= 0; i--) {
        if (idUp == -2) {
            break;
        }
        resultMinX[i] = tabl[idUp][i].Xm;
        idUp = tabl[idUp][i].upMin;
    }
    minX = resultMinX;
    return resultMinX;
}




typeT svp::AlgEuclid(typeT a, typeT b, typeT& x, typeT& y) {
    if (b == 0) {
        x = 1;
        y = 0;
        return a;
    }
    typeT x1 = 0, y1 = 0;
    typeT d1;

    d1 = AlgEuclid(b, a % b, x1, y1);
    x = y1;
    y = x1 - ((a / b) * y1);
    return d1;
}

typeT svp::SearchMinOneLevel(typeT _b, int _m) {
    typeT result = 1, cur;
    typeT sumP = 0, sumB = 0;
    int flag = 0;
    for (int i = sizeE; i < n; i++) {
        if (transpP[_m][i] != modS(_b, S[i])) {
            flag = 1;
            break;
        }
    }
    if (flag == 0) {
        return 1;
    }
    for (int i = sizeE; i < n; i++) {
        typeT bModS = modS(_b, S[i]);
        if (transpP[_m][i] == 0 && bModS != 0) {
            return delta;
        }
        typeT koef = S[n - 1] / S[i];
        sumP = modS(sumP + transpP[_m][i] * koef, S[n - 1]);
        sumB = modS(sumB + bModS * koef, S[n - 1]);
    }
    typeT x, y;
    typeT d = AlgEuclid(sumP, S[n - 1], x, y);
    result = sumB * x;
    return result;
}

typeT svp::RecursiveSearchMin(typeT _b, int _m) {
    if (_b == 8 && _m == 4) {
        int t = 1;
    }
    if (tabl[_b][_m].upMin != -1) {
        return tabl[_b][_m].minNormX;
    }

    if (_m == 0) {
        tabl[_b][0].Xm = SearchMinOneLevel(_b, 0);
        if (tabl[_b][0].Xm == delta || tabl[_b][0].Xm == 0) {
            return -1;
        }
        tabl[_b][0].upMin = _b;
        tabl[_b][0].minNormX = pow(abs(tabl[_b][0].Xm), powP);

        return tabl[_b][0].minNormX;
    }
    elemTabl minElem;
    typeT min = RecursiveSearchMin(_b, _m - 1);
    minElem.minNormX = min;
    minElem.upMin = _b;
    minElem.Xm = 0;

    ////
    typeT XmOneLevel = SearchMinOneLevel(_b, _m);
    if (XmOneLevel != delta && XmOneLevel != 0) {
        typeT minOneLevel = pow(abs(XmOneLevel), powP);
        if (minOneLevel <= min || min == -1) {
            min = minOneLevel;
            minElem.minNormX = min;
            minElem.upMin = -2;
            minElem.Xm = XmOneLevel;
        }
    }
    ////
    typeT delta2 = (delta / 2) + 1;

    for (typeT l = 1; l < delta2; l++) {
        typeT bPx = -1;
        //test new b
        vector<typeT> tmpNewB(n, 0);
        vector<typeT> tmpTestRealB(n, 0);
        for (int i = sizeE; i < n; i++) {
            tmpNewB[i] = modS(_b - l * transpP[_m][i], S[i]);
        }
        for (int i = sizeE; i < n; i++) {
            tmpTestRealB[i] = modS(tmpNewB[n - 1], S[i]);
        }
        for (int i = sizeE; i < n; i++) {
            if (tmpNewB[i] != tmpTestRealB[i]) {
                bPx = -2;
                break;
            }
        }
        if (bPx == -2) {
            continue;
        }
        //
        bPx = tmpNewB[n - 1];
        if (bPx == _b) {
            break;
        }
        typeT tmpMin = RecursiveSearchMin(bPx, _m - 1);
        if (tmpMin == -1) {
            continue;
        }
        tmpMin += static_cast<typeT>(pow(l, powP));
        if (tmpMin < min || min == -1) {
            min = tmpMin;
            minElem.minNormX = min;
            minElem.upMin = bPx;
            minElem.Xm = l;
        }

    }

    /////-delta

    for (typeT l = -1; l > -delta2; l--) {
        typeT bPx = -1;
        //test new b
        vector<typeT> tmpNewB(n, 0);
        vector<typeT> tmpTestRealB(n, 0);
        for (int i = sizeE; i < n; i++) {
            tmpNewB[i] = modS(_b - l * transpP[_m][i], S[i]);
        }
        for (int i = sizeE; i < n; i++) {
            tmpTestRealB[i] = modS(tmpNewB[n - 1], S[i]);
        }
        for (int i = sizeE; i < n; i++) {
            if (tmpNewB[i] != tmpTestRealB[i]) {
                bPx = -2;
                break;
            }
        }
        if (bPx == -2) {
            continue;
        }
        //
        bPx = tmpNewB[n - 1];
        if (bPx == _b) {
            break;
        }
        typeT tmpMin = RecursiveSearchMin(bPx, _m - 1);
        if (tmpMin == -1) {
            continue;
        }
        tmpMin += static_cast<typeT>(pow(-l, powP));
        if (tmpMin < min || min == -1) {
            min = tmpMin;
            minElem.minNormX = min;
            minElem.upMin = bPx;
            minElem.Xm = l;
        }

    }

    /////

    if (min == -1) {
        return -1;
    }
    tabl[_b][_m].minNormX = minElem.minNormX;
    tabl[_b][_m].upMin = minElem.upMin;
    tabl[_b][_m].Xm = minElem.Xm;

    return tabl[_b][_m].minNormX;
}




typeT svp::modS(typeT _b, typeT _s) {
    if (_b == 0 || _s == 1)
        return 0;
    typeT res = _b - ((_b / _s) * _s);
    if (res == 0) {
        return 0;
    }
    if (_b > 0) {
        return res;
    }
    else {
        return _s + res;
    }
}

vector<typeT> svp::multK(vector<typeT> _vector, typeT _koef) {
    int size = _vector.size();
    vector<typeT> result = _vector;
    for (int j = 0; j < size; j++) {
        result[j] *= _koef;
    }
    return result;
}

vector<typeT> svp::differenceVec(vector<typeT> _vec1, vector<typeT> _vec2, typeT _koef) {
    int size1 = _vec1.size();
    int size2 = _vec2.size();
    if (size1 == size2) {
        vector<typeT> result(size1, 0);
        for (int j = 0; j < size1; j++) {
            result[j] = _vec1[j] - (_koef * _vec2[j]);
        }
        return result;
    }
    return vector<typeT>(1, 0);
}

SmithStruct svp::SmithNormalForm(vector<vector<typeT>> _matrix) {
    int row = _matrix.size();
    int col = _matrix[0].size();

    vector<vector<typeT>> result = _matrix;
    int k = 0, flag1 = 3, flag2 = 0;
    typeT determinant = 1;
    ///
    vector<vector<typeT>> resultL(col, vector<typeT>(row, 0));
    vector<vector<typeT>> resultR(col, vector<typeT>(row, 0));
    for (int i = 0; i < row && i < col; i++) {
        resultL[i][i] = 1;
        resultR[i][i] = 1;
    }
    ///
    while (k < row && k < col) {
        if (flag1 == 3) {
            int minI, minJ;
            typeT min;
            flag2 = 0;
            for (int i = k; i < row; i++) {
                for (int j = k; j < col; j++) {
                    if (result[i][j] != 0) {
                        flag2 = 1;
                        minI = i;
                        minJ = j;
                        min = abs(result[i][j]);
                        break;
                    }
                }
                if (flag2 == 1) {
                    break;
                }
            }
            if (flag2 == 0) {
                break;
            }
            for (int i = minI; i < row; i++) {
                for (int j = minJ; j < col; j++) {
                    if (result[i][j] != 0 && abs(result[i][j]) < min) {
                        minI = i;
                        minJ = j;
                        min = abs(result[i][j]);
                    }
                }
            }

            if (k != minJ) {//swap col
                determinant *= -1;
                for (int i = 0; i < row; i++) {
                    typeT tmp = result[i][k];
                    result[i][k] = result[i][minJ];
                    result[i][minJ] = tmp;
                }
            }
            ///
            if (k != minI) {
                resultL[minI].swap(resultL[k]);//col
            }
            if (k != minJ) {
                for (int i = 0; i < col; i++) {//row
                    typeT tmp = resultR[i][k];
                    resultR[i][k] = resultR[i][minJ];
                    resultR[i][minJ] = tmp;
                }
            }
            ///

            if (k != minI) {//swap row
                determinant *= -1;
                result[minI].swap(result[k]);
            }
            if (result[k][k] < 0) {
                result[k] = multK(result[k], -1);
                determinant *= -1;
                ///
                resultL[k] = multK(resultL[k], -1);
                ///
            }
            flag1 = 4;
        }
        if (flag1 == 4) {
            for (int i = k + 1; i < row; i++) {
                if (result[i][k] != 0) {
                    typeT koef = result[i][k] / result[k][k];
                    result[i] = differenceVec(result[i], result[k], koef);
                    ///
                    resultL[i] = differenceVec(resultL[i], resultL[k], koef);
                    ///
                    flag1 = 3;
                    break;
                }
            }
            if (flag1 == 4) {
                flag1 = 5;
            }
        }
        if (flag1 == 5) {
            for (int j = k + 1; j < col; j++) {
                if (result[k][j] != 0) {
                    typeT koef = result[k][j] / result[k][k];
                    for (int i = 0; i < row; i++) {
                        result[i][j] -= result[i][k] * koef;
                    }
                    ///
                    for (int i = 0; i < col; i++) {
                        resultR[i][j] -= resultR[i][k] * koef;
                    }
                    ///
                    flag1 = 3;
                    break;
                }
            }
            if (flag1 == 5) {
                flag1 = 6;
            }
        }
        if (flag1 == 6) {
            for (int i = k + 1; i < row; i++) {
                for (int j = k + 1; j < col; j++) {
                    if (result[i][j] % result[k][k] != 0) {
                        for (int l = 0; l < row; l++) {
                            result[l][k] += result[l][j];
                        }
                        ///
                        for (int l = 0; l < col; l++) {
                            resultR[l][k] += resultR[l][j];
                        }
                        ///
                        flag1 = 4;
                        break;
                    }
                }
                if (flag1 == 4) {
                    break;
                }
            }
            if (flag1 == 6) {
                k++;
                flag1 = 3;
            }
        }
    }

    if (row == col) {
        for (int i = 0; i < row; i++) {
            determinant *= result[i][i];
        }
    }

    SmithStruct res;
    res.L = resultL;
    res.S = result;
    res.R = resultR;
    res.determinant = determinant;
    return res;
}
