#pragma once
#include <vector>
#define typeT int64_t
using namespace std;

struct SmithStruct {
    vector<vector<typeT>> L;//Левая унимодальная матрица
    vector<vector<typeT>> S;//Нормальная форма Смита
    vector<vector<typeT>> R;//Правая унимодальная матрица
    typeT determinant;      //Детерминант матрицы
};

struct elemTabl {
    typeT minNormX; //Минимальная норма (сумма верхнего минимума и текущего Xm)
    typeT upMin;    //Номер верхнего минимума
    typeT Xm;       //Текущий Xm
};

class svp {
private:
    vector<vector<elemTabl>> tabl;  //Таблица для хранения промежуточных минимумов
    vector<vector<typeT>> transpP;  //Транспонированная матрица P
    vector<typeT> S;                //Диагональ матрицы S
    typeT delta;                    //Детерминант матрицы
    int n;                          //Размерность матрицы
    int sizeE;                      //Размерность единичной подматрицы
    int powP;                       //Степень нормы
    vector<typeT> minX;             //Кратчайший вектор

public:
    svp();
    ~svp();
    svp(vector<vector<typeT>> _matrix, int _powP = 2);  //Предобработка данных
    vector<typeT> StartSearch();                        //Запуск поиска минимального вектора

    vector<typeT> getMinX();
    typeT getMin();
    typeT getDelta();
    int getN();
    int getSizeE();
    int getPowP();
    void setPowP(int _newP);

private:
    typeT modS(typeT _b, typeT _s);//b по модулю s
    vector<typeT> multK(vector<typeT> _vector, typeT _koef);//произведение вектора на число
    vector<typeT> differenceVec(vector<typeT> _vec1, vector<typeT> _vec2, typeT _koef = 1);//разность векторов
    SmithStruct SmithNormalForm(vector<vector<typeT>> _matrix);//Приведение матрицы к нормальной форма Смита

    typeT AlgEuclid(typeT a, typeT b, typeT& x, typeT& y);  //Алгоритм Евклида
    typeT SearchMinOneLevel(typeT _b, int _m);              //Поиск минимума имея один столбец
    typeT RecursiveSearchMin(typeT _b, int _m);             //Рекурсивный поиск минимума
};
