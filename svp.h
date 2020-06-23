#pragma once
#include <vector>
#define typeT int64_t
using namespace std;

struct SmithStruct {
    vector<vector<typeT>> L;//����� ������������ �������
    vector<vector<typeT>> S;//���������� ����� �����
    vector<vector<typeT>> R;//������ ������������ �������
    typeT determinant;      //����������� �������
};

struct elemTabl {
    typeT minNormX; //����������� ����� (����� �������� �������� � �������� Xm)
    typeT upMin;    //����� �������� ��������
    typeT Xm;       //������� Xm
};

class svp {
private:
    vector<vector<elemTabl>> tabl;  //������� ��� �������� ������������� ���������
    vector<vector<typeT>> transpP;  //����������������� ������� P
    vector<typeT> S;                //��������� ������� S
    typeT delta;                    //����������� �������
    int n;                          //����������� �������
    int sizeE;                      //����������� ��������� ����������
    int powP;                       //������� �����
    vector<typeT> minX;             //���������� ������

public:
    svp();
    ~svp();
    svp(vector<vector<typeT>> _matrix, int _powP = 2);  //������������� ������
    vector<typeT> StartSearch();                        //������ ������ ������������ �������

    vector<typeT> getMinX();
    typeT getMin();
    typeT getDelta();
    int getN();
    int getSizeE();
    int getPowP();
    void setPowP(int _newP);

private:
    typeT modS(typeT _b, typeT _s);//b �� ������ s
    vector<typeT> multK(vector<typeT> _vector, typeT _koef);//������������ ������� �� �����
    vector<typeT> differenceVec(vector<typeT> _vec1, vector<typeT> _vec2, typeT _koef = 1);//�������� ��������
    SmithStruct SmithNormalForm(vector<vector<typeT>> _matrix);//���������� ������� � ���������� ����� �����

    typeT AlgEuclid(typeT a, typeT b, typeT& x, typeT& y);  //�������� �������
    typeT SearchMinOneLevel(typeT _b, int _m);              //����� �������� ���� ���� �������
    typeT RecursiveSearchMin(typeT _b, int _m);             //����������� ����� ��������
};
