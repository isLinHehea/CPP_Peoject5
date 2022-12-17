#ifndef __MATRIX_HPP__
#include "matrix.hpp"
#endif

int main()
{
    Matrix <unsigned char> matrix_char(2,2);
    Matrix <short> matrix_short(2,2);
    Matrix <int> matrix_int(2,2);
    Matrix <float> matrix_float(2,2);
    Matrix <double> matrix_double(2,2);
    unsigned char data1[] = {'d','f','g','h'};
    short data2[] = {0,1,2,3};
    int data3[] = {0,1,2,3};
    float data4[] = {0,1,2,3};
    double data5[] = {0,1,2,3};
    matrix_char.value(data1);
    matrix_short.value(data2);
    matrix_int.value(data3);
    matrix_float.value(data4);
    matrix_double.value(data5);
    cout<<matrix_char;
    cout<<matrix_short;
    cout<<matrix_int;
    cout<<matrix_float;
    cout<<matrix_double;


    // Matrix<int> matrixA(5, 2);
    // Matrix<int> matrixA(5, 5);
    // Matrix<float> *matrixA = new Matrix<float>(3, 3);
    // Matrix<int> matrixB(2, 2);
    // Matrix<int> matrixB(2, 5);
    // int dataA[] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9};
    // int dataA[] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25};
    // float dataA[] = {1, 2, 3, 4, 5, 6, 7, 8, 10};
    // int dataB[] = {9, 8, 7, 6, 5, 4, 3, 2, 1, 0};
    // matrixA.value(dataA);
    // Matrix<int> *matrixE = matrixA.ROI(2, 2, 2, 2);
    // matrixB.value(dataB);
    // Matrix<int> matrixC  = matrixA + matrixB;
    // Matrix<int> matrixD  = matrixA - matrixB;
    // Matrix<int> matrixE  = matrixA * matrixB;
    // Matrix<int> matrixC  = matrixA + 9;
    // Matrix<int> matrixD  = matrixA - 9;
    // Matrix<int> matrixE  = matrixA * 9;
    // *matrixA = *matrixB;
    // matrixA->reshape(2, 5);
    // matrixB->flatten();
    // Matrix<int> *matrixC = matrixA->inverse();
    // Matrix<float> *matrixD = new Matrix<float>(3, 3);
    // float dataD[] = {1, 2, 3, 4, 5, 6, 7, 8, 10};
    // matrixD->value(dataD);

    // cout << matrixA;
    // cout << *matrixE;
    // cout << matrixC;
    // cout << matrixD;
    // cout << matrixE;
    // cout << (*matrixA == *matrixB) << endl;
    // cout << *(matrixD->inverse());
    // cout << *(matrixA->transpose());
    // cout << matrixA->determinant() << endl;
    // cout <<*(matrixA->inverse());

    // return 0;
}
// Matrix<float> matrixA(5, 2);
// Matrix<float> *matrixB = new Matrix<float>(2, 5);
// float dataA[] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9};
// float dataB[] = {9, 8, 7, 6, 5, 4, 3, 2, 1, 0};
// matrixA.setData(dataA);
// matrixB->setData(dataB);
// matrixA.print();
// matrixB->print();
// Matrix<float> matrixC(matrixA);
// Matrix<float> matrixD = matrixA;d
// Matrix<float> matrixE = Matrix<float>(matrixA);
// Matrix<float> *matrixF = new Matrix<float>(matrixA);
// Matrix<float> *matrixG = matrixB;
// Matrix<float> *matrixH(matrixB);
// Matrix<float> *matrixI = new Matrix<float>(*matrixB);
// Matrix<float> matrixJ(*matrixB);
// matrixC.print();
// matrixD.print();
// matrixE.print();
// matrixF->print();
// matrixG->print();
// matrixH->print();
// matrixI->print();
// matrixJ.print();
// cin >> matrixA;
// cout << matrixA;
// unsigned char, short, int, float, double
// int v1 = 1;
// unsigned char v2 = 'a';
// double v3 = 1;
// float v4 = 1.1;
// short v5 = 1;
// string v6 = "hh";
// cout << typeid(v1).name() << endl;
// cout << typeid(v2).name() << endl;
// cout << typeid(v3).name() << endl;
// cout << typeid(v4).name() << endl;
// cout << typeid(v5).name() << endl;
// cout << typeid(v6).name() << endl;
// Matrix<float> matrixA(2, 5);
// Matrix<float> *matrixD = new Matrix<float>(2, 5);
// cout << matrixD;

// Matrix<unsigned char> matrixA(5, 2);
// Matrix<float> matrixB(2, 2);
// Matrix<float> matrixA(2, 5);
// Matrix<float> matrixB(2, 5);
// Matrix<float> *matrixD = new Matrix<float>(2, 5);
// Matrix<float> matrixC(&matrixA);

// float dataA[] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9};
// float dataB[] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9};
// matrixA.value(dataA);
// matrixB.value(dataB);
// Matrix<float> matrixC(2, 5);
// matrixC = matrixA + matrixB;
// matrixA += matrixC;
// cout << matrixC;
// cout << matrixA;

// if (matrixA != matrixB)
//     cout << 1 << endl;
// else
//     cout << 2 << endl;

// Matrix<float> matrixA(5, 2);
// Matrix<float> matrixB(2, 5);
// float dataA[] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9};
// float dataB[] = {9, 8, 7, 6, 5, 4, 3, 2, 1, 0};
// matrixA.value(dataA);
// matrixB.value(dataB);
// Matrix<float> matrixC = matrixA * matrixB;
// cout << matrixC;
// Matrix<int> matrixA(5, 2);
// Matrix<int> matrixB(2, 5);
// int dataA[] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9};
// int dataB[] = {9, 8, 7, 6, 5, 4, 3, 2, 1, 0};
// matrixA.value(dataA);
// matrixB.value(dataB);
// Matrix<int> *matrixC = matrixA.subMatrix(1, 1, 2, 2);
// cout << matrixA;
// cout << matrixB;
// cout << *matrixC;
// Matrix<int> matrixA(2, 2);
// Matrix<int> matrixB(2, 2);
// int dataA[] = {1, 2, 3, 4};
// int dataB[] = {1, 2, 3, 4};
// matrixA.value(dataA);
// matrixB.value(dataB);
// Matrix<int> matrixC = matrixA * matrixB;
// cout << matrixC;
// Matrix<short> matrixA(2, 2);
// Matrix<short> matrixB(2, 2);
// short dataA[] = {1, 2, 3, 4};
// short dataB[] = {1, 2, 3, 4};
// matrixA.value(dataA);
// matrixB.value(dataB);
// Matrix<short> matrixC = matrixA * matrixB;
// cout << matrixC;
// cout << *matrixB;
// cout << *matrixC;
// cout << (*matrixA == *matrixB) << endl;
// cout << *(matrixD->inverse());
// cout << matrixD->determinant() << endl;

// return 0;
// }
// Matrix<float> matrixA(5, 2);
// Matrix<float> *matrixB = new Matrix<float>(2, 5);
// float dataA[] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9};
// float dataB[] = {9, 8, 7, 6, 5, 4, 3, 2, 1, 0};
// matrixA.setData(dataA);
// matrixB->setData(dataB);
// matrixA.print();
// matrixB->print();
// Matrix<float> matrixC(matrixA);
// Matrix<float> matrixD = matrixA;d
// Matrix<float> matrixE = Matrix<float>(matrixA);
// Matrix<float> *matrixF = new Matrix<float>(matrixA);
// Matrix<float> *matrixG = matrixB;
// Matrix<float> *matrixH(matrixB);
// Matrix<float> *matrixI = new Matrix<float>(*matrixB);
// Matrix<float> matrixJ(*matrixB);
// matrixC.print();
// matrixD.print();
// matrixE.print();
// matrixF->print();
// matrixG->print();
// matrixH->print();
// matrixI->print();
// matrixJ.print();
// cin >> matrixA;
// cout << matrixA;
// unsigned char, short, int, float, double
// int v1 = 1;
// unsigned char v2 = 'a';
// double v3 = 1;
// float v4 = 1.1;
// short v5 = 1;
// string v6 = "hh";
// cout << typeid(v1).name() << endl;
// cout << typeid(v2).name() << endl;
// cout << typeid(v3).name() << endl;
// cout << typeid(v4).name() << endl;
// cout << typeid(v5).name() << endl;
// cout << typeid(v6).name() << endl;
// Matrix<float> matrixA(2, 5);
// Matrix<float> *matrixD = new Matrix<float>(2, 5);
// cout << matrixD;

// Matrix<unsigned char> matrixA(5, 2);
// Matrix<float> matrixB(2, 2);
// Matrix<float> matrixA(2, 5);
// Matrix<float> matrixB(2, 5);
// Matrix<float> *matrixD = new Matrix<float>(2, 5);
// Matrix<float> matrixC(&matrixA);

// float dataA[] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9};
// float dataB[] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9};
// matrixA.value(dataA);
// matrixB.value(dataB);
// Matrix<float> matrixC(2, 5);
// matrixC = matrixA + matrixB;
// matrixA += matrixC;
// cout << matrixC;
// cout << matrixA;

// if (matrixA != matrixB)
//     cout << 1 << endl;
// else
//     cout << 2 << endl;

// Matrix<float> matrixA(5, 2);
// Matrix<float> matrixB(2, 5);
// float dataA[] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9};
// float dataB[] = {9, 8, 7, 6, 5, 4, 3, 2, 1, 0};
// matrixA.value(dataA);
// matrixB.value(dataB);
// Matrix<float> matrixC = matrixA * matrixB;
// cout << matrixC;
// Matrix<int> matrixA(5, 2);
// Matrix<int> matrixB(2, 5);
// int dataA[] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9};
// int dataB[] = {9, 8, 7, 6, 5, 4, 3, 2, 1, 0};
// matrixA.value(dataA);
// matrixB.value(dataB);
// Matrix<int> *matrixC = matrixA.subMatrix(1, 1, 2, 2);
// cout << matrixA;
// cout << matrixB;
// cout << *matrixC;
// Matrix<int> matrixA(2, 2);
// Matrix<int> matrixB(2, 2);
// int dataA[] = {1, 2, 3, 4};
// int dataB[] = {1, 2, 3, 4};
// matrixA.value(dataA);
// matrixB.value(dataB);
// Matrix<int> matrixC = matrixA * matrixB;
// cout << matrixC;
// Matrix<short> matrixA(2, 2);
// Matrix<short> matrixB(2, 2);
// short dataA[] = {1, 2, 3, 4};
// short dataB[] = {1, 2, 3, 4};
// matrixA.value(dataA);
// matrixB.value(dataB);
// Matrix<short> matrixC = matrixA * matrixB;
// cout << matrixC;