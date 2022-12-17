#include <iostream>
#include <cstring>
#include <cstdio>
using namespace std;

#ifndef __MATRIX_H__
#define __MATRIX_H__

template <typename T>
class Matrix
{
    // friend function
    template <typename T_>
    friend ostream &operator<<(ostream &os, const Matrix<T_> &matrix);
    template <typename T_>
    friend istream &operator>>(istream &is, Matrix<T_> &matrix);
    template <typename T_>
    friend Matrix operator+(const Matrix<T_> &matrixLeft, const Matrix<T_> &matrixRight);
    template <typename T_>
    friend Matrix operator-(const Matrix<T_> &matrixLeft, const Matrix<T_> &matrixRight);
    template <typename T_>
    friend Matrix operator*(const Matrix<T_> &matrixLeft, const Matrix<T_> &matrixRight);

private:
    size_t row;
    size_t column;
    size_t size;
    // will do it
    // size_t channel;
    T *data;
    size_t *dataref;
    string type;

public:
    Matrix();
    Matrix(size_t row, size_t column);
    Matrix(const Matrix &matrix);
    Matrix &operator=(const Matrix &matrix);
    size_t getRow();

    template <typename T_>
    bool operator==( Matrix<T_> &matrix);
    Matrix &operator+=(const Matrix &matrix);
    Matrix &operator-=(const Matrix &matrix);
    Matrix &operator+(const Matrix &matrix);
    Matrix &operator-(const Matrix &matrix);
    Matrix &operator*(const Matrix &matrix);
    ~Matrix();  

    bool setData(const T *data);
    bool print();
};
// Constructor
template <class T>
Matrix<T>::Matrix()
{
    row = column = size = 0;
    data = nullptr;
    dataref = nullptr;
    type = typeid(T).name();
}

template <class T>
Matrix<T>::Matrix(size_t row, size_t column)
{
    this->row = row;
    this->column = column;
    size = this->row * this->column;
    data = new T[size * sizeof(T)]{};
    dataref = new size_t(1);
    type = typeid(T).name();
}

template <class T>
Matrix<T>::Matrix(const Matrix &matrix)
{
    row = matrix.row;
    column = matrix.column;
    size = row * column;
    data = matrix.data;
    dataref = matrix.dataref;
    *dataref++;
    type = matrix.type;
}

template <class T>
Matrix<T> &Matrix<T>::operator=(const Matrix &matrix)
{
    if (--*dataref == 0)
    {
        delete[] data;
        delete dataref;
    }
    row = matrix.row;
    column = matrix.column;
    size = row * column;
    data = matrix.data;
    dataref = matrix.dataref;
    *dataref++;
    type = matrix.type;
    return *this;
}

// SetData
template <class T>
bool Matrix<T>::setData(const T *data)
{
    if (this->data == NULL)
    {
        cerr << "In File " << __FILE__ << ", Line " << __LINE__ << ", Function " << __FUNCTION__ << "():\nError: The data of matrix to be set data is NULL." << endl;
        return false;
    }
    else if (data == NULL)
    {
        cerr << "In File " << __FILE__ << ", Line " << __LINE__ << ", Function " << __FUNCTION__ << "():\nError: The data to be set is NULL." << endl;
        return false;
    }
    else
        memcpy(this->data, data, sizeof(T) * this->size);
    return true;
}

// Print
template <class T>
bool Matrix<T>::print()
{
    if (data == NULL)
    {
        cerr << "In File " << __FILE__ << ", Line " << __LINE__ << ", Function " << __FUNCTION__ << "():\nError: The data of matrix to be printed is NULL." << endl;
        return NULL;
    }
    else
    {
        cout << "Matrix Type " << type << endl;
        for (size_t i = 0; i < row; i++)
        {
            cout << "[";
            for (size_t j = 0; j < column; j++)
            {
                cout << data[i * column + j] << ", ";
            }
            cout << "\b\b]" << endl;
        }
        return true;
    }
}

// Destructor
template <class T>
Matrix<T>::~Matrix()
{
    row = column = size = 0;
    if (data == NULL)
        cout << "In File " << __FILE__ << ", Line " << __LINE__ << ", Function " << __FUNCTION__ << "():\nWarning: The data of matrix to be destructed is NULL." << endl;
    else if (--*dataref == 0)
    {
        delete[] data;
        delete dataref;
    }
}

// -----------------------------------------------------------------------------------------------------------
template <class T>
template <class T_>
bool Matrix<T>::operator==(Matrix<T_> &matrix)
{
    if (type != matrix.type)
        return false;
    else if (row != matrix.getRow() || column != matrix.column)
        return false;
    else
    {
        for (size_t i = 0; i < size; i++)
            if (data[i] != matrix.data[i])
                return false;
        return true;
    }
}

template <class T>
Matrix<T> &Matrix<T>::operator+=(const Matrix &matrix)
{
    if (type != matrix.type)
        cout << "In File " << __FILE__ << ", Line " << __LINE__ << ", Function " << __FUNCTION__ << "():\nError: The data of matrix to be added is NULL" << endl;
    else if (row != matrix.row || column != matrix.column)
        cout << "In File " << __FILE__ << ", Line " << __LINE__ << ", Function " << __FUNCTION__ << "():\nError: The two matrices to be added have the different size." << endl;
    else
    {
        for (size_t i = 0; i < size; i++)
            data[i] += matrix.data[i];
    }
    return *this;
}

template <class T>
Matrix<T> &Matrix<T>::operator-=(const Matrix &matrix)
{
    for (size_t i = 0; i < size; i++)
        data[i] -= matrix.data[i];
    return *this;
}

template <class T>
Matrix<T> &Matrix<T>::operator+(const Matrix &matrix)
{
    Matrix<T> *matrixAfter = new Matrix<T>(this->row, this->column);
    for (size_t i = 0; i < size; i++)
        matrixAfter[i] = data[i] + matrix.data[i];
    return matrixAfter;
}

template <class T>
Matrix<T> &Matrix<T>::operator-(const Matrix &matrix)
{
    Matrix<T> *matrixAfter = new Matrix<T>(this->row, this->column);
    for (size_t i = 0; i < size; i++)
        matrixAfter[i] = data[i] - matrix.data[i];
    return matrixAfter;
}

template <class T>
Matrix<T> &Matrix<T>::operator*(const Matrix &matrix)
{
    Matrix<T> *matrixAfter = new Matrix<T>(this->row, this->column);
    for (size_t i = 0; i < size; i++)
        matrixAfter[i] = data[i] * matrix.data[i];
    return matrixAfter;
}

template <class T>
Matrix<T> operator+(const Matrix<T> &matrixLeft, const Matrix<T> &matrixRight)
{
    return matrixLeft + matrixRight;
}

template <typename T>
ostream &operator<<(ostream &os, const Matrix<T> &matrix)
{
    if (matrix.data == NULL)
    {
        cerr << "In File " << __FILE__ << ", Line " << __LINE__ << ", Function " << __FUNCTION__ << "():\nError: The data of matrix to be printed is NULL.\n";
        return os;
    }
    else
    {
        string str;
        str = "Matrix Type " + matrix.type + "\n";
        for (size_t i = 0; i < matrix.row; i++)
        {
            str += "[";
            for (size_t j = 0; j < matrix.column; j++)
            {
                str += to_string(matrix.data[i * matrix.column + j]);
                str += ",";
            }
            str += "\b]\n";
        }
        os << str;
        return os;
    }
}

template <typename T>
istream &operator>>(istream &is, Matrix<T> &matrix)
{
    char temp;
    cout << "Do you change the rows, columns and channels of the matrix?[Y/N]" << endl;
    cin >> temp;
    if (temp == 'Y')
    {
        cout << "Please enter the rows, columns, channels of the matrix(Integer):" << endl;
        int row, column, channel;
        cin >> row >> column >> channel;
        matrix.row = row;
        matrix.column = column;
        matrix.row = channel;
    }
    if (matrix.data == NULL)
    {
        cerr << "In File " << __FILE__ << ", Line " << __LINE__ << ", Function " << __FUNCTION__ << "():\nError: The data of matrix to be valued is NULL." << endl;
        return is;
    }
    else
    {
        for (size_t i = 0; i < matrix.row; i++)
        {
            for (size_t j = 0; j < matrix.column; j++)
                is >> matrix.data[i * matrix.column + j];
        }
    }
    return is;
}
#endif
