#include <iostream>
#include <cstring>
#include <cstdio>
#include <string.h>
using namespace std;

#ifndef __MATRIX_HPP__
#define __MATRIX_HPP__

#ifdef WITH_AVX2
#include <immintrin.h>
#endif

#ifdef _OPENMP
#include <omp.h>
#endif

template <typename T>
class Matrix
{
    // Friend function
    template <typename T_>
    friend ostream &operator<<(ostream &os, const Matrix<T_> &matrix);
    template <typename T_>
    friend istream &operator>>(istream &is, Matrix<T_> &matrix);
    template <typename T_>
    friend Matrix<T_> operator+(const Matrix<T_> &matrixLeft, const Matrix<T_> &matrixRight);
    template <typename T_>
    friend Matrix<T_> operator-(const Matrix<T_> &matrixLeft, const Matrix<T_> &matrixRight);
    template <typename T_>
    friend Matrix<T_> operator*(const Matrix<T_> &matrixLeft, const Matrix<T_> &matrixRight);
    template <typename T_>
    friend Matrix<T_> operator+(T_ scalar, const Matrix<T_> &matrixRight);
    template <typename T_>
    friend Matrix<T_> operator*(T_ scalar, const Matrix<T_> &matrixRight);

private:
    // Member variable
    size_t row;
    size_t column;
    size_t size;
    T *data;
    size_t *dataref;
    string type;

public:
    // Constructor and destructor
    Matrix();
    Matrix(size_t row, size_t column);
    Matrix(const Matrix &matrix);
    ~Matrix();

    size_t getRow() const;
    size_t getColumn() const;
    T *getData() const;
    string getType() const;
    size_t *getDataref() const;

    // Overloading of the operator
    Matrix<T> &operator=(const Matrix<T> &matrix);
    Matrix<T> &operator+(const Matrix<T> &matrix);
    Matrix<T> &operator-(const Matrix<T> &matrix);
    Matrix<T> &operator*(const Matrix<T> &matrix);
    Matrix<T> &operator+=(const Matrix<T> &matrix);
    Matrix<T> &operator-=(const Matrix<T> &matrix);
    Matrix<T> &operator+(T scalar);
    Matrix<T> &operator-(T scalar);
    Matrix<T> &operator*(T scalar);
    Matrix<T> &operator+=(T scalar);
    Matrix<T> &operator-=(T scalar);
    Matrix<T> &operator*=(T scalar);

    template <typename T_>
    bool operator==(const Matrix<T_> &matrix);
    template <typename T_>
    bool operator!=(const Matrix<T_> &matrix);

    // Member function
    bool value(const T *data);
    bool print();
    bool reshape(size_t row, size_t column);
    bool flatten();
    Matrix<T> *transpose();
    Matrix<T> *inverse();
    Matrix<T> *subMatrix(size_t rowStart, size_t colStart, size_t rowLength, size_t colLength);
    Matrix<T> *ROI(size_t rowStart, size_t colStart, size_t rowLength, size_t colLength);
    T determinant();

    Matrix<T> *multiplyWithScalar(T scalar);
    Matrix<T> *addWithScalar(T scalar);
    Matrix<T> *subtractWithScalar(T scalar);
};

template <typename T>
class MatrixROI : public Matrix<T>
{
    size_t rowStart;
    size_t colStart;
    size_t rowLength;
    size_t colLength;
    MatrixROI(const Matrix<T> &matrix, size_t rowStart, size_t colStart, size_t rowLength, size_t colLength)
    {
        this->rowLength = rowLength;
        this->rowStart = rowStart;
        this->colLength = colLength;
        this->colStart = colStart;
        this->row = matrix.row;
        this->column = matrix.column;
        size = this->row * this->column;
        this->data = matrix.data;
        this->dataref = matrix.dataref;
        *(this->dataref)++;
        this.type = matrix.type;
    }
    bool print()
    {
        if (this->data == NULL)
        {
            cerr << "In File " << __FILE__ << ", Line " << __LINE__ << ", Function " << __FUNCTION__ << "():\nError: The data of matrix to be printed is NULL." << endl;
            return false;
        }
        else
        {
            size_t start = rowStart * this->column + colStart;
            cout << "Matrix Type " << this->type << endl;
            for (size_t i = 0; i < rowLength; i++)
            {
                cout << "[";
                for (size_t j = 0; j < colLength; j++)
                {
                    cout << data[start + i * this->column + j] << ", ";
                }
                cout << "\b\b]" << endl;
            }
            return true;
        }
    }
};

// Get
template <typename T>
size_t Matrix<T>::getRow() const
{
    return row;
}
template <typename T>
size_t Matrix<T>::getColumn() const
{
    return column;
}
template <typename T>
T *Matrix<T>::getData() const
{
    return data;
}
template <typename T>
string Matrix<T>::getType() const
{
    return type;
}
template <typename T>
size_t *Matrix<T>::getDataref() const
{
    return dataref;
}

// Default constructor
template <typename T>
Matrix<T>::Matrix()
{
    row = column = size = 0;
    data = nullptr;
    dataref = nullptr;
    type = typeid(T).name();
}
// Constructor
template <typename T>
Matrix<T>::Matrix(size_t row, size_t column)
{
    this->row = row;
    this->column = column;
    size = this->row * this->column;
    data = new T[size * sizeof(T)]{};
    dataref = new size_t(1);
    type = typeid(T).name();
}
// Copy constructor
template <typename T>
Matrix<T>::Matrix(const Matrix<T> &matrix)
{
    row = matrix.row;
    column = matrix.column;
    size = row * column;
    data = matrix.data;
    dataref = matrix.dataref;
    *dataref++;
    type = matrix.type;
}
// Copy assignment
template <typename T>
Matrix<T> &Matrix<T>::operator=(const Matrix<T> &matrix)
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
// Destructor
template <typename T>
Matrix<T>::~Matrix()
{
    row = column = size = 0;
    if (data == NULL)
        cerr << "In File " << __FILE__ << ", Line " << __LINE__ << ", Function " << __FUNCTION__ << "():\nWarning: The data of matrix to be destructed is NULL." << endl;
    else if (--*dataref == 0)
    {
        delete[] data;
        delete dataref;
    }
}

// Value
template <typename T>
bool Matrix<T>::value(const T *data)
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
template <typename T>
bool Matrix<T>::print()
{
    if (data == NULL)
    {
        cerr << "In File " << __FILE__ << ", Line " << __LINE__ << ", Function " << __FUNCTION__ << "():\nError: The data of matrix to be printed is NULL." << endl;
        return false;
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
// Reshape
template <typename T>
bool Matrix<T>::reshape(size_t row, size_t column)
{
    if (size != row * column)
    {
        cerr << "In File " << __FILE__ << ", Line " << __LINE__ << ", Function " << __FUNCTION__ << "():\nError: The shape and the size of matrice to be reshaped do not match." << endl;
        return false;
    }
    else
    {
        this->column = column;
        this->row = row;
        return true;
    }
}
// Flatten
template <typename T>
bool Matrix<T>::flatten()
{
    row = size;
    column = 1;
    return true;
}

// Equal
template <typename T>
template <typename T_>
bool Matrix<T>::operator==(const Matrix<T_> &matrix)
{
    if (type != matrix.getType())
        return false;
    else if (row != matrix.getRow() || column != matrix.getColumn())
        return false;
    else
    {
        if (data == NULL || matrix.getData() == NULL)
        {
            cerr << "In File " << __FILE__ << ", Line " << __LINE__ << ", Function " << __FUNCTION__ << "():\nError: The data of the two matrices to be compared is NULL." << endl;
            return false;
        }
        else
            for (size_t i = 0; i < size; i++)
                if (data[i] != matrix.getData()[i])
                    return false;
        return true;
    }
}
// Not equal
template <typename T>
template <typename T_>
bool Matrix<T>::operator!=(const Matrix<T_> &matrix)
{
    if (type != matrix.getType())
        return true;
    else if (row != matrix.getRow() || column != matrix.getColumn())
        return true;
    else
    {
        if (data == NULL || matrix.getData() == NULL)
        {
            cerr << "In File " << __FILE__ << ", Line " << __LINE__ << ", Function " << __FUNCTION__ << "():\nError: The data of the two matrices to be compared is NULL." << endl;
            return false;
        }
        else
            for (size_t i = 0; i < size; i++)
                if (data[i] != matrix.getData()[i])
                    return true;
        return false;
    }
}
// Add
template <typename T>
Matrix<T> &Matrix<T>::operator+(const Matrix<T> &matrix)
{
    if (!(type == "s" || type == "i" || type == "f" || type == "d"))
    {
        cerr << "In File " << __FILE__ << ", Line " << __LINE__ << ", Function " << __FUNCTION__ << "():\nError: The type of the two matrices cannot be added." << endl;
    }
    else if (data == NULL || matrix.data == NULL)
    {
        cerr << "In File " << __FILE__ << ", Line " << __LINE__ << ", Function " << __FUNCTION__ << "():\nError: The data of the two matrices to be added is NULL." << endl;
    }
    else
    {
        if (row != matrix.row || column != matrix.column)
        {
            cerr << "In File " << __FILE__ << ", Line " << __LINE__ << ", Function " << __FUNCTION__ << "():\nError: The two matrices to be added have different sizes." << endl;
        }
        else
        {
            Matrix<T> *matrixAfter = new Matrix<T>(this->row, this->column);
            for (size_t i = 0; i < size; i++)
                matrixAfter->data[i] = data[i] + matrix.data[i];
            return *matrixAfter;
        }
    }
    return *this;
}
template <typename T>
Matrix<T> &Matrix<T>::operator+=(const Matrix<T> &matrix)
{
    if (!(type == "s" || type == "i" || type == "f" || type == "d"))
    {
        cerr << "In File " << __FILE__ << ", Line " << __LINE__ << ", Function " << __FUNCTION__ << "():\nError: The type of the two matrices cannot be added." << endl;
    }
    else if (data == NULL || matrix.data == NULL)
    {
        cerr << "In File " << __FILE__ << ", Line " << __LINE__ << ", Function " << __FUNCTION__ << "():\nError: The data of the two matrices to be added is NULL." << endl;
    }
    else
    {
        if (row != matrix.row || column != matrix.column)
        {
            cerr << "In File " << __FILE__ << ", Line " << __LINE__ << ", Function " << __FUNCTION__ << "():\nError: The two matrices to be added have different sizes." << endl;
        }
        else
        {
            for (size_t i = 0; i < size; i++)
                data[i] += matrix.data[i];
        }
    }
    return *this;
}
// Subtracted
template <typename T>
Matrix<T> &Matrix<T>::operator-(const Matrix<T> &matrix)
{
    if (!(type == "s" || type == "i" || type == "f" || type == "d"))
    {
        cerr << "In File " << __FILE__ << ", Line " << __LINE__ << ", Function " << __FUNCTION__ << "():\nError: The type of the two matrices cannot be subtracted." << endl;
    }
    else if (data == NULL || matrix.data == NULL)
    {
        cerr << "In File " << __FILE__ << ", Line " << __LINE__ << ", Function " << __FUNCTION__ << "():\nError: The data of the two matrices to be subtracted is NULL." << endl;
    }
    else
    {
        if (row != matrix.row || column != matrix.column)
        {
            cerr << "In File " << __FILE__ << ", Line " << __LINE__ << ", Function " << __FUNCTION__ << "():\nError: The two matrices to be subtracted have different sizes." << endl;
        }
        else
        {
            Matrix<T> *matrixAfter = new Matrix<T>(this->row, this->column);
            for (size_t i = 0; i < size; i++)
                matrixAfter->data[i] = data[i] - matrix.data[i];
            return *matrixAfter;
        }
    }
    return *this;
}
template <typename T>
Matrix<T> &Matrix<T>::operator-=(const Matrix<T> &matrix)
{
    if (!(type == "s" || type == "i" || type == "f" || type == "d"))
    {
        cerr << "In File " << __FILE__ << ", Line " << __LINE__ << ", Function " << __FUNCTION__ << "():\nError: The type of the two matrices cannot be subtracted." << endl;
        return NULL;
    }
    else if (data == NULL || matrix.data == NULL)
    {
        cerr << "In File " << __FILE__ << ", Line " << __LINE__ << ", Function " << __FUNCTION__ << "():\nError: The data of the two matrices to be subtracted is NULL." << endl;
        return NULL;
    }
    else
    {
        if (row != matrix.row || column != matrix.column)
        {
            cerr << "In File " << __FILE__ << ", Line " << __LINE__ << ", Function " << __FUNCTION__ << "():\nError: The two matrices to be subtracted have different sizes." << endl;
            return NULL;
        }
        else
        {
            for (size_t i = 0; i < size; i++)
                data[i] -= matrix.data[i];
            return *this;
        }
    }
}
// Multiply
template <typename T>
Matrix<T> &Matrix<T>::operator*(const Matrix<T> &matrix)
{
    if (!(type == "s" || type == "i" || type == "f" || type == "d"))
    {
        cerr << "In File " << __FILE__ << ", Line " << __LINE__ << ", Function " << __FUNCTION__ << "():\nError: The type of the two matrices cannot be multiplied." << endl;
    }
    else if (data == NULL || matrix.data == NULL)
    {
        cerr << "In File " << __FILE__ << ", Line " << __LINE__ << ", Function " << __FUNCTION__ << "():\nError: The data of the two matrices to be multiplied is NULL." << endl;
    }
    else
    {
        if (row != matrix.column || column != matrix.row)
        {
            cerr << "In File " << __FILE__ << ", Line " << __LINE__ << ", Function " << __FUNCTION__ << "():\nError: The sizes of two matrices to be multiplied do not match." << endl;
        }
        else
        {
#ifdef WITH_AVX2
            Matrix<T> *matrixA = new Matrix<T>(*this);
            Matrix<T> *matrixB = new Matrix<T>(matrix);
            Matrix<T> *matrixC = new Matrix<T>(this->row, matrix.column);
            if (this->column % 8 != 0)
            {
                matrixA = leftSupplementZero(*this);
                matrixB = rightSupplementZero(matrix);
            }
            Matrix *matrixB_trans = matrixB->transpose();
#ifdef FLOAT
            __m256 a, b;
            float sum[8] = {0.0};
            omp_set_num_threads(8);
#pragma omp parallel for
            for (size_t i = 0; i < matrixA->row; i++)
            {
                size_t matrixAIndexMinus = i * matrixA->column;
                size_t matrixCIndexMinus = i * matrixC->column;
                for (size_t j = 0; j < matrixB_trans->row; j++)
                {
                    float ans = 0.0;
                    __m256 c = _mm256_setzero_ps();
                    size_t matrixB_transIndexMinus = j * matrixB_trans->column;
                    for (size_t k = 0; k < matrixA->column; k += 8)
                    {
                        a = _mm256_loadu_ps(matrixA->data + matrixAIndexMinus + k);
                        b = _mm256_loadu_ps(matrixB_trans->data + matrixB_transIndexMinus + k);
                        c = _mm256_add_ps(c, _mm256_mul_ps(a, b));
                    }
                    _mm256_storeu_ps(sum, c);
                    for (size_t i = 0; i < 8; i++)
                        ans += sum[i];
                    matrixC->data[matrixCIndexMinus + j] = ans;
                }
            }
#endif

#ifdef INT
            __m256i a, b;
            omp_set_num_threads(8);
#pragma omp parallel for
            for (size_t i = 0; i < matrixA->row; i++)
            {
                size_t matrixAIndexMinus = i * matrixA->column;
                size_t matrixCIndexMinus = i * matrixC->column;
                for (size_t j = 0; j < matrixB_trans->row; j++)
                {
                    int ans = 0;
                    __m256i c = _mm256_setzero_si256();
                    __m256i d = _mm256_setzero_si256();
                    __m256i e = _mm256_setzero_si256();
                    size_t matrixB_transIndexMinus = j * matrixB_trans->column;
                    for (size_t k = 0; k < matrixA->column; k += 8)
                    {
                        a = _mm256_loadu_si256((__m256i *)(matrixA->data + matrixAIndexMinus + k));
                        b = _mm256_loadu_si256((__m256i *)(matrixB_trans->data + matrixB_transIndexMinus + k));
                        c = _mm256_mullo_epi32(a, b);
                        d = _mm256_add_epi32(c, d);
                    }
                    _mm256_store_si256(&e, c);
                    int *sum = (int *)&e;
                    for (size_t i = 0; i < 8; i++)
                        ans += sum[i];
                    matrixC->data[matrixCIndexMinus + j] = ans;
                }
            }
#endif
            return *matrixC;
#else
            printf("AVX2 is not supported");
#endif
        }
    }
    return *this;
}
template <typename T>
Matrix<T> *leftSupplementZero(const Matrix<T> &matrix)
{
    if (matrix.getData() == NULL)
    {
        cerr << "In File " << __FILE__ << ", Line " << __LINE__ << ", Function " << __FUNCTION__ << "():\nError: The data of the matrix to be supplemented is NULL." << endl;
        return NULL;
    }
    else
    {
        size_t remainder = matrix.getColumn() % 8;
        Matrix<T> *matrix_sup = new Matrix<T>(matrix.getRow(), matrix.getColumn() - remainder + 8);
        // matrix_sup->dataref = matrix.getDataref();
        for (size_t i = 0; i < matrix_sup->getRow(); i++)
        {
            size_t matrixIndexMinus = i * matrix.getColumn();
            size_t matrix_supIndexMinus = i * matrix_sup->getColumn();
            for (size_t j = 0; j < matrix_sup->getColumn(); j++)
            {
                if (j < matrix.getColumn())
                    matrix_sup->getData()[matrix_supIndexMinus + j] = matrix.getData()[matrixIndexMinus + j];
                else
                    matrix_sup->getData()[matrix_supIndexMinus + j] = 0.0;
            }
        }
        return matrix_sup;
    }
}
template <typename T>
Matrix<T> *rightSupplementZero(const Matrix<T> &matrix)
{
    if (matrix.getData() == NULL)
    {
        cerr << "In File " << __FILE__ << ", Line " << __LINE__ << ", Function " << __FUNCTION__ << "():\nError: The data of the matrix to be supplemented is NULL." << endl;
        return NULL;
    }
    else
    {
        size_t remainder = matrix.getRow() % 8;
        Matrix<T> *matrix_sup = new Matrix<T>(matrix.getRow() - remainder + 8, matrix.getColumn());
        // matrix_sup->dataref = matrix.getDataref();
        memcpy(matrix_sup->getData(), matrix.getData(), sizeof(T) * matrix.getRow() * matrix.getColumn());
        memset(matrix_sup->getData() + sizeof(T) * matrix.getRow() * matrix.getColumn(), 0.0, sizeof(T) * remainder * matrix.getColumn());
        return matrix_sup;
    }
}

template <typename T>
Matrix<T> &Matrix<T>::operator+(T scalar)
{
    if (!(type == "s" || type == "i" || type == "f" || type == "d"))
    {
        cerr << "In File " << __FILE__ << ", Line " << __LINE__ << ", Function " << __FUNCTION__ << "():\nError: The type of the matrix cannot be added." << endl;
    }
    else if (data == NULL)
    {
        cerr << "In File " << __FILE__ << ", Line " << __LINE__ << ", Function " << __FUNCTION__ << "():\nError: The data of the matrix to be added is NULL." << endl;
    }
    else
    {
        Matrix<T> *matrixAfter = new Matrix<T>(this->row, this->column);
        for (size_t i = 0; i < size; i++)
            matrixAfter->data[i] = data[i] + scalar;
        return *matrixAfter;
    }
    return *this;
}
template <typename T>
Matrix<T> &Matrix<T>::operator-(T scalar)
{
    if (!(type == "s" || type == "i" || type == "f" || type == "d"))
    {
        cerr << "In File " << __FILE__ << ", Line " << __LINE__ << ", Function " << __FUNCTION__ << "():\nError: The type of the matrix cannot be subtracted." << endl;
    }
    else if (data == NULL)
    {
        cerr << "In File " << __FILE__ << ", Line " << __LINE__ << ", Function " << __FUNCTION__ << "():\nError: The data of the matrix to be subtracted is NULL." << endl;
    }
    else
    {
        Matrix<T> *matrixAfter = new Matrix<T>(this->row, this->column);
        for (size_t i = 0; i < size; i++)
            matrixAfter->data[i] = data[i] - scalar;
        return *matrixAfter;
    }
    return *this;
}
template <typename T>
Matrix<T> &Matrix<T>::operator*(T scalar)
{
    if (!(type == "s" || type == "i" || type == "f" || type == "d"))
    {
        cerr << "In File " << __FILE__ << ", Line " << __LINE__ << ", Function " << __FUNCTION__ << "():\nError: The type of the matrix cannot be multiplied." << endl;
    }
    else if (data == NULL)
    {
        cerr << "In File " << __FILE__ << ", Line " << __LINE__ << ", Function " << __FUNCTION__ << "():\nError: The data of the matrix to be multiplied is NULL." << endl;
    }
    else
    {
        Matrix<T> *matrixAfter = new Matrix<T>(this->row, this->column);
        for (size_t i = 0; i < size; i++)
            matrixAfter->data[i] = data[i] * scalar;
        return *matrixAfter;
    }
    return *this;
}
template <typename T>
Matrix<T> &Matrix<T>::operator+=(T scalar)
{
    if (!(type == "s" || type == "i" || type == "f" || type == "d"))
    {
        cerr << "In File " << __FILE__ << ", Line " << __LINE__ << ", Function " << __FUNCTION__ << "():\nError: The type of the matrix cannot be added." << endl;
    }
    else if (data == NULL)
    {
        cerr << "In File " << __FILE__ << ", Line " << __LINE__ << ", Function " << __FUNCTION__ << "():\nError: The data of the matrix to be added is NULL." << endl;
    }
    else
    {
        for (size_t i = 0; i < size; i++)
            data[i] += scalar;
    }
    return *this;
}
template <typename T>
Matrix<T> &Matrix<T>::operator-=(T scalar)
{
    if (!(type == "s" || type == "i" || type == "f" || type == "d"))
    {
        cerr << "In File " << __FILE__ << ", Line " << __LINE__ << ", Function " << __FUNCTION__ << "():\nError: The type of the matrix cannot be subtracted." << endl;
    }
    else if (data == NULL)
    {
        cerr << "In File " << __FILE__ << ", Line " << __LINE__ << ", Function " << __FUNCTION__ << "():\nError: The data of the matrix to be subtracted is NULL." << endl;
    }
    else
    {
        for (size_t i = 0; i < size; i++)
            data[i] -= scalar;
    }
    return *this;
}
template <typename T>
Matrix<T> &Matrix<T>::operator*=(T scalar)
{
    if (!(type == "s" || type == "i" || type == "f" || type == "d"))
    {
        cerr << "In File " << __FILE__ << ", Line " << __LINE__ << ", Function " << __FUNCTION__ << "():\nError: The type of the matrix cannot be multiplied." << endl;
    }
    else if (data == NULL)
    {
        cerr << "In File " << __FILE__ << ", Line " << __LINE__ << ", Function " << __FUNCTION__ << "():\nError: The data of the matrix to be multiplied is NULL." << endl;
    }
    else
    {
        for (size_t i = 0; i < size; i++)
            data[i] *= scalar;
    }
    return *this;
}

// Transpose
template <typename T>
Matrix<T> *Matrix<T>::transpose()
{
    if (data == NULL)
    {
        cerr << "In File " << __FILE__ << ", Line " << __LINE__ << ", Function " << __FUNCTION__ << "():\nError: The data of matrix to be transposed is NULL." << endl;
    }
    else
    {

        Matrix<T> *matrixTrans = new Matrix<T>(column, row);
        for (size_t i = 0; i < row; i++)
        {
            size_t matrixIndexMinus = i * column;
            for (size_t j = 0; j < column; j++)
                matrixTrans->data[j * matrixTrans->column + i] = data[matrixIndexMinus + j];
        }
        return matrixTrans;
    }
    return this;
}
// SubMatrix
template <typename T>
Matrix<T> *Matrix<T>::subMatrix(size_t rowStart, size_t colStart, size_t rowLength, size_t colLength)
{
    if (data == NULL)
    {
        cerr << "In File " << __FILE__ << ", Line " << __LINE__ << ", Function " << __FUNCTION__ << "():\nError: The data of matrix to be gotten submatrix is NULL." << endl;
        return NULL;
    }
    else
    {
        if (column < colStart + colLength || row < rowStart + rowLength)
        {
            cerr << "In File " << __FILE__ << ", Line " << __LINE__ << ", Function " << __FUNCTION__ << "():\nError: The row and column of submatrix are out of range." << endl;
            return NULL;
        }
        else
        {
            Matrix<T> *matrixAfter = new Matrix<T>(rowLength, colLength);
            size_t dropLength = rowStart * column;
            for (size_t i = 0; i < rowLength; i++)
            {
                for (size_t j = 0; j < colLength; j++)
                {
                    matrixAfter->data[i * colLength + j] = data[dropLength + i * column + j + colStart];
                }
            }
            return matrixAfter;
        }
    }
}
// Inverse
template <typename T>
Matrix<T> *Matrix<T>::inverse()
{
    if (!(type == "s" || type == "i" || type == "f" || type == "d"))
    {
        cerr << "In File " << __FILE__ << ", Line " << __LINE__ << ", Function " << __FUNCTION__ << "():\nError: The type of the two matrices cannot be calculated inverse matrixd." << endl;
        return NULL;
    }
    else if (data == NULL)
    {
        cerr << "In File " << __FILE__ << ", Line " << __LINE__ << ", Function " << __FUNCTION__ << "():\nError: The data of matrix to be calculated inverse matrix is NULL." << endl;
        return NULL;
    }
    else
    {
        if (row != column)
        {
            cerr << "In File " << __FILE__ << ", Line " << __LINE__ << ", Function " << __FUNCTION__ << "():\nError: The data of matrix to be calculated inverse matrix is noninvertible matrix." << endl;
            return NULL;
        }
        else
        {
            if (column == 1)
            {
                T f = data[0];
                Matrix *matrixAfter = this->multiplyWithScalar(1. / (f * f));
                return matrixAfter;
            }
            else
            {
                Matrix *matrixAfter = new Matrix(row, column);
                for (int i = 0; i < matrixAfter->row; i++)
                {
                    for (int j = 0; j < matrixAfter->column; j++)
                    {
                        if ((i + j + 2) % 2 == 0)
                        {
                            Matrix<T> *matrix = matrixForDet(this, i, j);
                            matrixAfter->data[i * matrixAfter->column + j] = matrix->determinant();
                        }
                        else
                        {
                            Matrix<T> *matrix = matrixForDet(this, i, j);
                            T f = matrix->determinant();
                            if (f != 0.)
                                matrixAfter->data[i * matrixAfter->column + j] = -f;
                            else
                                matrixAfter->data[i * matrixAfter->column + j] = f;
                        }
                    }
                }
                return matrixAfter->transpose()->multiplyWithScalar(1. / this->determinant());
            }
        }
    }
}
// Determinant
template <typename T>
T Matrix<T>::determinant()
{
    if (!(type == "s" || type == "i" || type == "f" || type == "d"))
    {
        cerr << "In File " << __FILE__ << ", Line " << __LINE__ << ", Function " << __FUNCTION__ << "():\nError: The type of the matrice cannot be calculated determinant." << endl;
        return NULL;
    }
    else if (data == NULL)
    {
        cerr << "In File " << __FILE__ << ", Line " << __LINE__ << ", Function " << __FUNCTION__ << "():\nError: The data of the matrice to be calculated determinant is NULL." << endl;
        return NULL;
    }
    else
    {
        T result = 0.0;
        if (column == 1 && row == 1)
            result = data[0];
        else if (column == 2 && row == 2)
            result = data[0] * data[3] - data[1] * data[2];
        else
        {
            for (int i = 0; i < column; i++)
            {
                if (data[i] == 0)
                    result += 0.0;
                else
                {
                    if ((i + 2) % 2 == 0)
                    {
                        Matrix<T> *matrix = matrixForDet(this, 0, i);
                        result += data[i] * matrix->determinant();
                    }
                    else
                    {
                        Matrix<T> *matrix = matrixForDet(this, 0, i);
                        result += data[i] * (-1) * matrix->determinant();
                    }
                }
            }
        }
        return result;
    }
}
template <typename T>
Matrix<T> *matrixForDet(const Matrix<T> *matrix, size_t row, size_t column)
{
    if (matrix->getData() == NULL)
    {
        cerr << "In File " << __FILE__ << ", Line " << __LINE__ << ", Function " << __FUNCTION__ << "():\nError: The data of the matrix to be calculated the matrix for determinant is NULL." << endl;
        return NULL;
    }
    else
    {
        if (matrix->getColumn() < column || matrix->getRow() < row)
        {
            cerr << "In File " << __FILE__ << ", Line " << __LINE__ << ", Function " << __FUNCTION__ << "():\nError: The col or row is out of range." << endl;
            return NULL;
        }
        else
        {
            Matrix<T> *matrixAfter = new Matrix<T>(matrix->getRow() - 1, matrix->getColumn() - 1);
            size_t indexAns = 0;
            for (size_t i = 0; i < matrix->getRow(); i++)
            {
                if (i == row)
                    continue;
                for (size_t j = 0; j < matrix->getColumn(); j++)
                {
                    if (j == column)
                        continue;
                    else
                        matrixAfter->getData()[indexAns++] = matrix->getData()[i * matrix->getColumn() + j];
                }
            }
            return matrixAfter;
        }
    }
}

// Region of interest (ROI)
template <typename T>
Matrix<T> *Matrix<T>::ROI(size_t rowStart, size_t colStart, size_t rowLength, size_t colLength)
{
    if (data == NULL)
    {
        cerr << "In File " << __FILE__ << ", Line " << __LINE__ << ", Function " << __FUNCTION__ << "():\nError: The data of matrix to be gotten ROI is NULL." << endl;
        return NULL;
    }
    else
    {
        if (column < colStart + colLength || row < rowStart + rowLength)
        {
            cerr << "In File " << __FILE__ << ", Line " << __LINE__ << ", Function " << __FUNCTION__ << "():\nError: The row and column of ROI are out of range." << endl;
            return NULL;
        }
        else
        {
            MatrixROI<T> *matrixROI = new MatrixROI<T>(*this, rowStart, colStart, rowLength, colLength);
            return matrixROI;
        }
    }
}

// With scalar
template <typename T>
Matrix<T> *Matrix<T>::multiplyWithScalar(T scalar)
{
    if (!(type == "s" || type == "i" || type == "f" || type == "d"))
    {
        cerr << "In File " << __FILE__ << ", Line " << __LINE__ << ", Function " << __FUNCTION__ << "():\nError: The type of the two matrices cannot be multiplied with scalar." << endl;
        return NULL;
    }
    else if (data == NULL)
    {
        cerr << "In File " << __FILE__ << ", Line " << __LINE__ << ", Function " << __FUNCTION__ << "():\nError: The data of the matrice to be multiplied with scalar is NULL." << endl;
        return NULL;
    }
    else
    {
        Matrix<T> *matrixAfter = new Matrix<T>(this->row, this->column);
        for (int i = 0; i < row; i++)
        {
            for (int j = 0; j < column; j++)
                matrixAfter->data[i * column + j] = this->data[i * column + j] * scalar;
        }
        return matrixAfter;
    }
}
template <typename T>
Matrix<T> *Matrix<T>::addWithScalar(T scalar)
{
    if (!(type == "s" || type == "i" || type == "f" || type == "d"))
    {
        cerr << "In File " << __FILE__ << ", Line " << __LINE__ << ", Function " << __FUNCTION__ << "():\nError: The type of the two matrices cannot be added with scalar." << endl;
        return NULL;
    }
    else if (data == NULL)
    {
        cerr << "In File " << __FILE__ << ", Line " << __LINE__ << ", Function " << __FUNCTION__ << "():\nError: The data of the matrice to be added with scalar is NULL." << endl;
        return NULL;
    }
    else
    {
        Matrix<T> *matrixAfter = new Matrix<T>(this->row, this->column);
        for (int i = 0; i < row; i++)
        {
            for (int j = 0; j < column; j++)
                matrixAfter->data[i * column + j] = this->data[i * column + j] + scalar;
        }
        return matrixAfter;
    }
}
template <typename T>
Matrix<T> *Matrix<T>::subtractWithScalar(T scalar)
{
    if (!(type == "s" || type == "i" || type == "f" || type == "d"))
    {
        cerr << "In File " << __FILE__ << ", Line " << __LINE__ << ", Function " << __FUNCTION__ << "():\nError: The type of the two matrices cannot be subtracted with scalar." << endl;
        return NULL;
    }
    else if (data == NULL)
    {
        cerr << "In File " << __FILE__ << ", Line " << __LINE__ << ", Function " << __FUNCTION__ << "():\nError: The data of the matrice to be subtracted with scalar is NULL." << endl;
        return NULL;
    }
    else
    {
        Matrix<T> *matrixAfter = new Matrix<T>(this->row, this->column);
        for (int i = 0; i < row; i++)
        {
            for (int j = 0; j < column; j++)
                matrixAfter->data[i * column + j] = this->data[i * column + j] - scalar;
        }
        return matrixAfter;
    }
}

// Output and Input
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
                if (matrix.type=="h")
                {
                    string temp(1,matrix.data[i * matrix.column + j]);
                    str += temp;
                }
                else
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
    cout << "Will you change the rows and columns of the matrix?[Y/N]" << endl;
    cin >> temp;
    if (temp == 'Y')
    {
        cout << "Please enter the row and column of the matrix(Integer):" << endl;
        size_t row, column;
        cin >> row >> column;
        matrix.row = row;
        matrix.column = column;
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
template <typename T>
Matrix<T> operator+(const Matrix<T> &matrixLeft, const Matrix<T> &matrixRight)
{
    return matrixLeft + matrixRight;
}
template <typename T>
Matrix<T> operator-(const Matrix<T> &matrixLeft, const Matrix<T> &matrixRight)
{
    return matrixLeft - matrixRight;
}
template <typename T>
Matrix<T> operator*(const Matrix<T> &matrixLeft, const Matrix<T> &matrixRight)
{
    return matrixLeft * matrixRight;
}
template <typename T>
Matrix<T> operator+(T scalar, const Matrix<T> &matrixRight)
{
    return matrixRight + scalar;
}
template <typename T>
Matrix<T> operator*(T scalar, const Matrix<T> &matrixRight)
{
    return matrixRight * scalar;
}

#endif