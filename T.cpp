#include <iostream>
#include <cstdlib>
#include <cstring>
class MatrixAllocationError : std::exception {
};
class MatrixWrongSizeError : std::exception {
};
class MatrixIndexError : std::exception{
};

class MatrixIsDegenerateError : std::exception {
};
//=============== Matrix class ===============//
template <typename T>
class Matrix {
protected:
    int rowsCnt_;
    int colsCnt_;
    T** array_;
public:
    friend Matrix<T> operator+(const Matrix<T>& left, const Matrix<T>& right){
        if (left.rowsCnt_ != right.rowsCnt_ || left.colsCnt_ != right.colsCnt_){
            throw MatrixWrongSizeError();
        }
        Matrix<T> result(left.rowsCnt_, left.colsCnt_);
        for(int i = 0; i < left.rowsCnt_; ++i){
            for(int j = 0; j < left.colsCnt_; ++j){
                result.array_[i][j] = left.array_[i][j] + right.array_[i][j];
            }
        }
        return result;
    }
    friend Matrix<T> operator-(const Matrix<T>& left, const Matrix<T>& right) {
        if (left.rowsCnt_ != right.rowsCnt_ || left.colsCnt_ != right.colsCnt_) {
            throw MatrixWrongSizeError();
        }
        Matrix<T> result(left.rowsCnt_, left.colsCnt_);
        for(int i = 0; i < left.rowsCnt_; ++i) {
            for(int j = 0; j < left.colsCnt_; ++j){
                result.array_[i][j] = left.array_[i][j] - right.array_[i][j];
            }
        }
        return result;
    }
    friend Matrix<T> operator*(const Matrix<T>& left, const T digit) {
        Matrix<T> result(left.rowsCnt_, left.colsCnt_);
        for(int i = 0; i < left.rowsCnt_; ++i) {
            for(int j = 0; j < left.colsCnt_; ++j) {
                result.array_[i][j] = left.array_[i][j] * digit;
            }
        }
        return result;
    }
    friend Matrix<T> operator*(const Matrix<T>& left, const Matrix<T>& right) {
        if (left.colsCnt_ != right.rowsCnt_) {
            throw MatrixWrongSizeError();
        }
        Matrix<T> result(left.rowsCnt_, right.colsCnt_);
        for(int i = 0; i < result.rowsCnt_; ++i){
            for(int j = 0; j < result.colsCnt_; ++j) {
                for(int k = 0; k < left.colsCnt_; ++k) {
                    result.array_[i][j] += left.array_[i][k] * right.array_[k][j];
                }
            }
        }
        return result;
    }
    friend Matrix<T> operator*(const T digit, const Matrix<T>& right){
        Matrix<T> result(right.rowsCnt_, right.colsCnt_);
        for(int i = 0; i < right.rowsCnt_; ++i){
            for(int j = 0; j < right.colsCnt_; ++j){
                result.array_[i][j] = right.array_[i][j] * digit;
            }
        }
        return result;
    }
    friend std::istream&operator>>(std::istream& in, Matrix<T>& matrix) {
        for(int i = 0; i < matrix.rowsCnt_; ++i) {
            for(int j = 0; j < matrix.colsCnt_; ++j) {
                std::cin >> matrix.array_[i][j];
            }
        }
        return in;
    }
    friend std::ostream&operator<<(std::ostream& out, const Matrix<T>& matrix) {
        for(int i = 0; i < matrix.rowsCnt_; ++i) {
            for(int j = 0; j < matrix.colsCnt_; ++j) {
                std::cout << matrix.array_[i][j] << ' ';
                if (j == matrix.colsCnt_ - 1){
                    std::cout << std::endl;
                }
            }
        }
        return out;
    }
    Matrix<T>(const Matrix<T>& matrix) {
        rowsCnt_ = matrix.rowsCnt_;
        colsCnt_ = matrix.colsCnt_;
        array_ = new T*[rowsCnt_];
        for(int i = 0; i < rowsCnt_; ++i) {
            array_[i] = new T[colsCnt_];
        }
        for(int i = 0; i < rowsCnt_; ++i) {
            for(int j = 0; j < colsCnt_; ++j) {
                array_[i][j] = matrix.array_[i][j];
            }
        }
    }
    Matrix<T>(const int rows, const int columns) {
        rowsCnt_ = rows;
        colsCnt_ = columns;
        array_ = new T*[rowsCnt_];
        for(int i = 0; i < rowsCnt_; ++i) {
            array_[i] = new T[colsCnt_];
        }
        for(int i = 0; i < rowsCnt_; ++i) {
            for(int j = 0; j < colsCnt_; ++j){
                array_[i][j] = 0;
            }
        }
    }
    ~Matrix<T>() {
        for(int i = 0; i < rowsCnt_; ++i) {
            delete[] array_[i];
        }
        delete[] array_;
    }
    int getRowsNumber() const {
        return  rowsCnt_;
    }
    int getColumnsNumber() const {
        return colsCnt_;
    }
    Matrix<T>&operator=(const Matrix<T>& right) {
        if (&right == this) {
            return *this;
        }
        if(right.rowsCnt_ != getRowsNumber() || right.colsCnt_ != getColumnsNumber()) {
            for(int i = 0; i < getRowsNumber(); ++i) {
                delete[] array_[i];
            }
            delete[] array_;
            rowsCnt_ = right.getRowsNumber();
            colsCnt_ = right.getColumnsNumber();
            array_ = new T*[getRowsNumber()];
            for(int i = 0; i < getRowsNumber(); ++i) {
                array_[i] = new T[getColumnsNumber()];
            }
        }
        for(int i = 0; i < getRowsNumber(); ++i) {
            for(int j = 0; j < getColumnsNumber(); ++j) {
                array_[i][j] = right.array_[i][j];
            }
        }
        return *this;
    }
    Matrix<T> getTransposed() const {
        Matrix<T> result(getColumnsNumber(), getRowsNumber());
        for (int i = 0; i < getRowsNumber(); ++i) {
            for (int j = 0; j < getColumnsNumber(); ++j) {
                result.array_[j][i] = array_[i][j];
            }
        }
        return result;
    }
    Matrix<T>& transpose() {
        *this = this->getTransposed();
        return *this;
    }
    Matrix<T>&operator+=(const Matrix<T>& right) {
        Matrix<T> temp = *this;
        *this = temp + right;
        return *this;
    }
    Matrix<T>&operator-=(const Matrix<T>& right) {
        Matrix<T> temp = *this;
        *this = temp - right;
        return *this;
    }
    Matrix<T>&operator*=(const T& n) {
        for(int i = 0; i < getRowsNumber(); ++i) {
            for(int j = 0; j <getColumnsNumber(); ++j) {
                array_[i][j] *= n;
            }
        }
        return *this;
    }
    Matrix<T>&operator*=(const Matrix<T>& right) {
        Matrix<T> temp = *this;
        *this = temp * right;
        return *this;
    }
    T& operator()(const int& i, const int& j) {
        if(i >= getRowsNumber() || j >= getColumnsNumber()) {
            throw MatrixIndexError();
        }
        if (i < 0 || j < 0) {
            throw MatrixIndexError();
        }
        return array_[i][j];
    }
    T operator()(const int& i, const int& j) const {
        if(i >= getRowsNumber() || j >= getColumnsNumber()) {
            throw MatrixIndexError();
        }
        if (i < 0 || j < 0) {
            throw MatrixIndexError();
        }
        return array_[i][j];
    }
};
//=============== SquareMatrix class ===============//
template <typename T>
class SquareMatrix : public Matrix<T> {
public:
    int getSize() const {
        return this->rowsCnt_;
    }
    SquareMatrix<T>(const SquareMatrix<T>& matrix) {
        this->rowsCnt_ = matrix.rowsCnt_;
        this->colsCnt_ = this->rowsCnt_;
        for (int i = 0; i < this->rowsCnt_; ++i) {
            for (int j = 0; j < this->colsCnt_; ++j) {
                this->array_[i][j] = matrix.array_[i][j];
            }
        }
    }
    explicit SquareMatrix<T>(const int size) : Matrix<T>(size, size) {}
    ~SquareMatrix<T>() {
        for(int i = 0; i < this->rowsCnt_; ++i) {
            delete[] this->array_[i];
        }
        delete[] this->array_;
    }
    SquareMatrix<T> operator=(const SquareMatrix<T>& right) {
        static_cast<Matrix<T>>(*this) = static_cast<Matrix<T>>(right);
    }
    friend SquareMatrix<T> operator+(const SquareMatrix<T>& left, const SquareMatrix<T>& right) {
        return static_cast<SquareMatrix<T>>(static_cast<Matrix<T>>(left) + static_cast<Matrix<T>>(right));
    }
    friend SquareMatrix<T> operator-(const SquareMatrix<T>& left, const SquareMatrix<T>& right) {
        return static_cast<SquareMatrix<T>>(static_cast<Matrix<T>>(left) - static_cast<Matrix<T>>(right));
    }
    friend SquareMatrix<T> operator*(const SquareMatrix<T>& right, const T digit) {
        return static_cast<SquareMatrix<T>>(static_cast<Matrix<T>>(right) * digit);
    }
    friend SquareMatrix<T> operator*(const T digit, const SquareMatrix<T>& right) {
        return static_cast<SquareMatrix<T>>(digit * static_cast<Matrix<T>>(right));
    }
    friend SquareMatrix<T> operator*(const SquareMatrix<T>& left, const SquareMatrix<T>& right) {
        return static_cast<SquareMatrix<T>>(static_cast<Matrix<T>>(left) * static_cast<Matrix<T>>(right));
    }
    friend std::ostream&operator<<(std::ostream& out, SquareMatrix<T>& matrix) {
        std::cout << static_cast<Matrix<T>>(matrix);
    }
    friend std::istream&operator>>(std::istream& in, SquareMatrix<T>& matrix) {
        std::cin >> static_cast<Matrix<T>>(matrix);
    }
    T operator()(const int& i, const int& j) const {
        return static_cast<Matrix<T>>(*this)(i, j);
    }
    SquareMatrix<T>&operator+=(const SquareMatrix<T>& other) {
        return (*this) + other;
    }
    SquareMatrix<T>&operator-=(const SquareMatrix<T>& other) {
        return (*this) - other;
    }
    SquareMatrix<T>&operator*=(const SquareMatrix<T>& other) {
        return (*this) * other;
    }
    SquareMatrix<T>&operator*=(const T digit) {
        return  digit * (*this);
    }
    SquareMatrix<T> getTransposed() const {
        return static_cast<SquareMatrix<T>>(static_cast<Matrix<T>>(*this).getTransposed());
    }
    const T getDeterminant() const {
        SquareMatrix<T> matrix = *this;
        T determinant = 0;
        return determinant;
    }
    const SquareMatrix<T> GetInverse() const {
        if (this->Determinant() == 0) {
            throw MatrixIsDegenerateError();
        }
        int size_ = getSize();
        SquareMatrix<T> E(size_);
        SquareMatrix<T> answer = *this;
        for (int i = 0; i < size_; ++i) {
            for (int j = 0; j < size_; ++j) {
                if (i == j) {
                    E(i, j) += 1;
                }
            }
        }
        for(int i = 0; i < size_; ++i) {
            T tmp = answer.array_[i][i];
            for(int j = size_; j > 0; --j) {
                E(i, j - 1) /= tmp;
                answer(i, j - 1) /= tmp;
            }
            for(int j = 0; j < size_ ; ++j) {
                if (j != i) {
                    tmp = answer(j, i);
                    for (int k = size_; k > 0; --k) {
                        E(j, k - 1) -= E(i, k - 1) * tmp;
                        answer(j, k - 1) -= answer(i, k - 1) * tmp;
                    }
                }
            }
        }
        for (int i = 0; i < size_; ++i) {
            for (int j = 0; j < size_; ++j) {
                answer(i, j) = E(i, j);
            }
        }
    return answer;
    }
};
//================ class Rational ===============//
class BAN_Zero {};
class Rational {
    friend Rational operator+(const Rational &r1, const Rational &r2) {
        Rational result = r1;
        result += r2;
        return result;
    }
    friend Rational operator-(const Rational &r1, const Rational &r2) {
        Rational result = r1;
        result -= r2;
        return result;
    }
    friend Rational operator*(const Rational &r1, const Rational &r2) {
        Rational result = r1;
        result *= r2;
        return result;
    }
    friend Rational operator/(const Rational &r1, const Rational &r2) {
        Rational result = r1;
        result /= r2;
        return result;
    }
    friend std::ostream& operator<<(std::ostream& out, const Rational& r) {
        if (r.q == 1){
            std::cout << r.p;
        }
        else {
            std::cout << r.p << "/" << r.q;
        }
        return out;
    }
    friend std::istream& operator>>(std::istream& in, Rational& r){
        char s[15];
        std::cin >> s;
        int l = strlen(s);
        int i = 0;
        for (i = 0; i < l; ++i){
            if (s[i] == '/') break;
        }
        if (i == l) r = Rational(atoi(s));
        else r = Rational(atoi(s), atoi(s + i + 1));
        return in;
    }
    friend bool operator>(const Rational& r1, const Rational& r2) {
        return (r1.p * r2.q > r2.p * r1.q);
    }
    friend bool operator==(const Rational& r1, const Rational& r2) {
        return (r1.p * r2.q == r2.p * r1.q);
    }
    friend bool operator>=(const Rational& r1, const Rational& r2) {
        return (r1.p * r2.q >= r2.p * r1.q);
    }
    friend bool operator<(const Rational& r1, const Rational& r2) {
        return (r2 > r1);
    }
    friend bool operator!=(const Rational& r1, const Rational& r2) {
        return (r1.p * r2.q != r2.p * r1.q);
    }
    friend bool operator<=(const Rational& r1, const Rational& r2) {
        return (r1.p * r2.q <= r2.p * r1.q);
    }
private:
    int p;
    int q;
    void reduce() {
        int n = GCD(abs(p), abs(q));
        q /= n;
        p /= n;
        if (q < 0){
            q *= -1;
            p *= -1;
        }
        if (p == 0) q = 1;
    }
    int GCD(int q, int p) {
        if (q == 0) return p;
        if (p == 0) return q;
        if (p > q) return GCD(q, p % q);
        else return GCD(q % p, p);
    }
public:
    Rational(int _p = 0, int _q = 1): p(_p), q(_q) {
        if (q < 0) {
            q *= -1;
            p *= -1;
        }
        this->reduce();
    }

    int getNumerator() const {
        return p;
    }

    int getDenuminator() const {
        return q;
    }

    Rational& operator+=(const Rational &r) {
        p = p * r.q + r.p * q;
        q *= r.q;
        this->reduce();
        return *this;
    }
    Rational& operator-=(const Rational &r) {
        p = p * r.q - r.p * q;
        q *= r.q;
        this->reduce();
        return *this;
    }
    Rational operator-() const {
        return Rational(-p, q);
    }
    Rational& operator*=(const Rational &r) {
        p *= r.p;
        q *= r.q;
        this->reduce();
        return *this;
    }
    Rational& operator/=(const Rational &r) {
        if (r.p == 0){
            throw BAN_Zero();
            return *this;
        }
        p *= r.q;
        q *= r.p;
        this->reduce();
        return *this;
    }
    Rational& operator++() {
        (*this) += 1;
        return *this;
    }
    Rational operator++(int fictive) {
        Rational result = *this;
        (*this) += 1;
        return result;
    }
    Rational& operator--() {
        (*this) -= 1;
        return *this;
    }
    Rational operator--(int fictive) {
        Rational result = *this;
        (*this) -= 1;
        return result;
    }
    Rational operator+() const {
        return *this;
    }
    Rational& operator=(const Rational& r) {
        p = r.p;
        q = r.q;
        return *this;
    }
};
//=================== main() ===============//
using namespace std;
int main() {
    /*
    int m, n, p;
    Rational r;
    cin >> m >> n >> p >> r;

    Matrix<Rational> A(m, n);
    SquareMatrix<Rational> S(p);
    cin >> A >> S;

    try {
        cout << (A * S) * A.getTransposed() << endl;
    } catch (const MatrixWrongSizeError&) {
        cout << "A and S have not appropriate sizes for multiplication." << endl;
    }

    cout << (r * (S = S) * S).getSize() << endl;

    SquareMatrix<Rational> P(S);

    cout << (P * (S + S - 3 * P)).getDeterminant() << endl;

    const SquareMatrix<Rational>& rS = S;

    cout << rS.getSize() << ' ' << rS.getDeterminant() << ' ' << rS.getTrace() << endl;
    cout << (S = S) * (S + rS) << endl;
    cout << (S *= S) << endl;

    try {
        cout << rS.getInverse() << endl;
        cout << P.invert().getTransposed().getDeterminant() << endl;
        cout << P << endl;
    } catch (const MatrixIsDegenerateError&) {
        cout << "Cannot inverse matrix." << endl;
    }
*/
    return 0;
}
