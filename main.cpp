#include <iostream>
#include <vector>
#include <math.h>

using namespace std;

#define matrix_exception "Error: the dimensional problem occurred\n"
#define matrix_singular "Error: matrix A is singular\n"
#define COUT_MATRIX_OP(x) try {cout << x;} catch (const char *e) {cout << e;}

class Matrix {
protected:
    int n, m;
    vector<vector<double>> v;
public:
    Matrix ( int n, int m ) : n{n}, m{m} {
        this->v.resize(n);
        for ( auto &x : v ) { x.resize(m); };
    }

    vector<vector<double>> get_matrix() {
        return this->v;
    }

    double get_value(int row, int column) {
        return this->v[row-1][column-1];
    }

    void set_value(int row, int column, double value) {
        this->v[row - 1][column - 1] = value;
    }

    int get_rows() {
        return this->n;
    }
    int get_columns() {
        return this->m;
    }

    auto begin()  {
        return this->v.begin();
    }

    auto end() {
        return this->v.end();
    }


    Matrix operator+(const Matrix &A) {
        if ( this->n != A.n || this-> m != A.m ) throw matrix_exception;
        Matrix tmp(n, m);
        double i = 0;
        for ( auto &row : tmp ) {
            double j = 0;
            for ( auto &e : row ) {
                e = this->v[i][j] + A.v[i][j];
                ++j;
            }
            ++i;
        }
        return tmp;
    }

    Matrix operator-(const Matrix &A) {
        if ( this->n != A.n || this-> m != A.m ) throw matrix_exception;
        Matrix tmp(n, m);
        int i = 0;
        for ( auto &row : tmp ) {
            double j = 0;
            for ( auto &e : row ) {
                e = this->v[i][j] - A.v[i][j];
                ++j;
            }
            ++i;
        }
        return tmp;
    }

    Matrix operator*(Matrix A) {
        if ( this->m != A.n ) throw matrix_exception;
        Matrix tmp(this->n, A.m);
        for ( int i = 0; i < this->n; ++i ) {
            for ( int j = 0; j < A.m; ++j) {
                for ( int z = 0; z < this->m; ++z ) {
                    tmp.v[i][j] += this->v[i][z] * A.v[z][j];
                }
            }
        }
        return tmp;
    }

    Matrix T() {
        Matrix tmp(m, n);
        double i = 0;
        for ( auto &row : tmp ) {
            double j = 0;
            for ( auto &e : row ) {
                e = this->v[j][i];
                ++j;
            }
            ++i;
        }
        return tmp;
    }
};



// OUTPUT
void operator<<(basic_ostream<char> &out, Matrix A) {
    for ( const auto &row : A ) {
        for ( const auto &e : row ) {
            printf("%.4f", e);
            cout << (&row.back() == &e ? "" : " ") ;
        }
        out << '\n';
    }
}

// INPUT
void operator>>(basic_istream<char> &in, Matrix &A) {
    for ( auto &row : A ) {
        for ( auto &e : row ) {
            in >> e;
        }
    }
}


class SquareMatrix : public Matrix {
public:
    SquareMatrix(int n) : Matrix(n, n) {};

    SquareMatrix& operator=(const Matrix& other) {
        Matrix::operator=(other);
        return *this;
    }

    double determinant() const {
        SquareMatrix temp(*this);
        double det = 1.0;
        int step = 0;
        for (int i = 0; i < temp.get_rows(); ++i) {
            int pivot_row = i;
            for (int j = i + 1; j < temp.get_rows(); ++j) {
                if (abs(temp.v[j][i]) > abs(temp.v[pivot_row][i])) {
                    pivot_row = j;
                }
            }

            if (pivot_row != i) {
                swap(temp.v[pivot_row], temp.v[i]);
//                cout << "step #" << ++step << ": permutation\n";
//                cout << temp;
                det *= -1;
            }

            if (temp.v[i][i] == 0) {
                return 0;
            }

            for (int j = i + 1; j < temp.get_rows(); ++j) {
                double factor = temp.v[j][i] / temp.v[i][i];
                for (int k = i; k < temp.get_rows(); ++k) {
                    temp.v[j][k] -= factor * temp.v[i][k];
                }
                if ( factor != 0 ) {
//                    cout << "step #" << ++step << ": elimination\n";
//                    cout << temp;
                }
            }
            det *= temp.v[i][i];
        }

        return det;
    }

};

class IdentityMatrix : public SquareMatrix {
public:
    explicit IdentityMatrix(int n) : SquareMatrix(n) {
        for ( double i = 0; i < n; ++i ) {
            this->v[i][i] = 1;
        }
    };

    SquareMatrix& operator=(const Matrix& other) {
        Matrix::operator=(other);
        return *this;
    }
};

class EliminationMatrix : public IdentityMatrix {
protected:
    double r;
    double c;
public:
    EliminationMatrix(int r, int c, SquareMatrix &m) : IdentityMatrix(m.get_rows()), r(r), c(c) {
        for ( int i = c; i <= m.get_rows(); ++i ) {
            if ( m.get_value(i, c) != 0 and i != r ) {
                this->set_value(r,i,-m.get_value(r, c) / m.get_value(i,c));
                break;
            }
        }
    }
};

class PermutationMatrix : public IdentityMatrix {
protected:
    double r1, r2;
public:
    PermutationMatrix(SquareMatrix &m, double r1, double r2) : IdentityMatrix(m.get_rows()), r1(r1), r2(r2) {
        r1--, r2--;
        this->v[r1][r1] = 0;
        this->v[r2][r2] = 0;
        this->v[r1][r2] = 1;
        this->v[r2][r1] = 1;
    }
};


SquareMatrix getInverse(SquareMatrix &m) {
    if ( m.determinant() == 0 ) {
        throw matrix_singular;
    }
    int n = m.get_rows();
    Matrix augmented(n, n*2);
    SquareMatrix temp(m);
    IdentityMatrix I(n);
    for ( int r = 1; r <= n; r++ ) {
        for ( int c = 1; c <= n; c++ ) {
            augmented.set_value(r, c, m.get_value(r, c));
            augmented.set_value(r, c+n, I.get_value(r, c));
        }
    }
//    cout << "Augmented matrix:\n" << augmented;
//    cout << "Gaussian process:\n";
    int step = 0;
    /*
     * Upper triangle
     */
    for (int i = 1; i <= n; ++i) {
        int pivot_row = i;
        for (int j = i + 1; j <= n; ++j) {
            if (abs(augmented.get_value(j, i)) > abs(augmented.get_value(pivot_row, i))) {
                pivot_row = j;
            }
        }

        if (pivot_row != i) {
            PermutationMatrix Ppi(temp, pivot_row, i);
            temp = Ppi * temp;
            I = Ppi * I;
            augmented = (Ppi * augmented);
//            cout << "step #" << ++step << ": permutation\n";
//            cout << augmented;
        }

        if (augmented.get_value(i,i) == 0) {
            return 0;
        }

        for (int j = i + 1; j <= augmented.get_rows(); ++j) {
            if ( augmented.get_value(j, i) != 0 ) {
                EliminationMatrix Eij(j, i, temp);
                temp = Eij * temp;
                I = Eij * I;
                augmented = Eij * augmented;
//                cout << "step #" << ++step << ": elimination\n";
//                cout << augmented;
            }
        }
    }

    /*
     * Diagonal form
     */
    for (int c = augmented.get_rows(); c >= 1; --c ) { // for each column
        for (int i = c - 1; i >= 1; --i ) {
            if (augmented.get_value(i, c) != 0 ) {
                EliminationMatrix Eij(i, c, temp);
                temp = Eij * temp;
                I = Eij * I;
                augmented = Eij * augmented;
//                cout << "step #" << ++step << ": elimination\n";
//                cout << augmented;
            }
        }
    }
    /*
     * Diagonal normalization
     */
    IdentityMatrix dn(n);
    for (int i = 1; i <= n; ++i ) {
        dn.set_value(i,i,1/temp.get_value(i,i));
    }
    temp = dn * temp;
    I = dn * I;
    augmented = dn * augmented;
//    cout << "Diagonal normalization:\n" << augmented;
    return I;
}

class ColumnVector : public Matrix {
public:
    ColumnVector(int n) : Matrix(n, 1) {};

    ColumnVector& operator=(const Matrix& other) {
        Matrix::operator=(other);
        return *this;
    }

    Matrix operator*(Matrix &A) {
        try {
            return Matrix::operator*(A);
        } catch(const char *e) {
            return this->T() * A;
        }
    }
};

ColumnVector solvingSLAE(SquareMatrix &A, ColumnVector b) {
    if ( A.determinant() == 0 ) {
        throw matrix_singular;
    }
    int n = A.get_rows();
    SquareMatrix temp(A);
    cout << "Gaussian process:\n";
    int step = 0;
    /*
     * Upper triangle
     */
    for (int i = 1; i <= n; ++i) {
        int pivot_row = i;
        for (int j = i + 1; j <= n; ++j) {
            if (abs(temp.get_value(j, i)) > abs(temp.get_value(pivot_row, i))) {
                pivot_row = j;
            }
        }

        if (pivot_row != i) {
            PermutationMatrix Ppi(temp, pivot_row, i);
            temp = Ppi * temp;
            b = Ppi * b;
            cout << "step #" << ++step << ": permutation\n";
            cout << temp;
            cout << b.T();
        }

        if (temp.get_value(i,i) == 0) {
            return 0;
        }

        for (int j = i + 1; j <= n; ++j) {
            if ( temp.get_value(j, i) != 0 ) {
                EliminationMatrix Eij(j, i, temp);
                temp = Eij * temp;
                b = Eij * b;
                cout << "step #" << ++step << ": elimination\n";
                cout << temp;
                cout << b.T();
            }
        }
    }

    /*
     * Diagonal form
     */
    for (int c = n; c >= 1; --c ) { // for each column
        for (int i = c - 1; i >= 1; --i ) {
            if (temp.get_value(i, c) != 0 ) {
                EliminationMatrix Eij(i, c, temp);
                temp = Eij * temp;
                b = Eij * b;
                cout << "step #" << ++step << ": elimination\n";
                cout << temp;
                cout << b.T();
            }
        }
    }
    /*
     * Diagonal normalization
     */
    IdentityMatrix dn(n);
    for (int i = 1; i <= n; ++i ) {
        dn.set_value(i,i,1/temp.get_value(i,i));
    }
    temp = dn * temp;
    b = dn * b;
    cout << "Diagonal normalization:\n" << temp;
    cout << b.T();
    return b;
}


int main() {
    // Input format from "PA. Task #7. Least square approximation"
//    int m;
//    cin >> m;
//    vector<double> t(m), b(m);
//    for ( int i = 0; i < m; ++i ) {
//        cin >> t[i] >> b[i];
//    }
//    int n;
//    cin >> n;
    int m = 16;
    vector<double> t = {  1,   2,   3,   4,   5,   6,   7,   9,  11,   15,  18,  19,  21,  26,  30,  43};
    vector<double> b = {0.1, 0.2, 0.5, 1.1, 1.7, 2.5, 3.4, 5.0, 6.5, 8.2, 7.6, 7.4, 5.9, 4.2, 3.1, 2.1};
    int n = 4; // the degree of the approximation polynomial n

    Matrix A(m, n+1);
    for ( int i = 1; i <= m; ++i ) {
        for ( int j = 1; j <= n+1; j++ ) {
            A.set_value(i, j, pow(t[i-1], j-1));
        }
    }
    cout << "A:\n" << A;
    Matrix ata = A.T() * A;
    cout << "A_T*A:\n" << ata;
    auto* ata_sq = (SquareMatrix *) &ata;
    auto ata_sq_inv = getInverse(*ata_sq);
    cout << "(A_T*A)^-1:\n" << ata_sq_inv;
    Matrix b_matrix(m, 1);
    for ( int i = 0; i < m; ++i ) {
        b_matrix.set_value(i+1, 1, b[i]);
    }
    cout << "A_T*b:\n" << A.T() * b_matrix;
    Matrix x = (*(Matrix *)(&ata_sq_inv)) * (A.T() * b_matrix);
    cout << "x~:\n" << x;

    string equation;
    for ( int power = n; power >= 0; power-- ) {
        double coef = x.get_value(power+1, 1);
        if ( coef >= 0 ) equation += "+ ";
        equation += to_string(coef) + " * x**" + to_string(power) + " ";
    }
    equation.pop_back();

    cout << equation;

    FILE* pipe = popen("gnuplot -persist", "w");

    if ( pipe != NULL ) {
        fprintf(pipe, "plot %s title 'line' with lines, '-' with points pointtype 7 pointsize 1 title 'point'\n", equation.c_str());

        for (int i = 0; i < m; i++) {
            fprintf(pipe, "%f\t%f\n", t[i], b[i]);
        }

        fprintf(pipe, "e\n");
        fflush(pipe);
        pclose(pipe);
    } else {
        cout << "gnuplot problem";
    }


    return 0;
}
