#include <stdio.h>
#include <stdlib.h>
#include <vector>
#include <cmath>
#include <climits>
#include <algorithm>
#include <limits>
#include <iostream>
using namespace std;



void release_memory(int       m,
                    double ** A,
                    double *  b,
                    double *  c) {
    if (A) {
        int i;
        for (i = 0; i < m; i++) {
            if (A[i]) {
                free(A[i]);
            } else {
                /* If A[i] is NULL, then the next entries won't have been
                   allocated and their pointer values will be undefined. */
                break;
            }
        }
        free (A);
    }

    if (b) {
        free(b);
    }
    if (c) {
        free(c);
    }
}

/* The function reads an LP instance from filename. The file
   format is expected to be exactly as in the problem specification.
   On return, *m and *n will be the number of rows and columns,
   respectively; A, b and c will be the the matrix, right-hand side
   vector and objective function vector, respectively. Note that
   the order of the data in the input file is c, then b, then A.
   Memory will be allocated for A, b and c; it is the caller's
   responsibility to release the memory later.
*/
int read_LP(const char * filename,
            int *        m,
            int *        n,
            double ***   A,
            double **    b,
            double **    c)
{
    FILE *fp;

    if (!(fp = fopen(filename, "r")))
    {
        fprintf(stderr, "Could not open \"%s\".\n", filename);
        return EXIT_FAILURE;
    } else {
        int i;
        int j;

        fscanf(fp, "%d %d\n", m, n);

        /* Memory allocation. */
        if ((*A = (double**) malloc(*m * sizeof(double*))) &&
            (*b = (double*)  malloc(*m * sizeof(double)))  &&
            (*c = (double*)  malloc(*n * sizeof(double)))) {
            for (i = 0; i < *m; i++) {
                if (!((*A)[i] = (double*) malloc(*n * sizeof(double)))) {
                    fprintf(stderr, "Memory allocation failure.\n");
                    release_memory(i - 1, *A, *b, *c);
                    fclose(fp);
                    return EXIT_FAILURE;
                }
            }
        } else {
            fprintf(stderr, "Memory allocation failure.\n");
            release_memory(0, *A, *b, *c);
            fclose(fp);
            return EXIT_FAILURE;
        }

        /* Copying the values into A, b and c. */
        for (j = 0; j < *n; j++) {
            fscanf(fp, "%lf", *c + j);
        }
        for (i = 0; i < *m; i++) {
            fscanf(fp, "%lf", *b + i);
        }
        for (i = 0; i < *m; i++) {
            for (j = 0; j < *n; j++) {
                /* fscanf(fp, "%lf", (*A)[i]+j); */
                fscanf(fp, "%lf", *(*A+i)+j);
            }
        }

        fclose(fp);
        return EXIT_SUCCESS;
    }
}


void test_it(const char * lp_file)
{
    int       m = 0;
    int       n = 0;
    double ** A;
    double *  b;
    double *  c;

    if (read_LP(lp_file, &m, &n, &A, &b, &c) != EXIT_SUCCESS) {

        fprintf(stderr, "Error parsing \"%s\".\n", lp_file);

    } else {

        int i;
        int j;

        printf("A has %d row%s and %d column%s.\n",
               m, (m == 1 ? "" : "s"),
               n, (n == 1 ? "" : "s"));
        printf("c = [");
        for (j = 0; j < n; j++) {
            printf("%.1lf%s", c[j], (j == n-1 ? "]\n" : ", "));
        }
        printf("transpose(b) = [");
        for (i = 0; i < m; i++) {
            printf("%.1lf%s", b[i], (i == m-1 ? "]\n" : ", "));
        }
        for (i = 0; i < m; i++) {
            printf("A[%d] = [", i);
            for (j = 0; j < n; j++) {
                printf("%5.1lf%s", A[i][j], (j == n-1 ? "]\n" : ", "));
            }
        }

        release_memory(m, A, b, c);
    }
}


template <class T>
void swapElements(vector<T> *x, int indexA, int indexB) {
    T temp = (*x)[indexA];
    (*x)[indexA] = (*x)[indexB];
    (*x)[indexB] = temp;
}

void printVector(vector<double> x) {
    for (int i = 0; i < x.size(); i++) {
        cout << x[i] << ' ';
    }
    cout << endl;
}

void printMatrix(vector<vector<double> > A) {
    for (int i = 0; i < A.size(); i++) {
        printVector(A[i]);
    }
}

//find one possible solution x for Ax = b
//assume A is quadratic
vector<double> solve(vector<vector<double> > A, vector<double> b) {
    int n = A.size();

    vector<double> x(n);

    printMatrix(A);
    cout << "*x=" << endl;
    printVector(b);
    cout << "=============================" << endl;

    //bring into upper triangle form

    for (int i = 0; i < n - 1; i++) {
        printMatrix(A);
        cout << "*x=" << endl;
        printVector(b);
        cout << "=============================" << endl;
        //find pivot
        int maxIndex = i;
        for (int j = i; j < n; j++) {
            if (A[j][i] > abs(A[maxIndex][i])) {
                maxIndex = j;
            }
        }
        cout << "max index: " << maxIndex << endl;
        swapElements(&A, maxIndex, i);
        swapElements(&b, maxIndex, i);

        printMatrix(A);
        cout << "*x=" << endl;
        printVector(b);
        cout << "=============================" << endl;

        for (int j = i+1; j < n; j++) {
            double rowFactor = A[j][i]/A[i][i];
            for (int k = i+1; k < n; k++) {
                A[j][k] -= A[i][k]*rowFactor;
            }
            A[j][i] = 0;

            b[j] -= b[i]*rowFactor;
        }
    }

    cout << "Upper triangle done!" << endl;
    printMatrix(A);
    cout << "*x=" << endl;
    printVector(b);
    cout << "=============================" << endl;

    //backsubstitution
    for (int i = n-1; i >= 0; i--) {
        x[i] = b[i];
        for(int j = i+1; j < n; j++) {
            x[i] -= A[i][j]*x[j];
        }

        x[i] /= A[i][i];
    }

    return x;
}

//return A_B or x_B or c_B, etc.
template <class T>
vector<T> indexMatrix(vector<T> x, vector<int> B) {
    vector<T> indexedVector(B.size());
    for (int i = 0; i < B.size(); i++) {
        indexedVector[i] = x[B[i]];
    }

    return indexedVector;
}

vector<vector<double> > transpose(vector<vector<double > > A) {
    vector<vector<double> > ATranspose(A[0].size());
    for (int i = 0; i < A[0].size();i++) {
        vector<double> row(A.size());
        for (int j = 0; j < A.size(); j++) {
            row[j] = A[j][i];
        }
        ATranspose[i] = row;
    }

    return ATranspose;
}

double dotProduct(vector<double> x, vector<double> y) {
    double sum = 0;
    for (int i = 0; i < x.size(); i++) {
        sum += x[i]*y[i];
    }
}

vector<double> matrixVectorMultiplication(vector<vector<double> > A, vector<double> x) {
    vector<double> result(A[0].size());
    for (int i = 0; i < A.size(); i++) {
        result[i] = dotProduct(A[i], x);
    }

    return result;
}

//returns x + factor*y
vector<double> add(vector<double> x, vector<double> y, double factor) {
    vector<double> result(x.size());
    for (int i = 0; i < x.size(); i++) {
        result[i] = x[i] - y[i];
    }

    return result;
}

void calcN(vector<int> B, vector<int>* N, int n) {
    int j = 0;
    for (int i = 0; i < n; i++) {
        //B not contains i
        if(find(B.begin(), B.end(), i) == B.end()) {
            (*N)[j] = i;
            j++;
        }
    }

}

//checks if x <= a componentwise
bool lessThan(vector<double> x, double a) {
    for (int i = 0; i < x.size(); i++) {
        if (x[i] > a) {
            return false;
        }
    }
    return true;
}

void removeValue(vector<double> *x, double value) {
    int firstIndex;
    for (int i = 0; i < (*x).size(); i++) {
        if ((*x)[i] == value) {
            firstIndex = i;
            break;
        }
    }
    (*x).erase((*x).begin()+firstIndex);
}



vector<double> simplex(vector<vector<double> > A, vector<double> b, vector<double> c, vector<int> B) {
    //calc extra A transpose once.
    vector<vector<double> > ATranspose = transpose(A);

    int m = A.size();
    int n = A[0].size();

    int STILLSEARCHING = 0;
    int OPTIMAL = 1;
    int UNBOUNDED = -1;

    int result = STILLSEARCHING;

    vector<int> N(n-m);
    calcN(B,&N,n);

    vector<double> x(n);
    vector<double> x_B = solve(A_B)

    while(result == STILLSEARCHING) {
        //maxIndex for update
        int j,gamma,i;

        //BTRAN
        vector<vector<double> > AT_B = indexMatrix(ATranspose, B);
        vector<vector<double> > A_B = transpose(AT_B);
        vector<double> c_B = indexMatrix(c,B);
        vector<double> y = solve(AT_B, c_B);
        //PRICING
        vector<double> c_N = indexMatrix(c,N);
        vector<vector<double> > AT_N = indexMatrix(ATranspose, N);
        vector<double> z = add(c_N, matrixVectorMultiplication(AT_N, y), -1);
        if(lessThan(z, 0)) {
            //B is optimal
            result = OPTIMAL;
            break;
        } else {
            //Dantzig's Rule
            int maxIndex = -1;
            double maxValue = -numeric_limits<double>::infinity();
            for (int k = 0; k < z.size(); k++) {
                if (z[k] > maxValue) {
                    maxValue = z[k];
                    maxIndex = N[k];
                }
            }
            j = maxIndex;
        }
        //FTRAN
        vector<double> jthColumnOfA = ATranspose[j];
        vector<double> w = solve(A_B, jthColumnOfA);
        //RATIOTEST
        if(lessThan(w,0)) {
            result == UNBOUNDED;
            break;
        } else {
            double min = INT_MAX;
            //noch ohne rule. erstmal nur finde erstes. ist wahrscheinlich bullshit
            for(int k = 0; k < m; k++) {
                if (w[k] > 0 && x[B[k]]/w[k] < min) {
                    i = k;
                    min = x[B[k]]/w[k];
                }
            }
            gamma = min;


        }
        //UPDATE
        for (int k = 0; k < B.size(); k++) {
            x[B[k]] -= gamma*w[k];
        }
        x[j] = gamma;
        //N := N \{j} union {B_i}
        for (int k = 0; k < N.size(); k++) {
            if (N[k] == j) {
                N[k] = B[i];
            }
        }
        B[i] = j;
    }


}

int main(int argc, const char * argv[]) {
    vector<double> row1(6);
    row1[0] = 1;
    row1[1] = 4;
    row1[2] = 7;
    row1[3] = 1;
    row1[4] = 0;
    row1[5] = 0;

    vector<double> row2(6);
    row2[0] = -2;
    row2[1] = 6;
    row2[2] = 8;
    row2[3] = 0;
    row2[4] = 1;
    row2[5] = 0;

    vector<double> row3(6);
    row3[0] = 8;
    row3[1] = -6;
    row3[2] = 10;
    row3[3] = 0;
    row3[4] = 0;
    row3[5] = 1;

    vector<vector<double> > A(3);

    A[0] = row1;
    A[1] = row2;
    A[2] = row3;

    vector<double> b(3);
    b[0] = -5;
    b[1] = -12;
    b[2] = 6;

    vector<double> x = solve(A,b);

    printVector(x);

    //cout << x[0] << " " << x[1] << " " << " " << x[2] << endl;

    vector<double> c = A[0];
    c[3] = 0;

    vector<int> B(3);
    B[0] = 3;
    B[1] = 4;
    B[2] = 5;

    simplex(A, b, c, B);

    return EXIT_SUCCESS;
}