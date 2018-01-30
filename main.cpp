#include <stdio.h>
#include <stdlib.h>
#include <vector>
#include <cmath>
#include <climits>
#include <algorithm>
#include <limits>
#include <iostream>
using namespace std;

//typedef vector<vector<double> > RowIndexedMatrix;
typedef vector<vector<double> > ColumnIndexedMatrix;

const int STILLSEARCHING = 0;
const int OPTIMAL = 1;
const int UNBOUNDED = -1;

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
                }
            }
        } else {
            fprintf(stderr, "Memory allocation failure.\n");
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

       // release_memory(m, A, b, c);
    }
}

//Speichert 2D-array in Spaltenindizierte Matrix

vector<vector<double> > array_to_vector_M(double **A, int m,int n){

    vector<vector<double> > M(n);
    //cout « M.size() « endl;

    for(int i=0;i<n;i++){
        vector<double> col(m);

        for (int j=0;j<m;j++){
            col[j] = A[j][i];
        }
        M[i] = col;
    }
    return M;
}

vector<double> array_to_vector(double *c, int n) {
    vector<double> c_vec(n);

    for (int i = 0; i < n; i++) {
        c_vec[i] = c[i];
    }

    return c_vec;
}
//Skalar mit Vektor multiplizieren
vector<double> scalar_vec(vector<double> v,double scalar){
    vector<double> result(v.size());
    for(int i=0;i<v.size();i++){
        result[i] = v[i]*scalar;
    }
    return result;
}


//Variablen splitten für Matrix
vector<vector<double> > split_variables(vector<vector<double> > M){

    vector<vector<double> > M_new(2*M.size());

    for(int i=0;i<M.size();i++){
        M_new[2*i] = M[i];
        M_new[2*i+1] = scalar_vec(M[i],-1.0);
    }
    return M_new;
}
//Variablen splitten für Vektor
vector<double> split_variables(vector<double> c){

    vector<double> c_new(2*c.size());

    for(int i=0;i<c.size();i++){
        c_new[2*i] = c[i];
        c_new[2*i+1] = -c[i];
    }
    return c_new;
}

//liefert den i-ten dim-dimensionalen Einheitsvektor
vector<double> unit_vector(int dim, int i){

    vector<double> u(dim);
    u[i] = 1;
    return u;
}

//Slackvariablen in Matrix einführen (Identitätsmatrix hinten dranklatschen)
vector<vector<double> > slack(vector<vector<double> > M){

    for(int i=0;i<M[0].size();i++){
        M.push_back(unit_vector(M[0].size(),i));
    }

    return M;
}
//Einträge der Slackvariablen in c auf 0 setzen
vector<double>  slack(vector<double> c, int n){

    for(int i=0;i<n;i++){
        c.push_back(0);
    }
    return c;


}

void get_standard(vector<vector<double> > *A, vector<double> *c, int n){

    *A = split_variables(*A);
    *c = split_variables(*c);
    *A = slack(*A);
    *c = slack(*c,(*A)[0].size());

}


//Liefert Matrix des Hilfs-LPs zur Bestimmung der initialen Basis
vector<vector<double> > aux_lp(vector<vector<double> > A, vector<double> b, vector<double> * aux_b){

    vector<vector<double> > M(A.size());
    //ZEILEN MIT b<0 NEGIEREN
    for(int j=0;j<A.size();j++){
        vector<double> col;
        for(int i=0;i<A[0].size();i++){

            if(b[i]>=0){
                col.push_back(A[j][i]);
                (*aux_b)[i] = b[i];
            }
            else{
                col.push_back(-A[j][i]);
                (*aux_b)[i] = - b[i];
            }
        }
        M[j] = col;
    }

    for(int i=0;i<M[0].size();i++){
        M.push_back(unit_vector(M[0].size(),i));
    }

    return M;
}

//Liefert c für das Hilfs-LP
vector<double> aux_lp(vector<double> c, int m){
    vector<double> c_new(c.size());

    for(int i=0;i<m;i++){
        c_new.push_back(-1);
    }
    return c_new;
}

void swapElements(vector<double> *x, int indexA, int indexB) {
    double temp = (*x)[indexA];
    (*x)[indexA] = (*x)[indexB];
    (*x)[indexB] = temp;
}

void swapRows(ColumnIndexedMatrix *A, int rowA, int rowB, vector<int> indexSet) {
    vector<double> tempRowA(indexSet.size());
    for (int i = 0; i < indexSet.size(); i++) {
        tempRowA[i] = (*A)[indexSet[i]][rowA];
    }

    for (int i = 0; i < indexSet.size(); i++) {
        (*A)[indexSet[i]][rowA] = (*A)[indexSet[i]][rowB];
        (*A)[indexSet[i]][rowB] = tempRowA[i];
    }
}

template<typename T>
void printVector(vector<T> x) {
    for (int i = 0; i < x.size(); i++) {
        cout << x[i] << ' ';
    }
    cout << endl;
}

void printMatrix(ColumnIndexedMatrix A) {
    for (int i = 0; i < A[0].size(); i++) {
        for (int j = 0; j < A.size(); j++) {
            cout << A[j][i] << ' ';
        }
        cout << endl;
    }
}

void transformIntoUpperTriangleMatrix(ColumnIndexedMatrix *A, vector<double> *b, vector<int> indexSet) {
    int n = (*A)[0].size();

    vector<double> x(n);

    //  printMatrix(A);
    //cout << "*x=" << endl;
    //printVector(b);
    // cout << "=============================" << endl;

    //bring into upper triangle form

    for (int i = 0; i < n - 1; i++) {
        int ithColumn = indexSet[i];
        // printMatrix(A);
        /*cout << "*x=" << endl;
        printVector(b);
        cout << "=============================" << endl;
        *///find pivot
        int maxIndex = i;
        for (int j = i; j < n; j++) {
            if ((*A)[ithColumn][j] > abs((*A)[ithColumn][maxIndex])) {
                maxIndex = j;
            }
        }
        //cout << "max index: " << maxIndex << endl;
        swapRows(A, maxIndex, i, indexSet);
        swapElements(b, maxIndex, i);

        /*printMatrix(A);
        cout << "*x=" << endl;
        printVector(b);
        cout << "=============================" << endl;
        */
        for (int j = i+1; j < n; j++) {
            double rowFactor = (*A)[ithColumn][j]/(*A)[ithColumn][i];
            for (int p = i+1; p < n; p++) {
                int k = indexSet[p];
                (*A)[k][j] -= (*A)[k][i]*rowFactor;
            }
            (*A)[ithColumn][j] = 0;

            (*b)[j] -= (*b)[i]*rowFactor;
            /*printMatrix(A);
            cout << "*x=" << endl;
            printVector(b);
            cout << "=============================" << endl;*/
        }
    }

}

//finds one possible solution x for Ax = b
//assumes A is quadratic
vector<double> solve(ColumnIndexedMatrix A, vector<double> b, vector<int> indexSet) {
   int n = indexSet.size();

    vector<double> x(n);
    transformIntoUpperTriangleMatrix(&A, &b, indexSet);
  /*  cout << "Upper triangle done!" << endl;
    printMatrix(A);
    cout << "*x=" << endl;
    printVector(b);
    cout << "=============================" << endl;*/
    //backsubstitution
    for (int i = n-1; i >= 0; i--) {
        x[i] = b[i];
        for(int j = i+1; j < n; j++) {
            int k = indexSet[j];
            x[i] -= A[k][i]*x[j];
        }

        x[i] /= A[indexSet[i]][i];
    }

    return x;
}

//returns A_B or x_B or c_B, etc. (is a copy not a reference to the old matrix/vector)
template <class T>
vector<T> indexMatrix(vector<T> x, vector<int> B) {
    vector<T> indexedVector(B.size());
    for (int i = 0; i < B.size(); i++) {
        indexedVector[i] = x[B[i]];
    }

    return indexedVector;
}

//returns the transpose of A_B as a new matrix. Does not change A.
ColumnIndexedMatrix transpose(vector<vector<double > > A, vector<int> indexSet) {
    ColumnIndexedMatrix ATranspose(indexSet.size());
    for (int i : indexSet) {
        vector<double> row(A.size());
        for (int j = 0; j < A.size(); j++) {
            row[j] = A[j][i];
        }
        ATranspose[i] = row;
    }

    return ATranspose;
}

//returns the transpose of A as a new matrix. Does not change A.
ColumnIndexedMatrix transpose(vector<vector<double > > A) {
    ColumnIndexedMatrix ATranspose(A[0].size());
    for (int i = 0; i < A[0].size(); i++) {
        vector<double> column(A.size());
        for (int j = 0; j < A.size(); j++) {
            column[j] = A[j][i];
        }
        ATranspose[i] = column;
    }

    return ATranspose;
}

//returns the dot product of x,y. aka x^Ty
double dotProduct(vector<double> x, vector<double> y) {
    double sum = 0;
    for (int i = 0; i < x.size(); i++) {
        sum += x[i]*y[i];

    }



    return sum;
}

//returns x*A_B as matrix vector product
vector<double> vectorMatrixMultiplication(ColumnIndexedMatrix A, vector<double> x, vector<int> indexSet) {
    vector<double> result(indexSet.size());
    for (int i = 0; i < indexSet.size(); i++) {
        for (int j = 0; j < A[0].size(); j++) {
            result[i] += A[indexSet[i]][j]*x[j];
        }
    }
    return result;
}

//returns x + factor*y
vector<double> add(vector<double> x, vector<double> y, double factor, vector<int> indexSet) {
    vector<double> result(indexSet.size());
    for (int i : indexSet) {
        result[i] = x[i] + factor*y[i];
    }

    return result;
}

//"calculates" N for the first time, for given B
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
        if (abs(x[i]) < pow(0.1,50)) {
            x[i] = 0;
        }
        if (x[i] > a) {
            return false;
        }
    }
    return true;
}

vector<int> vectorFromZeroToN(int n) {
    vector<int> fromOneToN(n);
    for (int i = 0; i < n; i++) {
        fromOneToN[i] = i;
    }

    return fromOneToN;
}

ColumnIndexedMatrix arrayToColumnIndexedMatrix(double **A, int m,int n){

    ColumnIndexedMatrix M(n);

    for(int i=0;i<n;i++){
        vector<double> column(m);

        for (int j=0;j<m;j++){
            column[j] = A[j][i];
        }
        M[i] = column;
    }
    return M;
}

/*
 * solves LP Ax = b, x >= 0, max c^Tx with given initial bases. Returns solution in optX, optB and result
 * Also returns the arbitrarily increasable x_j if LP is unbounded
 */
void simplex(ColumnIndexedMatrix A, vector<double> b, vector<double> c, vector<int> B, vector<double> * optX, vector<int> * optB, int * endJ, int * endResult) {
    //calc extra A transpose once.
    ColumnIndexedMatrix ATranspose = transpose(A);

    int n = A.size();
    int m = A[0].size();

    vector<int> allIndices = vectorFromZeroToN(n);



    int result = STILLSEARCHING;

    vector<int> N(n-m);
    calcN(B,&N,n);

    vector<double> x(n);
    vector<double> x_B = solve(A, b, B);

    for (int i = 0; i < B.size(); i++) {
        x[B[i]] = x_B[i];
    }
    int muh = 0;
    while(result == STILLSEARCHING) {
       /* cout << "B: ";
        printVector(B);
        cout << "N: ";
        printVector(N);*/
        //maxIndex for update
        int j,gamma,i;
        //cout << "current value: " <<  dotProduct(c, x) << "with x: ";
        //printVector(x);
        //BTRAN
        ColumnIndexedMatrix A_B = indexMatrix(A,B);
        ColumnIndexedMatrix AT_B = transpose(A_B);
        vector<double> y = solve(transpose(indexMatrix(A, B)), indexMatrix(c,B), vectorFromZeroToN(B.size()));
       /* for (int p = 0; p < y.size(); p++) {
            if (abs(y[p]) < pow(0.1, 50)) {
                y[p] = 0;
            }
        }*/
        //PRICING
        vector<double> A_Ny = vectorMatrixMultiplication(A, y, N);
        vector<double> c_N = indexMatrix(c,N);
        vector<double> z = add(c_N, A_Ny, -1, vectorFromZeroToN(N.size()));
        /*for (int p = 0; p < z.size(); p++) {
            if (abs(z[p]) < pow(0.1, 50)) {
                z[p] = 0;
            }
        }*/
        //printVector(z);
        if(lessThan(z, 0)) {
            //B is optimal
            result = OPTIMAL;
            //printVector(x);
            //cout << dotProduct(c, x) << endl;
            *endResult = OPTIMAL;
            *optX = x;
            *optB = B;
            break;
        } else {
            //Bland's rule
            int minIndex = INT_MAX;
            //double maxValue = -numeric_limits<double>::infinity();
            for (int k = 0; k < z.size(); k++) {
                if (z[k] > 0) {
                    if (N[k] < minIndex) {
                        minIndex = N[k];
                    }
                }
            }
            j = minIndex;
        }
        //FTRAN
        vector<double> jthColumnOfA = A[j];
        vector<double> w = solve(A, jthColumnOfA, B);
        //RATIOTEST (with Bland's rule)
        if(lessThan(w,0)) {
            result = UNBOUNDED;
            *endResult = UNBOUNDED;
            *optX = x;
            *optB = B;
            *endJ = j;
            break;
        } else {
            double min = INT_MAX;
            double minBi = INT_MAX;
            //vector<int> argmin;
            for(int k = 0; k < m; k++) {
                if (w[k] > 0 && x[B[k]]/w[k] <= min) {
                   /* if (x[B[k]]/w[k] < min) {
                        argmin.clear();
                    }*/
                 //   argmin.push_back(k);

                    if (B[k] < minBi) {
                        i = k;
                        minBi = B[k];
                    }
                    min = x[B[k]]/w[k];
                }
            }

            //findLexicographicSmallest(A, argmin, B, w)

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
                break;
            }
        }
        B[i] = j;
        /*for (int k : N) {
            x[k] = 0;
        }*/
        stable_sort(B.begin(), B.end());
        stable_sort(N.begin(), N.end());
    }

  //  cout << "Result: " << result << " after " << muh+1 << " iterations." << endl;


}

bool hasZeroRow(ColumnIndexedMatrix A) {
    bool foundZeroRow;
    for (int i = 0; i < A[0].size(); i++) {
        foundZeroRow = true;
        for (int j = 0; j < A.size(); j++) {
            if (A[j][i] != 0) {
                foundZeroRow = false;
            }
        }
        if (foundZeroRow) {
            return true;
        }
    }
    return false;
}

vector<int> extendToBasis(ColumnIndexedMatrix A, vector<int> B_x, int n) {
    vector<int> possibleBasisVectors;
    calcN(B_x, &possibleBasisVectors, n);

    for (int i : possibleBasisVectors) {
        vector<int> tempB_x = B_x;
        tempB_x.push_back(i);
        vector<double> unusedVector(A[0].size());
        ColumnIndexedMatrix tempA = indexMatrix(A, tempB_x);
        transformIntoUpperTriangleMatrix(&tempA, &unusedVector, vectorFromZeroToN(tempB_x.size()));
        if (!hasZeroRow(tempA)) {
            B_x.push_back(i);
        }
    }
}

void solveInequationalLPCompletely(const char * fileName) {
    int       m = 0;
    int       n = 0;
    double ** A1;
    double *  b1;
    double *  c1;

    if (read_LP(fileName, &m, &n, &A1, &b1, &c1) != EXIT_SUCCESS) {

        fprintf(stderr, "Error parsing \"%s\".\n", fileName);
        return;
    }

    ColumnIndexedMatrix A = array_to_vector_M(A1, m, n);
    vector<double> c(n);
    c.assign(c1, c1 + n);
    vector<double> b(m);
    b.assign(b1, b1 + m);

    get_standard(&A,&c,n);

  /*  cout << "A ===================" << endl;
    printMatrix(A);
    cout << "c ==================" << endl;
    printVector(c);*/

    vector<int> B(m);
    for (int i = 0; i < m; i++) {
        B[i] = A.size()-m+i;
    }

    //b >= 0?
    if (!lessThan(scalar_vec(b, -1), 0)) {
      /*  cout << "use auxiliary LP!" << endl;*/
        vector<double> aux_b(b.size());
        ColumnIndexedMatrix auxA = aux_lp(A, b, &aux_b);
        /*printMatrix(auxA);
        cout << "=================" << endl;*/
        vector <double> auxC = aux_lp(c, m);


//        B.clear();
        for (int i = 0; i < m; i++) {
            B[i] = auxA.size()-m+i;
        }

        vector<int> auxOptB(m);
        vector<double> auxOptX;
        int endJ;
        int auxResult;

        simplex(auxA, aux_b, auxC, B, &auxOptX, &auxOptB, &endJ, &auxResult);

        if (dotProduct(auxOptX, auxC) < 0) {
            //original LP is infeasible!
            cout << "Given LP is infeasible, since auxiliary LP has no solution with value 0!" << endl;
            return;
        }

        //checkifBisDegenerated
        bool degenerated = false;
        for (int i : auxOptB) {
            if (auxOptB[i] >= A.size()) {
                degenerated = true;
                break;
            }
        }

        if (degenerated) {
            bool done = false;
            while (!done) {

                //pivot step
                int i, j, B_i, gamma;
                gamma = 0;


                vector<int> auxOptN;

                calcN(auxOptB, &auxOptN, auxA.size());

                for (int j2 : auxOptN) {
                    //try basis change
                    if (j2 < n) {
                        vector<double> w = solve(auxA, A[j], auxOptB);
                        bool foundWi = false;
                        bool noArtificialVariables = true;
                        for (int i2 : auxOptB) {
                            if (i2 > n) {
                                noArtificialVariables = false;
                                if (w[i2] != 0) {
                                    B_i = auxOptB[i2];
                                    i = i2;
                                    foundWi = true;
                                }
                            }
                        }
                        if (noArtificialVariables) {
                            //found feasible basis!
                            done = true;
                            break;
                        } else if(!foundWi) {
                            //can extend to feasible basis!
                            auxOptB = extendToBasis(A, auxOptB, n);
                            done = true;
                            break;
                        } else {
                            auxOptX[j] = gamma;
                            //N := N \{j} union {B_i}
                            for (int k = 0; k < auxOptN.size(); k++) {
                                if (auxOptN[k] == j) {
                                    auxOptN[k] = B_i;
                                    break;
                                }
                            }
                            auxOptB[i] = j;
                            break;
                        }
                    }
                }

            }
        }

        B = auxOptB;


        //feasible basis B gleich so und so
    }
    //solve original LP with now known feasible basis B
    vector<int> optB(B.size());
    vector<double> optX(A.size());
    int result = 0;
    int endJ = 0; //for unboundness
    simplex(A,b,c,B, &optX, &optB, &endJ, &result);

    if (result == UNBOUNDED) {
        cout << "LP is unbounded! Since: " << endl;
        printVector(optX);
        cout << " can increased in " << endJ << " component arbitrary high, to reach any objective value without validating any constraints!"<< endl;
    } else  { //result = OPTIMAL
        //recalc solution of slacked Ax = b, x >= 0 version back to Ax <= 0, without nonnegativity constraints.
        vector<double> optSolutionForInitalLP(n);
        for (int i = 0; i < n; i++) {
            optSolutionForInitalLP[i] = optX[2*i] - optX[2*i+1];
        }
        cout << "Optimal feasible basic solution for LP is: " << endl;
        printVector(optSolutionForInitalLP);
        cout << "with value: " << dotProduct(c, optX) << "." <<endl;
    }
}

int main(int argc, const char * argv[]) {
    if (argc < 2) {
        fprintf(stderr, "Usage:  %s  <lp file>\n", argv[0]);
        return EXIT_FAILURE;
    }

    for (int i = 1; i <= argc-1; i++) {
        cout << "================================================================" << endl;
        cout << i << ". " << endl;
        solveInequationalLPCompletely(argv[i]);
    }



    return EXIT_SUCCESS;
}