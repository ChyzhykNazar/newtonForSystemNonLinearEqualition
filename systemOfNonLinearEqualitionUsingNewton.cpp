/*Знайти розв'язок системи нелінійних алгебраїчних рівнянь методом Ньютона*/
/*To find solution for system of nonlinear equation using Newton method*/
#include <iostream>

using namespace std;

void methodGaus(double** a, double* b, double* x, int n) {
    double maxValue, s, d;
    int maxIndex;
    int p = 0;
    while (p < n)
    {
        // Search for the string with the maximum value
        maxValue = abs(a[p][p]);
        maxIndex = p;
        for (int i = p + 1; i < n; i++)
        {
            if (abs(a[i][p]) > maxValue)
            {
                maxValue = abs(a[i][p]);        //remember the greatest value
                maxIndex = i;                   //remember the index
            }
        }
        //permutation of lines
        if (maxValue < 0.00001)
        {
            //there are no non-zero diagonal elements
            cout << "The solution cannot be obtained because of the zero column ";
            cout << maxIndex << " matrix A" << endl;
            return;
        }
        for (int j = 0; j < n; j++)
        {
            swap(a[maxIndex][j], a[p][j]);      //rearrange the rows in the matrix
        }
        swap(b[maxIndex], b[p]);                //rearrange the coefficient
        p++;
    }
    // Applying Gauss Elimination
    for (int i = 0; i < n; i++) {
        for (int j = i + 1; j < n; j++) {
            s = a[j][i] / a[i][i];    
            for (int q = i; q < n; q++) {
                a[j][q] = a[j][q] - s * a[i][q];        //subtraction in the array
            }
            b[j] = b[j] - s * b[i];                     //subtraction in the result
            //printFullMatrix(a, b, n);
            //cout << endl;
        }
    }
    // Obtaining Solution by Back Substitution Method
    for (int i = n - 1; i >= 0; i--) {
        d = 0;
        s = 0;
        for (int j = i + 1; j < n; j++) {
            s = a[i][j] * x[j];     //find each multiplier  
            d = d + s;              //sum all the members of the matrix
        }
        x[i] = (b[i] - d) / a[i][i];    //finding roots
    }

    for (int i = 0; i < n; i++) {
        cout << "x[" << i << "]=" << x[i] << " " << endl << endl;       //displaying solution
    }
}

//The derivative of the first function in x
double df1_dx(double* x)
{
    return 8 + (4 * cos(8 * x[0] + x[1] - 36)) / 5;
}

//The derivative of the first function in y
double df1_dy(double* x)
{
    return cos(8 * x[0] + x[1] - 36) / 10;
}

//The derivative of the second function in x
double df2_dx(double* x)
{
    return -(2 * cos(8 * x[0] + x[1] - 36)) / 25;
}

//the derivative of the second function in y
double df2_dy(double* x)
{
    return -(cos(8 * x[0] + x[1] - 36)) / 100 + 1;
}

//the first function
double f1(double* x)
{
    return 8 * x[0] - 36 + (sin(8 * x[0] + x[1] - 36)) / 10;
}

//the second funtion
double f2(double* x)
{
    return x[1] - (sin(8 * x[0] + x[1] - 36)) / 100;
}

void newton(double eps) {
    cout << "___________________________Newton's method________________________________" << endl;
    int n = 2;
    int iter = 0;
    double** a = new double* [n];
    for (int i = 0; i < n; i++) {
        a[i] = new double[n];
    }
    double* b = new double[n];
    double* x = new double[n];
    double* xk = new double[n];
    double* p = new double[n];
    double norm = 0;
    // Choose the initial approximation
    xk[0] = 1; xk[1] = 1;
    do
    {
        iter++;
        for (int i = 0; i < n; i++)
        {
            x[i] = xk[i];
        }
        // Fill the matrix A with derivatives of functions
        a[0][0] = df1_dx(x);
        a[0][1] = df1_dy(x);
        a[1][0] = df2_dx(x);
        a[1][1] = df2_dy(x);
        // fill the matrix B
        b[0] = -f1(x); b[1] = -f2(x);
        // use the Gaussian method to solve the system of linear equation
        methodGaus(a, b, p, n);
        //  get an approximation by adding the previous and received
        for (int i = 0; i < n; i++)
        {
            xk[i] = x[i] + p[i];
        }
        norm = abs(xk[0] - x[0]);
        for (int i = 0; i < n; i++) {
            if (abs(xk[i] - x[i]) > norm) {
                norm = abs(xk[i] - x[i]);
            }
            x[i] = xk[i];
        }
    } while (norm >= eps);

    for (int i = 0; i < n; i++)
    {
        cout << "x[" << i << "] = " << xk[i] << endl << endl;
    }
    cout << "Iter: " << iter << endl << endl;
};

int main()
{
    double eps = 0.01;
    newton(eps);
}
