#include <iostream>
#include <iomanip>
#include "Eigen/Eigen"


using namespace std;
using namespace Eigen;



// return true if the matrix A is invertible, else it returns false
//
bool CheckSystem(const Matrix2d& A)
{
    if(abs(A.determinant()) < 1e-16)
    {
        return false; //determinant is smaller than eps, matrix is singular
    }
    return true;
}

//compute the systems and return a vector of the relatives errors [errLU, errQR]
//
Vector2d SolveSystem(const Matrix2d& A, const Vector2d& b)
{
    Vector2d sol = Vector2d::Zero();
    Vector2d x_solution(-1.0e+0, -1.0e+0);

    Vector2d x_LU = A.fullPivLu().solve(b); //solution with PA = LU factorization
    Vector2d x_QR = A.householderQr().solve(b); //solution with A = QR factorization

    //cout<< scientific << setprecision(16)<< x_solution << endl;
    //cout <<scientific<<setprecision(16)<< x_LU<<endl;
    //cout <<scientific<<setprecision(16)<< x_QR <<endl;

    double err_LU = (x_LU - x_solution).norm() / x_solution.norm();
    double err_QR = (x_QR - x_solution).norm() / x_solution.norm();

    sol(0) = err_LU;
    sol(1) = err_QR;

    return sol;
}




Matrix2d X1 {
    {5.547001962252291e-01, -3.770900990025203e-02},
    { 8.320502943378437e-01, -9.992887623566787e-01}
};
Vector2d b1(-5.169911863249772e-01, 1.672384680188350e-01);

Matrix2d X2 {
    {5.547001962252291e-01, -5.540607316466765e-01},
    {8.320502943378437e-01, -8.324762492991313e-01}
};
Vector2d b2(-6.394645785530173e-04, 4.259549612877223e-04);

Matrix2d X3 {
    {5.547001962252291e-01, -5.547001955851905e-01},
    {8.320502943378437e-01, -8.320502947645361e-01}
};
Vector2d b3(-6.400391328043042e-10, 4.266924591433963e-10);


int main()
{
    Vector2d results; // vettore con i risulati, contiene [errRel_LP, errRel_QR]
    if (CheckSystem(X1)){
        results = SolveSystem(X1,b1);
        cout << scientific << setprecision(16)<< "Relative Error LP: "<< results(0) << ", Relative Error QR: "<< results(1) << endl;}
    else{
        cout << "Matrix is singular"<< endl;}


    Vector2d results2;
    if (CheckSystem(X2)){
        results2 = SolveSystem(X2,b2);
        cout << scientific<< setprecision(16)<< "Relative Error LP: "<< results2(0) << ", Relative Error QR: "<< results2(1) << endl;}
    else{
        cout << "Matrix is singular"<< endl;}


    Vector2d results3;
    if (CheckSystem(X3)){
        results3 = SolveSystem(X3,b3);
        cout << scientific<< setprecision(16)<< "Relative Error LP: "<< results3(0) << ", Relative Error QR: "<< results3(1) << endl;}
    else{
        cout << "Matrix is singular"<< endl;}

}

