#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <cassert>
#include <cmath>
#include <cholesky.h>
#include <triangular_solve.h>
#include <Eigen/Dense>

int main(int argc, char* argv[])
{
    // Output precision
    int p=9;
    // Matrix dimension
    int N=9;
    // Matrix specification
    Eigen::MatrixXd A = Eigen::MatrixXd::Zero(N,N); 
    Eigen::MatrixXd A1 = Eigen::MatrixXd::Zero(N,N); 
    for (int i=0; i<N; i++) {
        A1(i,i) = 4.0;
        if (i>=1) A1(i,i-1) = -1.0;
        if (i>=2) A1(i,i-2) = 3.0;
        if (i<=6) A1(i,i+2) = -2.0;
    }
    // Could be AT+A or AT*A
    A = A1.transpose()+A1;
    // Vector specification
    Eigen::VectorXd b = Eigen::VectorXd::Zero(N);
    for (int i=0; i<N; i++) {
        b(i) = (i*i-9.0)/(i+5.0); 
    }

    // Cholesky decomposition
    Eigen::MatrixXd U = cholesky(A);
    Eigen::MatrixXd L = U.transpose();

    // Perform forward substitution and backward substitution
    // Using Eigen's library method: Eigen::VectorXd x = A.llt().solve(b);
    Eigen::VectorXd y = forward_subst(L, b);
    Eigen::VectorXd x = backward_subst(U, y);

    // Compute residual error
    double R = (A*x-b).norm();
    std::cout << "Residual error =, "
              << std::setprecision(p)
              << R
              << std::endl;

    // Compute L-inverse, U-inverse, and A-inverse
    Eigen::MatrixXd L_inv = Eigen::MatrixXd::Zero(N,N);
    Eigen::MatrixXd U_inv = Eigen::MatrixXd::Zero(N,N);
    Eigen::MatrixXd A_inv = Eigen::MatrixXd::Zero(N,N);
    // L-inverse
    for (int i=0; i<N; i++) {
        Eigen::VectorXd ei = Eigen::VectorXd::Zero(N);
        ei(i) = 1;
        Eigen::VectorXd x = forward_subst(L, ei);
        L_inv.col(i) = x;
    }
    // U-inverse
    for (int i=0; i<N; i++) {
        Eigen::VectorXd ei = Eigen::VectorXd::Zero(N);
        ei(i) = 1;
        Eigen::VectorXd x = backward_subst(U, ei);
        U_inv.col(i) = x;
    }
    // A-inverse
    for (int i=0; i<N; i++) {
        Eigen::VectorXd ei = Eigen::VectorXd::Zero(N);
        ei(i) = 1;
        Eigen::VectorXd y = forward_subst(L, ei);
        Eigen::VectorXd x = backward_subst(U, y);
        A_inv.col(i) = x;
    }
    
    // Compute residual error by matrix inverse
    Eigen::VectorXd v = A_inv*b;
    R = (b-A*v).norm();
    std::cout << "Residual error by inverse =, "
              << std::setprecision(p)
              << R
              << std::endl;

    std::cout << "Cholesky decomposition" << std::endl;
    for (int i=0; i<N; i++) {
        std::cout << " |, ";
        for (int j=0; j<N; j++) {
            std::cout << std::right 
                      << std::setw(p+3) 
                      << std::fixed 
                      << std::setprecision(p) 
                      << A(i,j) 
                      << ", ";
        }
        std::cout << " |, =, |, ";
        for (int j=0; j<N; j++) {
            std::cout << std::right 
                      << std::setw(p+3) 
                      << std::fixed 
                      << std::setprecision(p) 
                      << L(i,j) 
                      << ", ";
        }
        std::cout << " |, ";
        for (int j=0; j<N; j++) {
            std::cout << std::right 
                      << std::setw(p+3) 
                      << std::fixed 
                      << std::setprecision(p) 
                      << U(i,j) 
                      << ", ";
        }
        std::cout << " |, " << std::endl;
    }

    std::cout << "Linear solve - x:" << std::endl;
    for (int i=0; i<N; i++) {
        std::cout << std::right 
                  << std::setw(p+3) 
                  << std::fixed 
                  << std::setprecision(p) 
                  << x(i)
                  << std::endl; 
    }
    
    std::cout << "Inverse relations:" << std::endl;
    for (int i=0; i<N; i++) {
        std::cout << " |, ";
        for (int j=0; j<N; j++) {
            std::cout << std::right 
                      << std::setw(p+3) 
                      << std::fixed 
                      << std::setprecision(p) 
                      << A_inv(i,j) 
                      << ", ";
        }
        std::cout << " |, =, |, ";
        for (int j=0; j<N; j++) {
            std::cout << std::right 
                      << std::setw(p+3) 
                      << std::fixed 
                      << std::setprecision(p) 
                      << U_inv(i,j) 
                      << ", ";
        }
        std::cout << " |, ";
        for (int j=0; j<N; j++) {
            std::cout << std::right 
                      << std::setw(p+3) 
                      << std::fixed 
                      << std::setprecision(p) 
                      << L_inv(i,j) 
                      << ", ";
        }
        std::cout << " |, " << std::endl;
    }
    
    std::cout << "Linear solve by inversion- x:" << std::endl;
    for (int i=0; i<N; i++) {
        std::cout << std::right 
                  << std::setw(p+3) 
                  << std::fixed 
                  << std::setprecision(p) 
                  << v(i)
                  << std::endl; 
    }
    
    // Verify accuracy of matrix inverse
    std::cout << "L * L-inverse residual error =, "
              << std::scientific
              << std::setprecision(p)
              << (L*L_inv-Eigen::MatrixXd::Identity(N,N)).norm() 
              << std::endl;
    std::cout << "U * U-inverse residual error =, "
              << std::scientific
              << std::setprecision(p)
              << (U*U_inv-Eigen::MatrixXd::Identity(N,N)).norm() 
              << std::endl;
    std::cout << "A * A-inverse residual error =, "
              << std::scientific
              << std::setprecision(p)
              << (A*A_inv-Eigen::MatrixXd::Identity(N,N)).norm() 
              << std::endl;

    return 0;
}

