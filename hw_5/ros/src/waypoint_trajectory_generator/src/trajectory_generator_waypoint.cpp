#include "trajectory_generator_waypoint.h"
#include <stdio.h>
#include <ros/ros.h>
#include <ros/console.h>
#include <iostream>
#include <fstream>
#include <string>
#include <math.h>

using namespace std;    
using namespace Eigen;

TrajectoryGeneratorWaypoint::TrajectoryGeneratorWaypoint(){}
TrajectoryGeneratorWaypoint::~TrajectoryGeneratorWaypoint(){}

//define factorial function, input i, output i!
int TrajectoryGeneratorWaypoint::Factorial(int x)
{
    return x*(x-1)*(x-2)*(x-3);
}

MatrixXd TrajectoryGeneratorWaypoint::blkdiag(const MatrixXd& a, const MatrixXd b)
{
    MatrixXd bdm = MatrixXd::Zero(a.rows() + 8, a.cols() + 8);
    bdm.block(0, 0, a.rows(), a.cols()) = a;
    bdm.block(a.rows(), a.cols() , b.rows(), b.cols()) = b;
    return bdm;
}

MatrixXd TrajectoryGeneratorWaypoint::getCoeff(double t)
{
    MatrixXd coeff = MatrixXd::Zero(4,8);
    coeff << 1, 1*t, 1*pow(t,2), 1*pow(t,3), 1*pow(t,4),  1*pow(t,5),  1*pow(t,6),   1*pow(t,7),
             0,   1, 2*pow(t,1), 3*pow(t,2), 4*pow(t,3),  5*pow(t,4),  6*pow(t,5),   7*pow(t,6),
             0,   0,          2,    6*pow(t,1), 12*pow(t,2), 20*pow(t,3), 30*pow(t,4),  42*pow(t,5),
             0,   0,          0,                6, 24*pow(t,1), 60*pow(t,2), 120*pow(t,3), 210*pow(t,4);
    return coeff;
}

/*
    STEP 2: Learn the "Closed-form solution to minimum snap" in L5, then finish this PolyQPGeneration function
    variable declaration: input       const int d_order,                    // the order of derivative
                                      const Eigen::MatrixXd &Path,          // waypoints coordinates (3d)
                                      const Eigen::MatrixXd &Vel,           // boundary velocity
                                      const Eigen::MatrixXd &Acc,           // boundary acceleration
                                      const Eigen::VectorXd &Time)          // time allocation in each segment
                          output      MatrixXd PolyCoeff(m, 3 * p_num1d);   // position(x,y,z), so we need (3 * p_num1d) coefficients
*/

Eigen::MatrixXd TrajectoryGeneratorWaypoint::PolyQPGeneration(
            const int d_order,                    // the order of derivative
            const Eigen::MatrixXd &Path,          // waypoints coordinates (3d)
            const Eigen::MatrixXd &Vel,           // boundary velocity
            const Eigen::MatrixXd &Acc,           // boundary acceleration
            const Eigen::MatrixXd &Jerk,           // boundary acceleration
            const Eigen::VectorXd &Time)          // time allocation in each segment
{
    // enforce initial and final velocity and accleration, for higher order derivatives, just assume them be 0;
    int p_order   = 2 * d_order - 1;              // the order of polynomial
    int p_num1d   = p_order + 1;                  // the number of variables in each segment

    int m = Time.size();                          // the number of segments
    MatrixXd PolyCoeff = MatrixXd::Zero(m, 3 * p_num1d);           // position(x,y,z), so we need (3 * p_num1d) coefficients
    VectorXd Px(p_num1d * m), Py(p_num1d * m), Pz(p_num1d * m);

    /*   Produce Mapping Matrix A to the entire trajectory, A is a mapping matrix that maps polynomial coefficients to derivatives.   */
    // Compute Q
    MatrixXd Q = MatrixXd::Zero(p_num1d, p_num1d);;
    for(int k=1; k<=m; k++)
    {
        MatrixXd Q_k = MatrixXd::Zero(p_num1d, p_num1d);
        for(int i = 0; i < p_order; i++) {
            for (int j = 0; j < p_order; j++) {
                if (i < 4 || j < 4)
                    continue;
                else {
                    Q_k(i, j) = Factorial(i) * Factorial(j) / (i + j  - 7) * pow(Time(k-1), i + j - 7);
                }
            }
        }
        if(k == 1)
            Q = Q_k;
        else
        {
            Q = blkdiag(Q,Q_k);
        }
    }

    cout << "---------------- Q Finish-------------" << endl;
    cout <<"Q size(row,col):" << Q.rows() <<  "," << Q.cols() << endl;

    // Compute M
    MatrixXd coeff1 = MatrixXd::Zero(p_num1d/2, p_num1d);
    MatrixXd coeff2 = MatrixXd::Zero(p_num1d/2, p_num1d);
    MatrixXd M = MatrixXd::Zero(p_num1d, p_num1d);
   for(int k=1; k<=m; k++)
   {
       MatrixXd M_k = MatrixXd::Zero(p_num1d, p_num1d);
       coeff1 = getCoeff(0);
       coeff2 = getCoeff(Time(k-1));
       M_k << coeff1,
              coeff2;
       if(k == 1)
           M = M_k;
       else
           M =  blkdiag(M, M_k);
   }

    cout << "------------------ M Finish --------------------------" << endl;
    cout <<"M size(row,col):" << M.rows() <<  "," << M.cols() << endl;

    // Compute Ct
    MatrixXd Ct = MatrixXd::Zero(p_num1d*m, p_num1d/2*m + p_num1d/2);
    Eigen::MatrixXd eye3 = MatrixXd::Identity(3, 3);
    Eigen::MatrixXd eye4 = MatrixXd::Identity(4, 4);
    int count = 0;
    for(int i = 1; i<=p_num1d*m/(p_num1d/2); i++)
    {
        if(i == 1)
        {
            Ct.block<4,4>(0,0) = eye4;
            count = count + 4;
            continue;
        }

        if (i == p_num1d*m/(p_num1d/2))
        {
            count++;
            Ct.block<4,4>(p_num1d*m-3-1,count-1) = eye4;
            continue;
        }

        if(i%2 == 0)
            count++;

        Ct((i-1)*4+1-1,count-1) = 1;
        Ct.block<3,3>((i-1)*4+2-1, 4*m+2-3*(m-2)+(int(i/2)-1)*3 -1) = eye3;
    }

    cout << "------------------ Ct Finish --------------------------" << endl;
    cout <<"Ct size(row,col):" << Ct.rows() <<  "," << Ct.cols() << endl;

    MatrixXd C = Ct.transpose();
    MatrixXd M_inv = M.inverse();
    MatrixXd R = C * M_inv.transpose() * Q * M_inv * Ct;

    MatrixXd R_pp = R.block(m+7, m+7, 3*(m-1), 3*(m-1));
    MatrixXd R_fp = R.block(0, m+7, m+7, 3*(m-1));

    /*   Produce the dereivatives in X, Y and Z axis directly.  */
    MatrixXd dFx = MatrixXd::Zero(p_num1d + m-1, 1);
    MatrixXd dFy = MatrixXd::Zero(p_num1d + m-1, 1);
    MatrixXd dFz = MatrixXd::Zero(p_num1d + m-1, 1);

    dFx(0) = Path(0,0);
    dFx(1) = 0;
    dFx(2) = 0;
    dFx(3) = 0;

    dFy(0) = Path(0,1);
    dFy(1) = 0;
    dFy(2) = 0;
    dFy(3) = 0;

    dFz(0) = Path(0,2);
    dFz(1) = 0;
    dFz(2) = 0;
    dFz(3) = 0;

    for(int i = 1; i <= m; i++)
    {
        dFx(3+i) = Path(i,0);
        dFy(3+i) = Path(i,1);
        dFz(3+i) = Path(i,2);
    }

    dFx(p_num1d + m-1-1) = 0;
    dFx(p_num1d + m-1-2) = 0;
    dFx(p_num1d + m-1-3) = 0;

    dFy(p_num1d + m-1-1) = 0;
    dFy(p_num1d + m-1-2) = 0;
    dFy(p_num1d + m-1-3) = 0;

    dFz(p_num1d + m-1-1) = 0;
    dFz(p_num1d + m-1-2) = 0;
    dFz(p_num1d + m-1-3) = 0;

    MatrixXd dpx = -R_pp.inverse() * R_fp.transpose() * dFx;
    MatrixXd dpy = -R_pp.inverse() * R_fp.transpose() * dFy;
    MatrixXd dpz = -R_pp.inverse() * R_fp.transpose() * dFz;
    /*   Produce the Minimum Snap cost function, the Hessian Matrix   */

    MatrixXd dFdPx = MatrixXd::Zero(p_num1d/2*m + 4, 1);
    MatrixXd dFdPy = MatrixXd::Zero(p_num1d/2*m + 4, 1);
    MatrixXd dFdPz = MatrixXd::Zero(p_num1d/2*m + 4, 1);
    dFdPx << dFx,
             dpx;
    dFdPy << dFy,
             dpy;
    dFdPz << dFz,
             dpz;

    MatrixXd PolyCoeffx = MatrixXd::Zero(p_num1d*m, 1);
    MatrixXd PolyCoeffy = MatrixXd::Zero(p_num1d*m, 1);
    MatrixXd PolyCoeffz = MatrixXd::Zero(p_num1d*m, 1);
    PolyCoeffx = M_inv * Ct * dFdPx;
    PolyCoeffy = M_inv * Ct * dFdPy;
    PolyCoeffz = M_inv * Ct * dFdPz;

    Map<MatrixXd> PolyCoeffx2(PolyCoeffx.data(),p_num1d,m);
    Map<MatrixXd> PolyCoeffy2(PolyCoeffy.data(),p_num1d,m);
    Map<MatrixXd> PolyCoeffz2(PolyCoeffz.data(),p_num1d,m);

    PolyCoeff << PolyCoeffx2.transpose(),PolyCoeffy2.transpose(),PolyCoeffz2.transpose();

    cout << "------------------ PolyCoeff Finish --------------------------" << endl;
    cout << "PolyCoeff sizes(row,col):" << PolyCoeff.rows() << "," << PolyCoeff.cols() <<endl;
    return PolyCoeff;
}
