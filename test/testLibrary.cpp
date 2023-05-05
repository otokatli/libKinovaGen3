#include "KinovaGen3.h"
#include <Eigen/Dense>
#include <iostream>


int main()
{
    KinovaGen3 robot;

    Eigen::Vector<double, 7> q0(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
    Eigen::Vector<double, 7> qp0(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0);

    Eigen::Vector<double, 3> x;
    Eigen::Matrix<double, 3, 3> R;
    Eigen::Matrix<double, 6, 7> J;
    Eigen::Matrix<double, 7, 7> M;
    Eigen::Vector<double, 7> C;
    Eigen::Vector<double, 7> G;

    robot.forwardKinematics(q0, x, R);
    std::cout << "End effector position:" << std::endl;
    std::cout << x << std::endl;
    std::cout << "End effector orientation:" << std::endl;
    std::cout << R << std::endl;

    robot.inverseKinematics();

    robot.jacobian(q0, J);
    std::cout << "Jacobian:" << std::endl;
    std::cout << J << std::endl;

    robot.massMatrix(q0, M);
    std::cout << "Mass matrix:" << std::endl;
    std::cout << M << std::endl;
    
    robot.coriolis(q0, qp0, C);
    std::cout << "Coriolis:" << std::endl;
    std::cout << C << std::endl;

    robot.gravity(q0, G);
    std::cout << "Gravity:" << std::endl;
    std::cout << G << std::endl;

}