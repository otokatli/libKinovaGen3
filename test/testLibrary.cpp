#include "KinovaGen3.h"
#include <Eigen/Dense>
#include <iostream>

int main()
{
    KinovaGen3 robot;

    Eigen::Vector<double, 7> q0(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
    Eigen::Vector<double, 7> qp0(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0);

    Eigen::Transform<double, 3, Eigen::Affine> T;
    Eigen::Matrix<double, 6, 7> J;
    Eigen::Matrix<double, 7, 7> M;
    Eigen::Vector<double, 7> C;
    Eigen::Vector<double, 7> G;

    T = robot.forwardKinematics(q0);
    std::cout << "End effector position:" << std::endl;
    std::cout << T.translation() << std::endl;
    std::cout << "End effector orientation:" << std::endl;
    std::cout << T.rotation() << std::endl;

    // robot.inverseKinematics();

    J = robot.jacobian(q0);
    std::cout << "Jacobian:" << std::endl;
    std::cout << J << std::endl;

    M = robot.massMatrix(q0);
    std::cout << "Mass matrix:" << std::endl;
    std::cout << M << std::endl;
    
    C = robot.coriolis(q0, qp0);
    std::cout << "Coriolis:" << std::endl;
    std::cout << C << std::endl;

    G = robot.gravity(q0);
    std::cout << "Gravity:" << std::endl;
    std::cout << G << std::endl;

}