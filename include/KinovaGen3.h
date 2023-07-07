#ifndef KINOVA_GEN_3_H_
#define KINOVA_GEN_3_H_

#include <Eigen/Dense>


template<class T>
class KinovaGen3
{
public:
    KinovaGen3();
    ~KinovaGen3();

    /**
     * Calculate the forward kinematics of the Kinova Gen3 robot
     * 
     * @param q Joint positions of the robot
     * @return The homogenous coordinate of the end-effector [R, x; 0, 1]
    */
    Eigen::Transform<T, 3, Eigen::Affine> forwardKinematics(const Eigen::Vector<T, 7>& q);

    /**
     * Calculate the inverse kinematics of the Kinova Gen3 robot at velocity level
     * 
     * @param q Joint positions of the robot
     * @param xp End-effector twist of the robot
     * @return Joint velocities of the robot
    */
    Eigen::Vector<T, 7> inverseKinematics(const Eigen::Vector<T, 7>& q,
                                          const Eigen::Vector<T, 6>& xp);

    /**
     * Calculate the Jacobian of the Kinova Gen3 robot
     * 
     * @param q Joint positions of the robot
     * @return The Jacobian of the robot evaluated at q
    */
    Eigen::Matrix<T, 6, 7> jacobian(const Eigen::Vector<T, 7>& q);

    /**
     * Calculate the mass matrix of the Kinova Gen3 robot
     * 
     * @param q Joint positions of the robot
     * @return The mass matrix of the robot evaluated at q
    */
    Eigen::Matrix<T, 7, 7> massMatrix(const Eigen::Vector<T, 7>& q);

    /**
     * Calculate the Coriolis matrix of the Kinova Gen3 robot
     * 
     * @param q Joint positions of the robot
     * @param qp Joint velocities of the robot
     * @return The Coriolis matrix of the robot evaluated at q and qp
    */
    Eigen::Vector<T, 7> coriolis(const Eigen::Vector<T, 7>& q,
                                      const Eigen::Vector<T, 7>& qp);

    /**
     * Calculate the gravity term of the Kinova Gen3 robot
     * 
     * @param q Joint positions of the robot
     * @return The gravity term of the robot evaluated at q
    */
    Eigen::Vector<T, 7> gravity(const Eigen::Vector<T, 7>& q);
};

#endif // KINOVA_GEN_3_H_