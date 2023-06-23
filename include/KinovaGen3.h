#include <Eigen/Dense>
#include <tuple>

class KinovaGen3
{
public:
    KinovaGen3();
    ~KinovaGen3();

    Eigen::Transform<double, 3, Eigen::Affine> forwardKinematics(const Eigen::Vector<double, 7> q);

    Eigen::Vector<double, 7> inverseKinematics(const Eigen::Vector<double, 7> q,
                                               const Eigen::Vector<double, 6> xp);

    Eigen::Matrix<double, 6, 7> jacobian(const Eigen::Vector<double, 7> q);

    Eigen::Matrix<double, 7, 7> massMatrix(const Eigen::Vector<double, 7> q);

    Eigen::Vector<double, 7> coriolis(const Eigen::Vector<double, 7> q,
                                      const Eigen::Vector<double, 7> qp);

    Eigen::Vector<double, 7> gravity(const Eigen::Vector<double, 7> q);
};