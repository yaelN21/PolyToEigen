#include "../include/LinearLS.h"


namespace Triangulation {
    Eigen::Vector3f Triangulation::LinearLS::triangulate(const Eigen::Vector2f& p0, const Eigen::Vector2f& p1) const
    {
        Eigen::Matrix<float, 4, 4> A = Eigen::Matrix<float, 4, 4>::Zero();
        A.row(0) = p0.x() * P0.row(2) - P0.row(0);
        A.row(1) = p0.y() * P0.row(2) - P0.row(1);
        A.row(2) = p1.x() * P1.row(2) - P1.row(0);
        A.row(3) = p1.y() * P1.row(2) - P1.row(1);

        Eigen::JacobiSVD<Eigen::Matrix<float, 4, 4>> svd(A, Eigen::ComputeFullU | Eigen::ComputeFullV);
        Eigen::Vector4f result = svd.matrixV().col(3);

        return Eigen::Vector3f(result[0] / result[3], result[1] / result[3], result[2] / result[3]);
	}

}