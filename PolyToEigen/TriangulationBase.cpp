#include "TriangulationBase.h"

namespace Triangulation {

	Triangulation::TriangulationBase::TriangulationBase(const Eigen::MatrixXd& P0, Eigen::MatrixXd& P1)
		: P0(P0), P1(P1) {}

	std::vector<Eigen::Vector3d> TriangulationBase::triangulate(const std::vector<Eigen::Vector2d>& p0, const std::vector<Eigen::Vector2d>& p1) const
	{
		std::vector<Eigen::Vector3d> result;
		for (size_t i = 0; i < p0.size(); ++i)
		{
			result.emplace_back(triangulate(p0[i], p1[i]));
		}
		return result;
	}

}