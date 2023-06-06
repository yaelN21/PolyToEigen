#include "TriangulationBase.h"

namespace Triangulation {

	Triangulation::TriangulationBase::TriangulationBase(const Eigen::MatrixXf& P0,const Eigen::MatrixXf& P1):
	 P0(P0), P1(P1) {
		//test
	}

	std::vector<Eigen::Vector3f> TriangulationBase::triangulate(const std::vector<Eigen::Vector2f>& p0, const std::vector<Eigen::Vector2f>& p1) const
	{
		std::vector<Eigen::Vector3f> result;
		for (size_t i = 0; i < p0.size(); ++i)
		{
			result.emplace_back(triangulate(p0[i], p1[i]));
		}
		return result;
	}

}