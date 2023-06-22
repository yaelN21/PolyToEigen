#include "TriangulationBase.h"

namespace Triangulation {

	Triangulation::TriangulationBase::TriangulationBase(const Eigen::MatrixXf& P0,const Eigen::MatrixXf& P1,
		const Eigen::MatrixXf& K0,const Eigen::MatrixXf& K1,const Eigen::MatrixXf& R0,const Eigen::MatrixXf& R1,
		const  Eigen::Vector3f& T0,const  Eigen::Vector3f& T1):
		P0(P0),P1(P1),K0(K0), K1(K1), R0(R0), R1(R1), T0(T0), T1(T1) {}

		Triangulation::TriangulationBase::TriangulationBase(const Eigen::MatrixXf& P0,const Eigen::MatrixXf& P1):
		P0(P0),P1(P1){}	


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