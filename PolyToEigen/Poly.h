#include <iostream>
#include <Eigen/Dense>
#include "TriangulationBase.h"

#include "PolyBase.h"

namespace Triangulation {

	/**
	 *	\brief	Performs non-linear triangulation of two image points by minimizing geometric error.
	 */
	class Poly : public PolyBase {
	public:
		using PolyBase::PolyBase;
		using PolyBase::TriangulationBase::triangulate;
	private:
		std::vector<float> PreparePolyCoeffs(const PolyParams& params) const override;
		std::vector<float> EvaluateRootsCosts(const Roots& roots, const PolyParams& params) const override;
	};
}