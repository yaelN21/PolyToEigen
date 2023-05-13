#include "TriangulationBase.h"

namespace Triangulation {

	class LinearLS : public TriangulationBase
	{
	public:
		using TriangulationBase::TriangulationBase;
		using TriangulationBase::triangulate;
		/**
		 *	\brief	Triangulate point using The Direct Linear Transformation Method.
		 *	\param	p0	Point in the image of the first camera.
		 *	\param	p1	Corresponding point in the image of the second camera.
		 *	\return	Triangulated point.
		 */
		Eigen::Vector3d triangulate(const Eigen::Vector2d& p0, const Eigen::Vector2d& p1) const override;
	};

}