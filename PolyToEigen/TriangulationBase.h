#include <iostream>
#include <Eigen/Dense>
#include <vector>


using Eigen::Matrix3f;
using Eigen::Vector3f;

using Eigen::MatrixXf;

namespace Triangulation {
		/**
	 * 	\class	TriangulationBase
	 *	\brief
	 */

#ifndef TRIANGULATIONBASE_H
#define TRIANGULATIONBASE_H

class TriangulationBase {
	public:
		/**
		 * 
		 *	\brief	Constructor
		 *	\param	P0	Camera matrix of the first camera.
		 *	\param 	P1	Camera matrix of the second camera.
		 *	\param 	K
		 */
		

		TriangulationBase(const Eigen::MatrixXf& P0, const Eigen::MatrixXf& P1);
		virtual ~TriangulationBase() {}
			/**
			 *	\brief	Triangulates image points.
			 *	\param	p0	Point in the image of the first camera.
			 *	\param	p1	Corresponding point in the image of the second camera.
			 *  \param 	K
			 *	\return	Triangulated point.
			 */
		virtual Eigen::Vector3f triangulate(const Eigen::Vector2f& p0, const Eigen::Vector2f& p1) const = 0;
			/**
			 *	\brief	Triangulate points using The Direct Linear Transformation Method.
			 *	\param	p0	Points in the image of the first camera.
			 *	\param	p1	Corresponding points in the image of the second camera.
			 *	\return	Triangulated points.
			 */
		std::vector<Eigen::Vector3f> triangulate(const std::vector<Eigen::Vector2f>& p0, const std::vector<Eigen::Vector2f>& p1) const;

	protected:
		const Eigen::MatrixXf P0, P1;
	};
#endif  // TRIANGULATIONBASE_H
}	// Triangulation

