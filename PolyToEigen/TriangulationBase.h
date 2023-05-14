#include <iostream>
#include <Eigen/Dense>
#include <vector>


using Eigen::Matrix3d;
using Eigen::Vector3d;

using Eigen::MatrixXd;

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
		 */
		

		TriangulationBase(const Eigen::MatrixXd& P0, const Eigen::MatrixXd& P1);
		virtual ~TriangulationBase() {}
			/**
			 *	\brief	Triangulates image points.
			 *	\param	p0	Point in the image of the first camera.
			 *	\param	p1	Corresponding point in the image of the second camera.
			 *	\return	Triangulated point.
			 */
		virtual Eigen::Vector3d triangulate(const Eigen::Vector2d& p0, const Eigen::Vector2d& p1) const = 0;
			/**
			 *	\brief	Triangulate points using The Direct Linear Transformation Method.
			 *	\param	p0	Points in the image of the first camera.
			 *	\param	p1	Corresponding points in the image of the second camera.
			 *	\return	Triangulated points.
			 */
		std::vector<Eigen::Vector3d> triangulate(const std::vector<Eigen::Vector2d>& p0, const std::vector<Eigen::Vector2d>& p1) const;

	protected:
		const Eigen::MatrixXd P0, P1;
	};
#endif  // TRIANGULATIONBASE_H
}	// Triangulation

