#pragma once
#include <iostream>
#include <Eigen/Dense>
#include <complex>
#include "TriangulationBase.h"
#include "LinearLS.h"
#include <vector>


using Eigen::Matrix3f;
using Eigen::Vector3f;

namespace Triangulation {
	class PolyBase : public TriangulationBase
	{
		typedef Eigen::Matrix3f Intrinsic;
		typedef Eigen::Matrix3f Fundamental;

		/**
		 *	\brief	Constructor
		 *	\param	F	Fundamental matrix.
		 */
		
		

		public:
		PolyBase(const Eigen::MatrixXf& P0, const Eigen::MatrixXf& P1, const Eigen::MatrixXf& K);
			/**
			 *	\brief	Constructor
			 *	\param	P0	Camera matrix of the first camera.
			 *	\param 	P1	Camera matrix of the second camera.
			  *	\param K
			 *	\param	F	Fundamental matrix.
			 */

		explicit PolyBase(const Fundamental& F, const Eigen::MatrixXf& K);
		/**
		 *	\brief	Constructor
		 *	\param	P0	Camera matrix of the first camera.
		 *	\param 	P1	Camera matrix of the second camera.
		 *\param 	K
		 */
		
		PolyBase(const Eigen::MatrixXf& P0, const Eigen::MatrixXf& P1, const Fundamental& F,const Eigen::MatrixXf& K);
			/**
			 *	\brief	Triangulates image points.
			 *	\param	p0	Point in the image of the first camera.
			 *	\param	p1	Corresponding point in the image of the second camera.
			  *	\param K
			 *	\return	Triangulated point.
			 */
		Eigen::Vector3f triangulate(const Eigen::Vector2f& p0, const Eigen::Vector2f& p1) const override;
	protected:
		typedef Eigen::Vector3f Line;
		typedef Eigen::Vector3f Epipole;
		typedef std::vector<std::complex<float> > Roots;
		typedef std::tuple<float, float, float, float, float, float> PolyParams;
		/**
		 *	\brief	Computes corrected correspondences, that minimize the geometric error.
		 *	\param	p0	Point in the image of the first camera.
		 *	\param	p1	Corresponding point in the image of the second camera.
		 *	\return	Corrected correspondences in homogeneous coordinates.
		 */
		std::pair<Eigen::Vector2f, Eigen::Vector2f> ComputeCorrectedCorrespondences(const Eigen::Vector2f& p0, const Eigen::Vector2f& p1) const;
		/**
		 *	\brief	Defines translation matrix that translate the given point to origin.
		 *	\param	p	Translated point.
		 *	\return	Translation matrix.
		 */
		Eigen::MatrixXf TranslateToOrigin(const Eigen::Vector2f& p) const;
		/**
		 *	\brief	Computes epipole e = (e1, e2, e3) such as eF = 0 and e1*e1 + e2*e2 = 1
		 *	\param	F	Fundamental matrix.
		 *	\return	Computed epipole
		 */
		Epipole ComputeLeftEpipole(const Eigen::MatrixXf& F) const;
		/**
		 *	\brief	Computes epipole e = (e1, e2, e3) such as Fe = 0 and e1*e1 + e2*e2 = 1
		 *	\param	F	Fundamental matrix.
		 *	\return	Computed epipole
		 */
		Epipole ComputeRightEpipole(const Eigen::MatrixXf& F) const;
		/**
		 *	\brief	Defines rotation matrix using given epipole.
		 *	\param	e	Epipole.
		 *	\return	Rotation matrix.
		 */
		Eigen::MatrixXf FormRotationMatrix(const Epipole& e) const;
		/**
		 *	\brief	Prepares polynomial coefficients.
		 *	\param	params	Polynomial coefficients params.
		 *	\return Polynomial coefficients.
		 */
		virtual std::vector<float> PreparePolyCoeffs(const PolyParams& params) const = 0;
		/**
		 *	\brief	Forms and solves 6th degree polynomial.
		 *	\param	params	Polynomial coefficients params.
		 *	\return Six roots.
		 */
		Roots Solve(const PolyParams& params) const;
		/**
		 *	\brief	Evaluates cost function for each root.
		 *	\param	roots	Six roots of 6th degree polynomial.
		 *	\param	params	Polynomial coefficients params.
		 *	\return	Array of cost for each root.
		 */
		virtual std::vector<float> EvaluateRootsCosts(const Roots& roots, const PolyParams& params) const = 0;
		/**
		 *	\brief	Evaluates cost function for roots.
		 *	\param	roots	Six roots of 6th degree polynomial.
		 *	\param	params	Polynomial coefficients params.
		 *	\return	Root that gives smallest cost function value.
		 */
		float EvaluateRoots(const Roots& roots, const PolyParams& params) const;
		/**
		 *	\brief	Construct two epipolar lines on which the corrected correspondences lie.
		 *	\param	t	Real root for which the cost function gives the smallest value.
		 *	\param	params	Polynomial coefficients params.
		 *	\return	Pair of epipolar lines.
		 */
		std::pair<Line, Line> ConstructLines(float t, const PolyParams& params) const;
		/**
		 *	\brief	Finds point on given line that is closest to the origin.
		 *	\param	l	Line on which the point lies.
		 *	\return	Point on line closest to the origin.
		 */
		Eigen::Vector3f FindPointOnLineClosestToOrigin(const Line& l) const;
		/**
		 *	\brief	Transfers point to original coordinates using R and T.
		 *	\param	p	A point to transfer.
		 *	\param	R	Rotational transformation.
		 *	\param	T	Translational transformation.
		 *	\return 	Point in original coordinates.
		 */
		Eigen::Vector3f TransferPointToOriginalCoordinates(const Eigen::Vector3f& p, const Eigen::MatrixXf& R, const Eigen::MatrixXf& T) const;
		/**
		 *	\brief	Sets origin of world coordinate system to first camera.
		 *	\param	P0	First camera projection matrix.
		 *	\param	P1	Second camera projection matrix.
		 *	\return	K0, K1, R and T - camera intrinsics and orientation of second camera in new world coordinates.
		 */

		std::tuple<Eigen::MatrixXf, Eigen::MatrixXf> SetOriginToCamera(const Eigen::MatrixXf& P0, const Eigen::MatrixXf& P1, const Eigen::MatrixXf& K) const;		/**
		 *	\brief	Computes fundamental matrix 	from camera projection matrices.
		 *	\param	P0	First camera projection matrix.
		 *	\param	P1	Second camera projection matrix.
		  *	\param	K	Second camera projection matrix.
		 *	\return	Computed fundamental matrix.
		 */
		Fundamental ComputeFundamentalMatrix(const Eigen::MatrixXf& P0, const Eigen::MatrixXf& P1,const Eigen::MatrixXf& K) const;
		/**
		 *	\brief	Returns the order of the polynomial with given coefficients (highest non-zero coeffs index).
		 *	\param	coeffs	Polynomial coefficients.
		 *	\return	Polynomial order
		 */
		Eigen::Matrix3f getIntrinsicMatrix(const Eigen::MatrixXf& P) const;
		int FindPolynomialOrder(const std::vector<float>& coeffs) const;
		/**
		 *	\brief	Returns canonic camera projection matrix of second camera computed from given fundamental matrix.
		 *	\param	F	Fundamental matrix.
		 *	\return	Canonic	camera projection matrix.
		 */
		Eigen::MatrixXf CameraProjectionMatrixFromFundamentalMatrix(const Fundamental& F) const;


		const Fundamental F;
		/// LinearLS method used for final triangulation of corrected points.
		const LinearLS LS;
	};
}