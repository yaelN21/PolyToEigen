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
 	   ~PolyBase(); 
		PolyBase(const Eigen::Matrix<float, 3, 4>& P0,const Eigen::Matrix<float, 3, 4>& P1,
		const Eigen::Matrix3f& K0,const Eigen::Matrix3f& K1,const Eigen::Matrix3f& R0,const Eigen::Matrix3f& R1,
		const  Eigen::Vector3f& T0,const  Eigen::Vector3f& T1);
			/**
			 *	\brief	Constructor
			 *	\param	P0	Camera matrix of the first camera.
			 *	\param 	P1	Camera matrix of the second camera.
			  *	\param K
			 *	\param	F	Fundamental matrix.
			 */

	
		explicit PolyBase(const Fundamental& F);
		/**
		 *	\brief	Constructor
		 *	\param	P0	Camera matrix of the first camera.
		 *	\param 	P1	Camera matrix of the second camera.
		 *\param 	K
		 */
		
		PolyBase(
		const Eigen::Matrix<float, 3, 4>& P0,const Eigen::Matrix<float, 3, 4>& P1,
		const Eigen::Matrix3f& K0,const Eigen::Matrix3f& K1,const Eigen::Matrix3f& R0,const Eigen::Matrix3f& R1,
		const  Eigen::Vector3f& T0,const  Eigen::Vector3f& T1, const PolyBase::Fundamental& F);
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
		std::pair<Eigen::Vector2f, Eigen::Vector2f> ComputeCorrectedCorrespondences(const Eigen::Vector2f& p0, const Eigen::Vector2f& p1) const;
		/**
		 *	\brief	Computes corrected correspondences, that minimize the geometric error.
		 *	\param	p0	Point in the image of the first camera.
		 *	\param	p1	Corresponding point in the image of the second camera.
		 *	\return	Corrected correspondences in homogeneous coordinates.
		 */
		Eigen::Matrix3f TranslateToOrigin(const Eigen::Vector2f& p) const;
		/**
		 *	\brief	Defines translation matrix that translate the given point to origin.
		 *	\param	p	Translated point.
		 *	\return	Translation matrix.
		 */
		Epipole ComputeLeftEpipole(const Eigen::Matrix3f& F) const;
		/**
		 *	\brief	Computes epipole e = (e1, e2, e3) such as eF = 0 and e1*e1 + e2*e2 = 1
		 *	\param	F	Fundamental matrix.
		 *	\return	Computed epipole
		 */
		Epipole ComputeRightEpipole(const Eigen::Matrix3f& F) const;
		/**
		 *	\brief	Computes epipole e = (e1, e2, e3) such as Fe = 0 and e1*e1 + e2*e2 = 1
		 *	\param	F	Fundamental matrix.
		 *	\return	Computed epipole
		 */
		Eigen::Matrix3f FormRotationMatrix(const Epipole& e) const;
		/**
		 *	\brief	Defines rotation matrix using given epipole.
		 *	\param	e	Epipole.
		 *	\return	Rotation matrix.
		 */
		virtual std::vector<float> PreparePolyCoeffs(const PolyParams& params) const = 0;
		/**
		 *	\brief	Prepares polynomial coefficients.
		 *	\param	params	Polynomial coefficients params.
		 *	\return Polynomial coefficients.
		 */
		Roots Solve(const PolyParams& params) const;
		/**
		 *	\brief	Forms and solves 6th degree polynomial.
		 *	\param	params	Polynomial coefficients params.
		 *	\return Six roots.
		 */
		virtual std::vector<float> EvaluateRootsCosts(const Roots& roots, const PolyParams& params) const = 0;
		/**
		 *	\brief	Evaluates cost function for each root.
		 *	\param	roots	Six roots of 6th degree polynomial.
		 *	\param	params	Polynomial coefficients params.
		 *	\return	Array of cost for each root.
		 */
		float EvaluateRoots(const Roots& roots, const PolyParams& params) const;
		/**
		 *	\brief	Evaluates cost function for roots.
		 *	\param	roots	Six roots of 6th degree polynomial.
		 *	\param	params	Polynomial coefficients params.
		 *	\return	Root that gives smallest cost function value.
		 */
		std::pair<Line, Line> ConstructLines(float t, const PolyParams& params) const;
		/**
		 *	\brief	Construct two epipolar lines on which the corrected correspondences lie.
		 *	\param	t	Real root for which the cost function gives the smallest value.
		 *	\param	params	Polynomial coefficients params.
		 *	\return	Pair of epipolar lines.
		 */
		Eigen::Vector3f FindPointOnLineClosestToOrigin(const Line& l) const;
		/**
		 *	\brief	Finds point on given line that is closest to the origin.
		 *	\param	l	Line on which the point lies.
		 *	\return	Point on line closest to the origin.
		 */
		Eigen::Vector3f TransferPointToOriginalCoordinates(const Eigen::Vector3f& p, const Eigen::Matrix3f& R, const Eigen::Matrix3f& T) const;
		/**
		 *	\brief	Transfers point to original coordinates using R and T.
		 *	\param	p	A point to transfer.
		 *	\param	R	Rotational transformation.
		 *	\param	T	Translational transformation.
		 *	\return 	Point in original coordinates.
		 */
		std::tuple<Eigen::Matrix3f, Eigen::Vector3f> SetOriginToCamera( const Eigen::Matrix<float, 3, 4>& P0,const Eigen::Matrix<float, 3, 4>& P1) const;
		/**
		 *	\brief	Sets origin of world coordinate system to first camera.
		 *	\param	P0	First camera projection matrix.
		 *	\param	P1	Second camera projection matrix.
		 *	\return	K0, K1, R and T - camera intrinsics and orientation of second camera in new world coordinates.
		 */
		Fundamental ComputeFundamentalMatrix( const Eigen::Matrix<float, 3, 4>& P0,const Eigen::Matrix<float, 3, 4>& P1) const;
	
		int FindPolynomialOrder(const std::vector<float>& coeffs) const;
			/**
		 *	\brief	Returns the order of the polynomial with given coefficients (highest non-zero coeffs index).
		 *	\param	coeffs	Polynomial coefficients.
		 *	\return	Polynomial order
		 */
		Eigen::MatrixXf CameraProjectionMatrixFromFundamentalMatrix(const Fundamental& F) const;
		/**
		 *	\brief	Returns canonic camera projection matrix of second camera computed from given fundamental matrix.
		 *	\param	F	Fundamental matrix.
		 *	\return	Canonic	camera projection matrix.
		 */
		const Fundamental F;
		/// LinearLS method used for final triangulation of corrected points.
		const LinearLS LS;
		
	};

}