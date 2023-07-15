#include "../include/PolyBase.h"
#include <Eigen/Dense>
#include <Eigen/QR>
#include <Eigen/Core>
#include <opencv2/opencv.hpp>




namespace Triangulation {
	//constructor for linear triangulation
	PolyBase::PolyBase(const PolyBase::Fundamental& F)
		: TriangulationBase(Eigen::Matrix<float, 3, 4>::Identity(), CameraProjectionMatrixFromFundamentalMatrix(F)),
		F(F), LS(P0, P1)
	{}

	//consturctor for poly triangulation - unkown fundamental matrix
	PolyBase::PolyBase(const Eigen::Matrix<float, 3, 4>& P0,const Eigen::Matrix<float, 3, 4>& P1,
		const Eigen::Matrix3f& K0,const Eigen::Matrix3f& K1,const Eigen::Matrix3f& R0,const Eigen::Matrix3f& R1,
		const  Eigen::Vector3f& T0,const  Eigen::Vector3f& T1)
		: TriangulationBase(P0,P1,K0, K1,R0, R1,T0, T1), F(ComputeFundamentalMatrix(P0,P1)), LS(P0, P1)
	
	{}
	
	//consturctor for poly triangulation - kown fundamental matrix
	PolyBase::PolyBase(const Eigen::Matrix<float, 3, 4>& P0,const Eigen::Matrix<float, 3, 4>& P1,
		const Eigen::Matrix3f& K0,const Eigen::Matrix3f& K1,const Eigen::Matrix3f& R0,const Eigen::Matrix3f& R1,
		const  Eigen::Vector3f& T0,const  Eigen::Vector3f& T1, const PolyBase::Fundamental& F)
		: TriangulationBase(P0,P1,K0, K1,R0, R1,T0, T1), F(F), LS(P0, P1)
	{}
 PolyBase::~PolyBase() {
        // Destructor implementation
        // Add any necessary cleanup code here
    }
   	/*
	void PolyBase::QRdecomposition(const Eigen::MatrixXf& A, Eigen::MatrixXf& Q, Eigen::MatrixXf& R)  const
	{
	
		//Credite to :https://ros-developer.com/2019/01/01/decomposing-projection-using-opencv-and-c/
		 //assert(A.channels() == 1);
    	assert(A.rows() >= A.cols());
   		 std::cout << "Assertions passed successfully." << std::endl;
		auto sign = [](float value) { return value >= 0 ? 1 : -1; };
		const auto totalRows = A.rows();
		const auto totalCols = A.cols();
		R = A;
		Q = Eigen::MatrixXf::Identity(totalRows, totalRows);
		for (int col = 0; col < A.cols(); ++col)
		{
			Eigen::MatrixXf matAROI = R.block(col, col, totalRows - col, totalCols - col);
			Eigen::MatrixXf y = matAROI.col(0);
			auto yNorm = y.norm();
			Eigen::MatrixXf e1 = Eigen::MatrixXf::Identity(y.rows(), 1);
			//Eigen::MatrixXf tmp1 =  yNorm * e1;
			Eigen::MatrixXf tmp1 = sign(y(0)) * yNorm * e1;
			Eigen::MatrixXf w = y +tmp1 ;
			Eigen::MatrixXf v = w.normalized();
			Eigen::MatrixXf vT = v.transpose();
			Eigen::MatrixXf I = Eigen::MatrixXf::Identity(matAROI.rows(), matAROI.rows());
			Eigen::MatrixXf I_2VVT = I - 2 * v * vT;
			Eigen::MatrixXf matH = Eigen::MatrixXf::Identity(totalRows, totalRows);
			Eigen::MatrixXf matHROI = matH.block(col, col, totalRows - col, totalRows - col);
			matHROI = I_2VVT;
			matH.block(col, col, matHROI.rows(), matHROI.cols()) = matHROI;
			R = matH * R;
			Q = Q * matH;
		}

	}
*/
/*
void PolyBase::decomposeProjectionMatrix(const Eigen::MatrixXf&P, Eigen::MatrixXf& K, Eigen::MatrixXf& R,  Eigen::MatrixXf&T) const
{
   MatrixXf Q_qr;
	MatrixXf R_qr;
	MatrixXf H_inf3x3 = P.block<3,3>(0,0);
	MatrixXf H_inf3x3_inv = H_inf3x3.inverse();
	QRdecomposition(H_inf3x3_inv,Q_qr,R_qr);
    K = R_qr.inverse();
	K = K / K(2, 2);
	K = K.triangularView<Eigen::Upper>();
    Eigen::MatrixXf rotationMatrix = Q_qr.inverse();
	rotationMatrix = rotationMatrix * -1;
	Eigen::VectorXf h3x1(3);
	h3x1 << P(0, 3), P(1, 3), P(2, 3);
    //t=-R*C, Q.inv()=R
    Eigen::VectorXf translation = (-Q_qr.inverse() * (-H_inf3x3.inverse() * h3x1));
	R= rotationMatrix;
	T = translation;

}
*/
std::tuple<Eigen::Matrix3f, Eigen::Vector3f> PolyBase::SetOriginToCamera( const Eigen::Matrix<float, 3, 4>& P0,const Eigen::Matrix<float, 3, 4>& P1) const{

	Eigen::Matrix<float, 4, 4> M = Eigen::Matrix<float, 4, 4>::Identity();
	M.block<3, 3>(0, 0) = R0.inverse();
	M(0, 3) = T0(0);
	M(1, 3) = T0(1) ;
	M(2, 3) = T0(2);

	Eigen::Matrix<float, 3, 4>  tmp = K1.inverse() * P1 * M;

	Eigen::Matrix3f R = tmp.block<3, 3>(0, 0);
	Eigen::Vector3f T = tmp.block<3, 1>(0, 3);

	 return std::make_tuple(R, T);
}


	PolyBase::Fundamental PolyBase::ComputeFundamentalMatrix( const Eigen::Matrix<float, 3, 4>& P0,const Eigen::Matrix<float, 3, 4>& P1) const
	{
		Eigen::Matrix3f R;
		Eigen::Vector3f T;

		std::tie(R,T) = SetOriginToCamera(P0,P1);
		/*
		Eigen::Vector3f A = K0 * R.transpose() * T; 
		Eigen::Matrix3f C = Eigen::Matrix3f::Zero();
		C(0, 1) = -A(2);
		C(0, 2) = A(1);
		C(1, 0) = A(2);
		C(1, 2) = -A(0);
		C(2, 0) = -A(1);
		C(2, 1) = A(0);
		Fundamental f= K1.inverse().transpose() * R * K0.transpose() * C;
		return f;
  		*/
		Eigen::Matrix3f skewSymmetric_T = Eigen::Matrix3f::Zero(3, 3);
		skewSymmetric_T(0, 1) = -T(2);
		skewSymmetric_T(0, 2) = T(1);
		skewSymmetric_T(1, 0) = T(2);
		skewSymmetric_T(1, 2) = -T(0);
		skewSymmetric_T(2, 0) = -T(1);
		skewSymmetric_T(2, 1) = T(0);
		Fundamental f = K1.inverse().transpose() * skewSymmetric_T * R * K0.inverse();
		return f;
		

	}

	Eigen::Vector3f PolyBase::triangulate(const Eigen::Vector2f& p0, const Eigen::Vector2f& p1) const
	{
		Eigen::Vector2f x0, x1;
		std::tie(x0, x1) = ComputeCorrectedCorrespondences(p0, p1);
		if (x0.isZero() || x1.isZero())
		{
				   Eigen::Vector3f zeroVector(0.0f, 0.0f,0.0f);
			    return zeroVector;
		}
		return LS.triangulate(x0, x1);
	}

	std::pair<Eigen::Vector2f, Eigen::Vector2f> PolyBase::ComputeCorrectedCorrespondences(const Eigen::Vector2f& p0, const Eigen::Vector2f& p1) const
	{
		
		Eigen::Matrix3f T0 = TranslateToOrigin(p0);
		Eigen::Matrix3f T1 = TranslateToOrigin(p1);

		Eigen::Matrix3f f = T1.transpose() * F * T0;
		
		Epipole e0 = ComputeRightEpipole(f);	
		Epipole e1 = ComputeLeftEpipole(f);

		Eigen::Matrix3f R0 = FormRotationMatrix(e0);
		Eigen::Matrix3f R1 = FormRotationMatrix(e1);

		f = R1 * f * R0.transpose();

		PolyParams params = std::make_tuple(f(1, 1), f(1, 2), f(2, 1), f(2, 2), e0.z(), e1.z());
		Roots roots = Solve(params);
		if (roots.empty())
		{
			   Eigen::Vector2f zeroVector(0.0f, 0.0f);
			   return std::make_pair(zeroVector, zeroVector);
		}

		float t = EvaluateRoots(roots, params);
		Line l0, l1;
		std::tie(l0, l1) = ConstructLines(t, params);

		Eigen::Vector3f x0 = FindPointOnLineClosestToOrigin(l0);
		Eigen::Vector3f x1 = FindPointOnLineClosestToOrigin(l1);

		x0 = TransferPointToOriginalCoordinates(x0, R0, T0);
		x1 = TransferPointToOriginalCoordinates(x1, R1, T1);

		return std::make_pair(Eigen::Vector2f(x0.x() / x0.z(), x0.y() / x0.z()), Eigen::Vector2f(x1.x() / x1.z(), x1.y() / x1.z()));
	}

	Eigen::Matrix3f PolyBase::TranslateToOrigin(const Eigen::Vector2f& p) const
	{
		Eigen::Matrix3f result = Eigen::Matrix3f::Identity(3, 3);
		result(0, 2) = p.x();
		result(1, 2) = p.y();
		
		return result;
	}

	PolyBase::Epipole PolyBase::ComputeLeftEpipole(const Eigen::Matrix3f& F) const
	{
		Eigen::MatrixXf W, U, VT;
		Eigen::JacobiSVD<Eigen::MatrixXf> svd(F.transpose(), Eigen::ComputeFullU | Eigen::ComputeFullV);
		U = svd.matrixU();
		VT = svd.matrixV().transpose();
		W = svd.singularValues();
		return Epipole(VT.row(2));
	}

	PolyBase::Epipole PolyBase::ComputeRightEpipole(const Eigen::Matrix3f& F) const
	{
		Eigen::MatrixXf W, U, VT;
		Eigen::JacobiSVD<Eigen::MatrixXf> svd(F, Eigen::ComputeFullU | Eigen::ComputeFullV);
		U = svd.matrixU();
		VT = svd.matrixV().transpose();
		W = svd.singularValues();
		return Epipole(VT.row(2));
	}

	Eigen::Matrix3f PolyBase::FormRotationMatrix(const Epipole& e) const
	{
		Eigen::Matrix3f result = Eigen::Matrix3f::Identity(3, 3);
		result(0, 0) = e.x();
		result(0, 1) = e.y();
		result(1, 0) = -e.y();
		result(1, 1) = e.x();
		return result;
	}

	//NOT USED
	int PolyBase::FindPolynomialOrder(const std::vector<float>& coeffs) const
	{
		for (int n = static_cast<int>(coeffs.size()) - 1; n >= 0; --n)
		{
			if (coeffs[n] != 0)
			{
				return n;
			}
		}
		return -1;
	}



PolyBase::Roots PolyBase::Solve(const PolyParams& params) const
{
	
    std::vector<float> coeffs = PreparePolyCoeffs(params);
    if (coeffs.size() <= 1)
    {
        return { 0 };
    }
	

    // Solve the polynomial equation and obtain the roots
       std::vector<cv::Vec2f> roots;
	 try {
        cv::solvePoly(coeffs, roots);
		}
	catch (const cv::Exception& e) {
        std::cout << "The polynomial equation cannot be solved." << std::endl;
    }

	Roots result(roots.size());
	for (size_t i = 0u; i < roots.size(); ++i)
	{
		result[i] = std::complex<float>(roots[i][0], roots[i][1]);
	}
	return result;

}


	float PolyBase::EvaluateRoots(const PolyBase::Roots& roots, const PolyParams& params) const
	{
		std::vector<float> costs = EvaluateRootsCosts(roots, params);
		return roots[std::min_element(costs.begin(), costs.end()) - costs.begin()].real();
	}

	std::pair<PolyBase::Line, PolyBase::Line> PolyBase::ConstructLines(float t, const PolyBase::PolyParams& params) const
	{
		float a, b, c, d, e, f;
		std::tie(a, b, c, d, e, f) = params;
		//Line l0(t * e, 1, -t); mistake- should be f
		Line l0(t * f, 1, -t);
		Line l1(-f * (c * t + d), a * t + b, c * t + d);
		return std::make_pair(l0, l1);
	}

	Eigen::Vector3f PolyBase::FindPointOnLineClosestToOrigin(const PolyBase::Line& l) const
	{
		return Eigen::Vector3f(-l[0] * l[2], -l[1] * l[2], l[0] * l[0] + l[1] * l[1]);
	}

	Eigen::Vector3f PolyBase::TransferPointToOriginalCoordinates(const Eigen::Vector3f& p, const Eigen::Matrix3f& R, const Eigen::Matrix3f& T) const
	{

		Eigen::Vector3f x;
		x(0) = p.x();
		x(1) = p.y();
		x(2) = p.z();
		x = T * R.transpose() * x;
		return Eigen::Vector3f(x(0), x(1), x(2));
	}

	Eigen::MatrixXf PolyBase::CameraProjectionMatrixFromFundamentalMatrix(const PolyBase::Fundamental& F) const
	{ 
		Eigen::MatrixXf e2;
		Eigen::JacobiSVD<Eigen::MatrixXf> svd(F.transpose(), Eigen::ComputeThinU | Eigen::ComputeThinV);
		Eigen::MatrixXf Z = svd.solve(e2);
		Eigen::Matrix3f e2x;
		e2x << 0, -Z(2), Z(1),
			Z(2), 0, -Z(0),
			-Z(1), Z(0), 0;
		Eigen::MatrixXf result = Eigen::MatrixXf::Identity(3, 4);
		Eigen::MatrixXf mul = e2x * F;

		result.block(0, 0, 3, 3) = mul;
		result(0, 3) = Z(0);
		result(1, 3) = Z(1);
		result(2, 3) = Z(2);
		return result;
	}
}
