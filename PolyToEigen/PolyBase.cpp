#include "PolyBase.h"
#include <Eigen/Dense>
#include <Eigen/QR>
#include <Eigen/Core>





namespace Triangulation {



	PolyBase::PolyBase(const PolyBase::Fundamental& F)
		: TriangulationBase(Eigen::Matrix<float, 3, 4>::Identity(), CameraProjectionMatrixFromFundamentalMatrix(F)),
		F(F), LS(P0, P1)
	{}


	PolyBase::PolyBase(const Eigen::MatrixXf& P0, const Eigen::MatrixXf& P1)
		: TriangulationBase(P0, P1), F(ComputeFundamentalMatrix(P0, P1)), LS(P0, P1)
	
	{}

	PolyBase::PolyBase(const Eigen::MatrixXf& P0, const Eigen::MatrixXf& P1, const PolyBase::Fundamental& F)
		: TriangulationBase(P0, P1), F(F), LS(P0, P1)
	{}

void PolyBase::decomposeProjectionMatrix(const Eigen::MatrixXf&P, Eigen::MatrixXf& K, Eigen::MatrixXf& R,  Eigen::MatrixXf&T) const
{
   MatrixXf Q_qr;
	MatrixXf R_qr;
	MatrixXf H_inf3x3 = P.block<3,3>(0,0);
	MatrixXf H_inf3x3_inv = H_inf3x3.inverse();
	QRdecomposition(H_inf3x3_inv,Q_qr,R_qr);
	//std::cout << "==============================Decomposing Using My Code==============================" << std::endl;
    K = R_qr.inverse();
	K = K / K(2, 2);
	K = K.triangularView<Eigen::Upper>();
   //std::cout << "Estimated Camera Matrix\n" << K << std::endl;
    Eigen::MatrixXf rotationMatrix = Q_qr.inverse();
	rotationMatrix = rotationMatrix * -1;
   // std::cout << "Estimated Camera Rotation\n" << rotationMatrix << std::endl;
   // std::cout << "Estimated Camera Translation" << std::endl;
	Eigen::VectorXf h3x1(3);
	h3x1 << P(0, 3), P(1, 3), P(2, 3);


    //t=-R*C, Q.inv()=R
    Eigen::VectorXf translation = (-Q_qr.inverse() * (-H_inf3x3.inverse() * h3x1));

	//Eigen::Vector4f translation4 = translation.homogeneous();
	//translation4 /= translation4[3];
   // std::cout << translation << std::endl;
	R= rotationMatrix;
	T = translation;

}


std::tuple<Eigen::MatrixXf, Eigen::MatrixXf,Eigen::MatrixXf, Eigen::MatrixXf> PolyBase::SetOriginToCamera(const Eigen::MatrixXf& P0, const Eigen::MatrixXf& P1) const{
Eigen::MatrixXf K0;
Eigen::MatrixXf K1;
Eigen::MatrixXf R0;
Eigen::MatrixXf  R1;
Eigen::MatrixXf T0;
Eigen::MatrixXf T1;
decomposeProjectionMatrix(P0,K0,R0,T0);
decomposeProjectionMatrix(P1,K1,R1,T1);
	/*
std::cout << "++++++++Final R0+++++++++++++++++" << std::endl;
	std::cout << R0 << std::endl;
	std::cout << "++++++++Final T0+++++++++++++++++" << std::endl;
	std::cout << T0 << std::endl;
	*/
	

Eigen::Matrix<float, 4, 4> M = Eigen::Matrix<float, 4, 4>::Identity();
M.block<3, 3>(0, 0) = R0.inverse();
M(0, 3) = T0(0) ;
M(1, 3) = T0(1) ;
M(2, 3) = T0(2);

//std::cout << "++++++++Final M+++++++++++++++++" << std::endl;
	//std::cout << M << std::endl;
Eigen::MatrixXf tmp = K1.inverse() * P1 * M;

Eigen::MatrixXf R = tmp.block<3, 3>(0, 0);
Eigen::MatrixXf T = tmp.block<3, 1>(0, 3);
	/*
	std::cout << "++++++++Final R+++++++++++++++++" << std::endl;
	std::cout << R << std::endl;
	std::cout << "++++++++Final T+++++++++++++++++" << std::endl;
	std::cout << T << std::endl;
	*/
	

 return std::make_tuple(K0, K1, R, T);

}


	void PolyBase::QRdecomposition(const Eigen::MatrixXf& A, Eigen::MatrixXf& Q, Eigen::MatrixXf& R)  const
	{
		/*
		Credite to :https://ros-developer.com/2019/01/01/decomposing-projection-using-opencv-and-c/
		*/
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


	PolyBase::Fundamental PolyBase::ComputeFundamentalMatrix(const  Eigen::MatrixXf& P0, const  Eigen::MatrixXf& P1) const
	{
		Eigen::MatrixXf R, T;
		//Eigen::MatrixXf
		Intrinsic K0;
		Intrinsic K1;
		std::tie(K0, K1, R,T) = SetOriginToCamera(P0, P1);
		Eigen::MatrixXf A = Eigen::MatrixXf(K0) * R.transpose() * T;

		Eigen::Matrix3f C = Eigen::Matrix3f::Zero();
		C(0, 1) = -A(2);
		C(0, 2) = A(1);
		C(1, 0) = A(2);
		C(1, 2) = -A(0);
		C(2, 0) = -A(1);
		C(2, 1) = A(0);

		return K1.inverse().transpose() * R * K0.transpose() * C;

	}

	Eigen::Vector3f PolyBase::triangulate(const Eigen::Vector2f& p0, const Eigen::Vector2f& p1) const
	{
		Eigen::Vector2f x0, x1;
		std::tie(x0, x1) = ComputeCorrectedCorrespondences(p0, p1);
		if (x0.isZero() || x1.isZero())
		{
				 return Eigen::Vector3f();
		}
		return LS.triangulate(x0, x1);
	}

	std::pair<Eigen::Vector2f, Eigen::Vector2f> PolyBase::ComputeCorrectedCorrespondences(const Eigen::Vector2f& p0, const Eigen::Vector2f& p1) const
	{
		Eigen::MatrixXf T0 = TranslateToOrigin(p0);
		Eigen::MatrixXf T1 = TranslateToOrigin(p1);
       		std::cout << "+++++++F+++++++" << std::endl;
	  	std::cout << F << std::endl;

		Eigen::MatrixXf f = T1.transpose() * Eigen::MatrixXf(F) * T0;
        //std::cout << f << std::endl;
		Epipole e0 = ComputeRightEpipole(f);	
			Epipole e1 = ComputeLeftEpipole(f);

		Eigen::MatrixXf R0 = FormRotationMatrix(e0);
		Eigen::MatrixXf R1 = FormRotationMatrix(e1);

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

	Eigen::MatrixXf PolyBase::TranslateToOrigin(const Eigen::Vector2f& p) const
	{
		Eigen::MatrixXf result = Eigen::MatrixXf::Identity(3, 3);
		result(0, 2) = p.x();
		result(1, 2) = p.y();
		
		return result;
	}

	PolyBase::Epipole PolyBase::ComputeLeftEpipole(const Eigen::MatrixXf& F) const
	{
		Eigen::MatrixXf W, U, VT;
		Eigen::JacobiSVD<Eigen::MatrixXf> svd(F.transpose(), Eigen::ComputeFullU | Eigen::ComputeFullV);
		U = svd.matrixU();
		VT = svd.matrixV().transpose();
		W = svd.singularValues();
		return Epipole(VT.row(2));
	}

	PolyBase::Epipole PolyBase::ComputeRightEpipole(const Eigen::MatrixXf& F) const
	{
		Eigen::MatrixXf W, U, VT;
		Eigen::JacobiSVD<Eigen::MatrixXf> svd(F, Eigen::ComputeFullU | Eigen::ComputeFullV);
		U = svd.matrixU();
		VT = svd.matrixV().transpose();
		W = svd.singularValues();
		return Epipole(VT.row(2));
	}

	Eigen::MatrixXf PolyBase::FormRotationMatrix(const Epipole& e) const
	{
		Eigen::MatrixXf result = Eigen::MatrixXf::Identity(3, 3);
		result(0, 0) = e.x();
		result(0, 1) = e.y();
		result(1, 0) = -e.y();
		result(1, 1) = e.x();
		return result;
	}

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


std::vector<Eigen::Vector2f> solvePoly(const std::vector<float>& coeffs) {
    int degree = coeffs.size() - 1;

    if (degree <= 0) {
        std::cerr << "Polynomial degree must be greater than zero!" << std::endl;
        return std::vector<Eigen::Vector2f>();
    }

    Eigen::MatrixXf A(degree, degree);
    Eigen::VectorXf b(degree);

    // Construct the companion matrix
	std::cout << "new solver" << std::endl;
    for (int i = 0; i < degree; i++) {
            A(i, 0) = -coeffs[i + 1] / coeffs[0];
        if (i > 0) {
            A(i, i) = 1;
        }
    }

    // Solve the eigenvalue problem
    Eigen::EigenSolver<Eigen::MatrixXf> eigensolver(A);
	
    if (eigensolver.info() != Eigen::Success) {
        std::cerr << "Failed to solve the polynomial equation!" << std::endl;
        return std::vector<Eigen::Vector2f>();
    }

    // Extract the eigenvalues as roots
    std::vector<Eigen::Vector2f> roots(degree);
    for (int i = 0; i < degree; i++) {
        const std::complex<float>& eigenvalue = eigensolver.eigenvalues()[i];
        roots[i] = Eigen::Vector2f(eigenvalue.real(), eigenvalue.imag());
    }

    return roots;
}


PolyBase::Roots PolyBase::Solve(const PolyParams& params) const
{
	
    std::vector<float> coeffs = PreparePolyCoeffs(params);
    if (coeffs.size() <= 1)
    {
        return { 0 };
    }
	

    // Solve the polynomial equation and obtain the roots

	std::vector<Eigen::Vector2f> roots =solvePoly(coeffs);

	Roots result(roots.size());
	for (size_t i = 0u; i < roots.size(); ++i)
	{
		result[i] = std::complex<float>(roots[i][0], roots[i][1]);
	}
	return result;

    //Eigen::VectorXf polynomial(coeffs.size());
    //for (size_t i = 0u; i < coeffs.size(); ++i)
    //{
     //   polynomial(i) = coeffs[i];
    //}

    //Eigen::EigenSolver<Eigen::MatrixXf> solver(polynomial.reverse().asDiagonal());

   // PolyBase::Roots result(solver.eigenvalues().size());
    //for (size_t i = 0u; i < solver.eigenvalues().size(); ++i)
    //{
      //  result[i] = solver.eigenvalues()[i].real();
    //}

    //return result;
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
		Line l0(t * e, 1, -t);
		Line l1(-f * (c * t + d), a * t + b, c * t + d);
		return std::make_pair(l0, l1);
	}

	Eigen::Vector3f PolyBase::FindPointOnLineClosestToOrigin(const PolyBase::Line& l) const
	{
		return Eigen::Vector3f(-l[0] * l[2], -l[1] * l[2], l[0] * l[0] + l[1] * l[1]);
	}

	Eigen::Vector3f PolyBase::TransferPointToOriginalCoordinates(const Eigen::Vector3f& p, const Eigen::MatrixXf& R, const Eigen::MatrixXf& T) const
	{
		Eigen::MatrixXf x = Eigen::MatrixXf::Identity(3, 1);
		x(0) = p.x();
		x(1) = p.y();
		x(2) = p.z();
		x = T * R.transpose() * x;
		return Eigen::Vector3f(x(0), x(1), x(2));
	}

	Eigen::MatrixXf PolyBase::CameraProjectionMatrixFromFundamentalMatrix(const PolyBase::Fundamental& F) const
	{ //not sure at all !!!!!!!!!!!!!!! maybe all what i did here is a mistake
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
