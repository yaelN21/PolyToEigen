#include "PolyBase.h"
#include <Eigen/Dense>
#include <Eigen/QR>
#include <Eigen/Core>





namespace Triangulation {



	PolyBase::PolyBase(const PolyBase::Fundamental& F,const Eigen::MatrixXf& K)
		: TriangulationBase(Eigen::Matrix<float, 3, 4>::Identity(), CameraProjectionMatrixFromFundamentalMatrix(F),K),
		F(F), LS(P0, P1,K)
	{}


	PolyBase::PolyBase(const Eigen::MatrixXf& P0, const Eigen::MatrixXf& P1,const Eigen::MatrixXf& K)
		: TriangulationBase(P0, P1,K), F(ComputeFundamentalMatrix(P0, P1,K)), LS(P0, P1,K)
	
	{}

	PolyBase::PolyBase(const Eigen::MatrixXf& P0, const Eigen::MatrixXf& P1, const PolyBase::Fundamental& F,const Eigen::MatrixXf& K)
		: TriangulationBase(P0, P1,K), F(F), LS(P0, P1,K)
	{}

std::tuple<Eigen::MatrixXf, Eigen::MatrixXf> PolyBase::SetOriginToCamera(const Eigen::MatrixXf& P0, const Eigen::MatrixXf& P1, const Eigen::MatrixXf& K) const{
	//Eigen::Matrix3f R0;
	Eigen::Matrix3f R1;
	Eigen::Vector4f T0;
	Eigen::Vector4f T1;
 
	
	
	Eigen::HouseholderQR<Eigen::MatrixXf> qr0(P0);
    //K0 = qr0.matrixQR().topLeftCorner<3, 3>().triangularView<Eigen::Upper>();
	
	//Eigen::ColPivHouseholderQR<Eigen::MatrixXf> qr1(P1);
    //K1 = qr1.matrixQR().topLeftCorner<3, 3>().triangularView<Eigen::Upper>();
	//K1 = getIntrinsicMatrix(P1);
	//std::cout << K1 << std::endl;

	Eigen::HouseholderQR<Eigen::MatrixXf> qr(P0.block<3,3>(0,0).inverse());
	Eigen::MatrixXf Q0 = qr.householderQ();
	Eigen::MatrixXf R0 = qr.matrixQR().triangularView<Eigen::Upper>();
	std::cout << Q0 << std::endl;
    std::cout << "++++++++++++++++++++++++++++++++" << std::endl;
	std::cout << Q0.transpose()*-1 << std::endl;
	std::cout << "++++++++++++++++++++++++++++++++" << std::endl;
	R0 = P0.block<3, 3>(0, 0);
	R1 = P1.block<3, 3>(0, 0);
	std::cout << Q0 << std::endl;
	T0 << P0(0, 3), P0(1, 3), P0(2, 3), P0(2, 3);
	T1 << P1(0, 3), P1(1, 3), P1(2, 3), P1(2, 3);
	//std::cout << R0 << std::endl;
	std::cout << "++++++++++++++++++++++++" << std::endl;
	//std::cout << T0 << std::endl;

    Eigen::Matrix4f M = Eigen::Matrix4f::Identity();
    M.block<3, 3>(0, 0) = R0.inverse();

    M(0, 3) = (T0(0) /T0(3));
    M(1, 3) = (T0(1) / T0(3));
    M(2, 3) = (T0(2) / T0(3));

    // K0.inv() * P0 * M - should be identity
    Eigen::Matrix<float, 3, 4> tmp = K.inverse() * P1 * M;

    Eigen::Matrix3f R = tmp.block<3, 3>(0, 0);
    Eigen::Vector3f T = tmp.block<3, 1>(0, 3);
	//std::cout << T << std::endl;
	std::cout << "+++++++++++++++++++++++++" << std::endl;
	//std::cout << R << std::endl;
    return std::make_tuple( R, T);
}
/*
 void PolyBase::getIntrinsicMatrix(const Eigen::MatrixXf &A, Eigen::MatrixXf &Q, Eigen::MatrixXf &R) const
{
 
    auto sign = [](float value) { return value >= 0 ? 1: -1; };
    const auto totalRows = A.rows;
    const auto totalCols = A.cols;
    R = A.clone();
    Q = Eigen::MatrixXf::eye ( totalRows, totalRows, A.type() );
    for ( int col = 0; col < A.cols; ++ col )
    {
        Eigen::MatrixXf matAROI = Eigen::MatrixXf ( R, cv::Range ( col, totalRows ), cv::Range ( col, totalCols ) );
        Eigen::MatrixXf y = matAROI.col ( 0 );
        auto yNorm = norm ( y );
        Eigen::MatrixXf e1 = Eigen::MatrixXf::eye ( y.rows, 1, A.type() );
        Eigen::MatrixXf w = y + sign(y.at<float>(0,0)) *  yNorm * e1;
        Eigen::MatrixXf v = w / norm( w );
        Eigen::MatrixXf vT; cv::transpose(v, vT );
        Eigen::MatrixXf I = Eigen::MatrixXf::eye( matAROI.rows, matAROI.rows, A.type() );
        Eigen::MatrixXf I_2VVT = I - 2 * v * vT;
        Eigen::MatrixXf matH = Eigen::MatrixXf::eye ( totalRows, totalRows, A.type() );
        Eigen::MatrixXf matHROI = Eigen::MatrixXf(matH, cv::Range ( col, totalRows ), cv::Range ( col, totalRows ) );
        I_2VVT.copyTo ( matHROI );
        R = matH * R;
        Q = Q * matH;
    }

}
*/

	void PolyBase::getIntrinsicMatrix(const Eigen::MatrixXf& A, Eigen::MatrixXf& Q, Eigen::MatrixXf& R)  const
	{
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
			Eigen::MatrixXf w = y + sign(y(0)) * yNorm * e1;
			Eigen::MatrixXf v = w.normalized();
			Eigen::MatrixXf vT = v.transpose();
			Eigen::MatrixXf I = Eigen::MatrixXf::Identity(matAROI.rows(), matAROI.rows());
			Eigen::MatrixXf I_2VVT = I - 2 * v * vT;
			Eigen::MatrixXf matH = Eigen::MatrixXf::Identity(totalRows, totalRows);
			Eigen::MatrixXf matHROI = matH.block(col, col, totalRows - col, totalRows - col);
			matHROI = I_2VVT;
			R = matH * R;
			Q = Q * matH;
		}
	}


	PolyBase::Fundamental PolyBase::ComputeFundamentalMatrix(const  Eigen::MatrixXf& P0, const  Eigen::MatrixXf& P1, const Eigen::MatrixXf& K) const
	{
		Eigen::MatrixXf R, T;
		Intrinsic K0 = K;
		Intrinsic K1 = K;
		std::tie(R, T) = SetOriginToCamera(P0, P1,K);

		Eigen::MatrixXf A = Eigen::MatrixXf(K0) * R.transpose() * T;

		Eigen::Matrix3f C = Eigen::Matrix3f::Zero();
		C(0, 1) = -A(2);
		C(0, 2) = A(1);
		C(1, 0) = A(2);
		C(1, 2) = -A(0);
		C(2, 0) = -A(1);
		C(2, 1) = A(0);

		return K1.transpose() * R * K0.transpose() * C;
	}

	Eigen::Vector3f PolyBase::triangulate(const Eigen::Vector2f& p0, const Eigen::Vector2f& p1) const
	{
		Eigen::Vector2f x0, x1;
		std::tie(x0, x1) = ComputeCorrectedCorrespondences(p0, p1);
		//std::cout << x0 << std::endl;
		//std::cout << x1 << std::endl;
		return LS.triangulate(x0, x1);
	}

	std::pair<Eigen::Vector2f, Eigen::Vector2f> PolyBase::ComputeCorrectedCorrespondences(const Eigen::Vector2f& p0, const Eigen::Vector2f& p1) const
	{
		Eigen::MatrixXf T0 = TranslateToOrigin(p0);
		Eigen::MatrixXf T1 = TranslateToOrigin(p1);

		Eigen::MatrixXf f = T1.transpose() * Eigen::MatrixXf(F) * T0;

		Epipole e0 = ComputeRightEpipole(f);
		Epipole e1 = ComputeLeftEpipole(f);

		Eigen::MatrixXf R0 = FormRotationMatrix(e0);
		Eigen::MatrixXf R1 = FormRotationMatrix(e1);

		f = R1 * f * R0.transpose();

		PolyParams params = std::make_tuple(f(1, 1), f(1, 2), f(2, 1), f(2, 2), e0.z(), e1.z());
		Roots roots = Solve(params);

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
