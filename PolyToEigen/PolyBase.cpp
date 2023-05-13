#include "PolyBase.h"
#include <Eigen/Dense>

namespace Triangulation {
	PolyBase::PolyBase(const PolyBase::Fundamental& F)
		: TriangulationBase(Eigen::Matrix<double, 3, 4>::Identity(), CameraProjectionMatrixFromFundamentalMatrix(F)),
		F(F), LS(P0, P1)
	{}


	PolyBase::PolyBase(const Eigen::MatrixXd& P0, const Eigen::MatrixXd& P1)
		: TriangulationBase(P0, P1), F(ComputeFundamentalMatrix(P0, P1)), LS(P0, P1)
	{}

	PolyBase::PolyBase(const Eigen::MatrixXd& P0, const Eigen::MatrixXd& P1, const PolyBase::Fundamental& F)
		: TriangulationBase(P0, P1), F(F), LS(P0, P1)
	{}

	std::tuple<PolyBase::Intrinsic, PolyBase::Intrinsic, Eigen::MatrixXd, Eigen::MatrixXd> PolyBase::SetOriginToCamera(const Eigen::MatrixXd& P0, const Eigen::MatrixXd& P1) const
	{

	}

	PolyBase::Fundamental PolyBase::ComputeFundamentalMatrix(const  Eigen::MatrixXd& P0, const  Eigen::MatrixXd& P1) const
	{
		Eigen::MatrixXd R, T;
		Intrinsic K0, K1;
		std::tie(K0, K1, R, T) = SetOriginToCamera(P0, P1);

		Eigen::MatrixXd A = Eigen::MatrixXd(K0) * R.transpose() * T;

		Eigen::Matrix3d C = Eigen::Matrix3d::Zero();
		C(0, 1) = -A(2);
		C(0, 2) = A(1);
		C(1, 0) = A(2);
		C(1, 2) = -A(0);
		C(2, 0) = -A(1);
		C(2, 1) = A(0);

		return K1.transpose() * R * K0.transpose() * C;
	}

	Eigen::Vector3d PolyBase::triangulate(const Eigen::Vector2d& p0, const Eigen::Vector2d& p1) const
	{
		Eigen::Vector2d x0, x1;
		std::tie(x0, x1) = ComputeCorrectedCorrespondences(p0, p1);
		//std::cout << x0 << std::endl;
		//std::cout << x1 << std::endl;
		return LS.triangulate(x0, x1);
	}

	std::pair<Eigen::Vector2d, Eigen::Vector2d> PolyBase::ComputeCorrectedCorrespondences(const Eigen::Vector2d& p0, const Eigen::Vector2d& p1) const
	{
		Eigen::MatrixXd T0 = TranslateToOrigin(p0);
		Eigen::MatrixXd T1 = TranslateToOrigin(p1);

		Eigen::MatrixXd f = T1.transpose() * Eigen::MatrixXd(F) * T0;

		Epipole e0 = ComputeRightEpipole(f);
		Epipole e1 = ComputeLeftEpipole(f);

		Eigen::MatrixXd R0 = FormRotationMatrix(e0);
		Eigen::MatrixXd R1 = FormRotationMatrix(e1);

		f = R1 * f * R0.transpose();

		PolyParams params = std::make_tuple(f(1, 1), f(1, 2), f(2, 1), f(2, 2), e0.z(), e1.z());
		Roots roots = Solve(params);

		double t = EvaluateRoots(roots, params);
		Line l0, l1;
		std::tie(l0, l1) = ConstructLines(t, params);

		Eigen::Vector3d x0 = FindPointOnLineClosestToOrigin(l0);
		Eigen::Vector3d x1 = FindPointOnLineClosestToOrigin(l1);

		x0 = TransferPointToOriginalCoordinates(x0, R0, T0);
		x1 = TransferPointToOriginalCoordinates(x1, R1, T1);

		return std::make_pair(Eigen::Vector2d(x0.x() / x0.z(), x0.y() / x0.z()), Eigen::Vector2d(x1.x() / x1.z(), x1.y() / x1.z()));
	}

	Eigen::MatrixXd PolyBase::TranslateToOrigin(const Eigen::Vector2d& p) const
	{
		Eigen::MatrixXd result = Eigen::MatrixXd::Identity(3, 3);
		result(0, 2) = p.x();
		result(1, 2) = p.y();
		return result;
	}

	PolyBase::Epipole PolyBase::ComputeLeftEpipole(const Eigen::MatrixXd& F) const
	{
		Eigen::MatrixXd W, U, VT;
		Eigen::JacobiSVD<Eigen::MatrixXf> svd(F.transpose(), Eigen::ComputeFullU | Eigen::ComputeFullV);
		U = svd.matrixU();
		VT = svd.matrixV().transpose();
		W = svd.singularValues();
		return Epipole(VT.row(2));
	}

	PolyBase::Epipole PolyBase::ComputeRightEpipole(const Eigen::MatrixXd& F) const
	{
		Eigen::MatrixXd W, U, VT;
		Eigen::JacobiSVD<Eigen::MatrixXf> svd(F, Eigen::ComputeFullU | Eigen::ComputeFullV);
		U = svd.matrixU();
		VT = svd.matrixV().transpose();
		W = svd.singularValues();
		return Epipole(VT.row(2));
	}

	Eigen::MatrixXd PolyBase::FormRotationMatrix(const Epipole& e) const
	{
		Eigen::MatrixXd result = Eigen::MatrixXd::Identity(3, 3);
		result(0, 0) = e.x();
		result(0, 1) = e.y();
		result(1, 0) = -e.y();
		result(1, 1) = e.x();
		return result;
	}

	int PolyBase::FindPolynomialOrder(const std::vector<double>& coeffs) const
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
	{ //not sure at all !!!!!!!!!!!!!!! maybe all what i did here is a mistake
		std::vector<double> coeffs = PreparePolyCoeffs(params);
		if (coeffs.size() <= 1)
		{
			return { 0 };
		}
		std::vector<Eigen::Vector2d> roots;
		cv::solvePoly(coeffs, roots); //roots is the output of solvePoly with parametr coeffs

		Roots result(roots.size());
		for (size_t i = 0u; i < roots.size(); ++i)
		{
			result[i] = std::complex<double>(roots[i][0], roots[i][1]);
		}
		return result;
	}

	double PolyBase::EvaluateRoots(const PolyBase::Roots& roots, const PolyParams& params) const
	{
		std::vector<double> costs = EvaluateRootsCosts(roots, params);
		return roots[std::min_element(costs.begin(), costs.end()) - costs.begin()].real();
	}

	std::pair<PolyBase::Line, PolyBase::Line> PolyBase::ConstructLines(double t, const PolyBase::PolyParams& params) const
	{
		double a, b, c, d, e, f;
		std::tie(a, b, c, d, e, f) = params;
		Line l0(t * e, 1, -t);
		Line l1(-f * (c * t + d), a * t + b, c * t + d);
		return std::make_pair(l0, l1);
	}

	Eigen::Vector3d PolyBase::FindPointOnLineClosestToOrigin(const PolyBase::Line& l) const
	{
		return Eigen::Vector3d(-l[0] * l[2], -l[1] * l[2], l[0] * l[0] + l[1] * l[1]);
	}

	Eigen::Vector3d PolyBase::TransferPointToOriginalCoordinates(const Eigen::Vector3d& p, const Eigen::MatrixXd& R, const Eigen::MatrixXd& T) const
	{
		Eigen::MatrixXd x = Eigen::MatrixXd::Identity(3, 1);
		x(0) = p.x();
		x(1) = p.y();
		x(2) = p.z();
		x = T * R.transpose() * x;
		return Eigen::Vector3d(x(0), x(1), x(2));
	}

	Eigen::MatrixXd PolyBase::CameraProjectionMatrixFromFundamentalMatrix(const PolyBase::Fundamental& F) const
	{ //not sure at all !!!!!!!!!!!!!!! maybe all what i did here is a mistake
		Eigen::MatrixXd e2;
		Eigen::JacobiSVD<Eigen::MatrixXd> svd(F.transpose(), Eigen::ComputeThinU | Eigen::ComputeThinV);
		Eigen::MatrixXd Z = svd.solve(e2);
		Eigen::Matrix3d e2x;
		e2x << 0, -Z(2), Z(1),
			Z(2), 0, -Z(0),
			-Z(1), Z(0), 0;
		Eigen::MatrixXd result = Eigen::MatrixXd::Identity(3, 4);
		Eigen::MatrixXd mul = e2x * F;

		//from here this is wrong for sure
		result(cv(0, 0, 3, 3)) = mul;
		result(0, 3) = Z(0);
		result(1, 3) = Z(1);
		result(2, 3) = Z(2);
		return result;
	}
}
