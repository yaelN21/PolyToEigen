#include <iostream>
#include <Eigen/Dense>
//#include <gtest/gtest.h>
#include "Poly.h"
#include <tuple>
//#include "PolyAbs.h"
#include "LinearLS.h"
#include <vector>
#include <Eigen/Core>

std::pair<Eigen::MatrixXd, Eigen::MatrixXd> SetupGeneralCameraConfiguration()
{
    Eigen::MatrixXd P0(3, 4);
    Eigen::MatrixXd P1(3, 4);
    Eigen::MatrixXd K(3, 3);

    P0 << 0.999701, 0.0174497, -0.017145, -500,
          -0.0171452, 0.999695, 0.0177517, -100,
          0.0174497, -0.0174524, 0.999695, -100;

    P1 << 0.99969, -0.0174497, 0.0177543, 500,
          0.0177543, 0.999695, -0.0171425, -100,
          -0.0174497, 0.0174524, 0.999695, -100;

    K << 7291.67, 0, 639.5,
         0, 7291.67, 511.5,
         0, 0, 1;

    P0 = K * P0;
    P1 = K * P1;
    
    return std::make_pair(P0, P1);
}

using Eigen::MatrixXd;
std::pair<Eigen::MatrixXd, Eigen::MatrixXd> SetupGeneralCameraConfiguration();


void EvaluateResult(const Eigen::Vector3d& result, const Eigen::Vector3d& expected_result);
int main(int argc, char *argv[])
{

	auto result_setup = SetupGeneralCameraConfiguration();
	const Eigen::MatrixXd P0 =result_setup.first;
	const Eigen::MatrixXd P1 = result_setup.second;
	//std::tie(P0, P1) = SetupGeneralCameraConfiguration();
	//std::pair<Eigen::MatrixXd, Eigen::MatrixXd> SetupPair;
	//SetupPair = SetupGeneralCameraConfiguration();
	//Eigen::MatrixXd P0 = SetupPair.first;
	//Eigen::MatrixXd P1  = SetupPair.second;
	//(P0, P1) = SetupGeneralCameraConfiguration();
	Triangulation::Poly p(P0, P1);
    //example of 2 points
    Eigen::Vector2d point1(1004.08, 511.5);
    Eigen::Vector2d point2(274.917, 511.5);

	Eigen::Vector3d result = p.triangulate(point1,point2);
	Eigen::Vector3d expected_result(500.0, 0.0, 10000.0);
	EvaluateResult(result, expected_result);
    std::cout << "yay" << std::endl;
    return 0;
}


void EvaluateResult(const Eigen::Vector3d& result, const Eigen::Vector3d& expected_result)
{   double expect_norm=  expected_result.norm();
	double tolerance = expect_norm * 0.001;
    Eigen::Vector3d dist = result - expected_result;
	double distance = dist.norm();
    std::cout  << distance << tolerance <<std::endl;
	//EXPECT_NEAR(distance, 0.0, tolerance);
}