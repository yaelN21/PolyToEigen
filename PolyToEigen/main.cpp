#include <iostream>
#include <Eigen/Dense>
//#include <gtest/gtest.h>
#include "Poly.h"
#include <tuple>
//#include "PolyAbs.h"
#include "LinearLS.h"
#include <vector>

 
using Eigen::MatrixXd;
std::pair<Eigen::MatrixXd, Eigen::MatrixXd> SetupGeneralCameraConfiguration();


void EvaluateResult(const Eigen::Vector3d& result, const Eigen::Vector3d& expected_result, double max_percentage_error = 0.001);
int main()
{
    Eigen::MatrixXd P0, P1;

	auto result_setup = SetupGeneralCameraConfiguration();
	P0 =result_setup.first;
	P1 = result_setup.second;
	//std::tie(P0, P1) = SetupGeneralCameraConfiguration();
	//std::pair<Eigen::MatrixXd, Eigen::MatrixXd> SetupPair;
	//SetupPair = SetupGeneralCameraConfiguration();
	//Eigen::MatrixXd P0 = SetupPair.first;
	//Eigen::MatrixXd P1  = SetupPair.second;
	//(P0, P1) = SetupGeneralCameraConfiguration();
	Triangulation::Poly p(P0, P1);
    //example of 2 points
    Eigen::Vector3d point1(146, 642.288);
    Eigen::Vector3d point2(1137.31, 385.201);
	Eigen::Vector3d result = p.triangulate(point1,point2);
	Eigen::Vector3d expected_result(0.0, 100.0, 10000.0);
	EvaluateResult(result, expected_result);
    
}


void EvaluateResult(const Eigen::Vector3d& result, const Eigen::Vector3d& expected_result, double max_percentage_error = 0.001)
{   double expect_norm=  expected_result.norm();
	double tolerance = expect_norm * max_percentage_error;
    Eigen::Vector3d dist = result - expected_result;
	double distance = dist.norm();
    std::cout  << distance << tolerance <<std::endl;
	//EXPECT_NEAR(distance, 0.0, tolerance);
}