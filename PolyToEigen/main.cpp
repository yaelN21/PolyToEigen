#include <iostream>
#include <Eigen/Dense>
//#include <gtest/gtest.h>
#include "Poly.h"
#include <tuple>
//#include "PolyAbs.h"
#include "LinearLS.h"
#include <vector>
#include <Eigen/Core>

std::pair<Eigen::MatrixXf, Eigen::MatrixXf> SetupGeneralCameraConfiguration()
{
  float scale =   0.000001;
  
    Eigen::MatrixXf P0(3, 4);
    Eigen::MatrixXf P1(3, 4);
    Eigen::MatrixXf K(3, 3);

    P0 << 0.999701, 0.0174497, -0.017145, -500*scale,
          -0.0171452, 0.999695, 0.0177517, -100*scale,
          0.0174497, -0.0174524, 0.999695, -100*scale;

    P1 << 0.99969, -0.0174497, 0.0177543, 500*scale,
          0.0177543, 0.999695, -0.0171425, -100*scale,
          -0.0174497, 0.0174524, 0.999695, -100*scale;
     
    K << 7291.67, 0, 639.5,
         0, 7291.67, 511.5,
         0, 0, 1;

    P0 = K * P0;
    P1 = K * P1;
   
    
    return std::make_pair(P0, P1);
}

using Eigen::MatrixXf;
std::pair<Eigen::MatrixXf, Eigen::MatrixXf> SetupGeneralCameraConfiguration();


void EvaluateResult(const Eigen::Vector3f& result, const Eigen::Vector3f& expected_result);
int main(int argc, char *argv[])
{
  
	auto result_setup = SetupGeneralCameraConfiguration();
	const Eigen::MatrixXf P0 =result_setup.first;
	const Eigen::MatrixXf P1 = result_setup.second;
	//std::tie(P0, P1) = SetupGeneralCameraConfiguration();
	//std::pair<Eigen::MatrixXf, Eigen::MatrixXf> SetupPair;
	//SetupPair = SetupGeneralCameraConfiguration();
	//Eigen::MatrixXf P0 = SetupPair.first;
	//Eigen::MatrixXf P1  = SetupPair.second;
	//(P0, P1) = SetupGeneralCameraConfiguration();
	Triangulation::Poly p(P0, P1);
    //example of 2 points
   float scale =  0.000001;
    Eigen::Vector2f point1(1004.08, 511.5);
    Eigen::Vector2f point2(274.917, 511.5);
   point1 *=scale;
   point2 *=scale;
	Eigen::Vector3f result = p.triangulate(point1,point2);
	Eigen::Vector3f expected_result(500.0, 0.0, 10000.0);
  expected_result*=scale;
    std::cout << "yay" << std::endl;
	//EvaluateResult(result, expected_result);

}


void EvaluateResult(const Eigen::Vector3f& result, const Eigen::Vector3f& expected_result)
{   float expect_norm=  expected_result.norm();
	float tolerance = expect_norm * 0.001;
    Eigen::Vector3f dist = result - expected_result;
	float distance = dist.norm();
    std::cout  << distance << tolerance <<std::endl;
	//EXPECT_NEAR(distance, 0.0, tolerance);
}