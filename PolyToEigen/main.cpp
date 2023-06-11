#include <iostream>
#include <Eigen/Dense>
//#include <gtest/gtest.h>
#include "Poly.h"
#include <tuple>
//#include "PolyAbs.h"
#include "LinearLS.h"
#include <vector>
#include <Eigen/Core>


std::tuple<Eigen::MatrixXf, Eigen::MatrixXf> SetupGeneralCameraConfiguration()
{
  //float scale =   0.000001;
    Eigen::MatrixXf P0(3, 4);
    Eigen::MatrixXf P1(3, 4);
    Eigen::MatrixXf K(3, 3);

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


std::pair<Eigen::MatrixXf, Eigen::MatrixXf> SetupSecondCameraRotatedLeftConfiguration()
{
    Eigen::MatrixXf P0(3, 4);
    Eigen::MatrixXf P1(3, 4);

    P0 << 0.999701, 0.0174497, -0.017145, 0,
          -0.0171452, 0.999695, 0.0177517, 0,
          0.0174497, -0.0174524, 0.999695, 0;

    P1 << 1, 0, 0, -1000,
          0, 1, 0, 0,
          0, 0, 1, 0;

    Eigen::MatrixXf K(3, 3);
    K << 7291.67, 0, 639.5,
         0, 7291.67, 511.5,
         0, 0, 1;

    P0 = K * P0;
    P1 = K * P1;

    return std::make_pair(P0, P1);
}


std::pair<Eigen::MatrixXf, Eigen::MatrixXf> SetupSecondCameraRotatedRightConfiguration()
{
    Eigen::Matrix<float, 3, 4> P0;
    P0 << 1, 0, 0, 0,
          0, 1, 0, 0,
          0, 0, 1, 0;

    Eigen::Matrix<float, 3, 4> P1;
    P1 << 0.999701f, 0.0174497f, -0.017145f, -1000.0f,
          -0.0171452f, 0.999695f, 0.0177517f, 0,
          0.0174497f, -0.0174524f, 0.999695f, 0;

    Eigen::Matrix<float, 3, 3> K = Eigen::Matrix<float, 3, 3>::Identity();
    K(0, 0) = 7291.67f;
    K(1, 1) = 7291.67f;
    K(0, 2) = 639.5f;
    K(1, 2) = 511.5f;

    P0 = K * P0;
    P1 = K * P1;

    return std::make_pair(P0, P1);
}

using Eigen::MatrixXf;
std::tuple<Eigen::MatrixXf, Eigen::MatrixXf> SetupGeneralCameraConfiguration();
std::pair<Eigen::MatrixXf, Eigen::MatrixXf> SetupSecondCameraRotatedRightConfiguration();
std::pair<Eigen::MatrixXf, Eigen::MatrixXf> SetupSecondCameraRotatedLeftConfiguration();
void EvaluateResult(const Eigen::Vector3f& result, const Eigen::Vector3f& expected_result);
void TestPoly(std::pair<Eigen::MatrixXf, Eigen::MatrixXf> result_setup,Eigen::Vector2f point1,Eigen::Vector2f point2, Eigen::Vector3f  expected_result);


void TestPoly(std::pair<Eigen::MatrixXf, Eigen::MatrixXf> result_setup,Eigen::Vector2f point1,Eigen::Vector2f point2, Eigen::Vector3f  expected_result)
{
    const Eigen::MatrixXf P0 = std::get<0>(result_setup);
	const Eigen::MatrixXf P1 = std::get<1>(result_setup);
    Triangulation::Poly p(P0, P1);
    Eigen::Vector3f result = p.triangulate(point1,point2);
	EvaluateResult(result,expected_result);

}
int main(int argc, char *argv[])
{
    //General
    std::cout << "Testing SetupGeneralCameraConfiguration:" << std::endl;
	auto result_setup = SetupGeneralCameraConfiguration();
    Eigen::Vector2f point1(146, 642.288);
    Eigen::Vector2f point2(1137.31, 385.201);
    Eigen::Vector3f expected_result(0.0, 100.0, 10000.0);
    TestPoly(result_setup,point1,point2,expected_result);
  //Left
   std::cout << "Testing SetupSecondCameraRotatedLeftConfiguration:" << std::endl;
    auto result_setup_L = SetupSecondCameraRotatedLeftConfiguration();
    Eigen::Vector2f point1_L(878.821, 634.619);
    Eigen::Vector2f point2_L(274.917, 511.5);
    Eigen::Vector3f expected_result_L(500.0, 0.0, 10000.0);
    TestPoly(result_setup_L,point1_L,point2_L,expected_result_L);
//Right
   std::cout << "Testing SetupSecondCameraRotatedRightConfiguration:" << std::endl;
    auto result_setup_R = SetupSecondCameraRotatedRightConfiguration();
    Eigen::Vector2f point1_R(1004.08, 511.5);
    Eigen::Vector2f point2_R(150.068, 634.618);
    Eigen::Vector3f expected_result_R(500.0, 0.0, 10000.0);
    TestPoly(result_setup_R,point1_R,point2_R,expected_result_R);
    std::cout << "yay" << std::endl;
    

}

void assertNear(const Eigen::Vector3f& distance, float tolerance)
{
    float norm = distance.norm();
    if (norm > tolerance)
    {
        std::cerr << "Distance is not within tolerance. Norm: " << norm << ", Tolerance: " << tolerance << std::endl;
        // You can choose to throw an exception or exit the program here
        // throw std::runtime_error("Distance is not within tolerance");
        // std::exit(EXIT_FAILURE);
    }
    else
    {
        std::cout << "Distance is within tolerance. Norm: " << norm << ", Tolerance: " << tolerance << std::endl;
    }
}

void EvaluateResult(const Eigen::Vector3f& result, const Eigen::Vector3f& expected_result)
{   float max_percentage_error = 0.001;
    float tolerance=  expected_result.norm()*max_percentage_error;
    Eigen::Vector3f dist = result - expected_result;
	 assertNear(dist, tolerance);
     
}

