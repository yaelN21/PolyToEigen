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
float scale =  1;
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
   
   
    return std::make_tuple(P0, P1);
}

using Eigen::MatrixXf;
std::tuple<Eigen::MatrixXf, Eigen::MatrixXf> SetupGeneralCameraConfiguration();


void EvaluateResult(const Eigen::Vector3f& result, const Eigen::Vector3f& expected_result);
int main(int argc, char *argv[])
{
  
	auto result_setup = SetupGeneralCameraConfiguration();
 
	const Eigen::MatrixXf P0 = std::get<0>(result_setup);
	const Eigen::MatrixXf P1 = std::get<1>(result_setup);
	Triangulation::Poly p(P0, P1);
    //example of 2 pointsexpected_result
    //float scale =  0.000001;
    float scale =  1;
    Eigen::Vector2f point1(146, 642.288);
    Eigen::Vector2f point2(1137.31, 385.201);
    point1 *=scale;
    point2 *=scale;
	Eigen::Vector3f result = p.triangulate(point1,point2);
	Eigen::Vector3f expected_result(0.0, 100.0, 10000.0);
    expected_result*=scale;
    std::cout << "yay" << std::endl;
	EvaluateResult(result,expected_result);

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
	float distance = dist.norm();
	 assertNear(dist, tolerance);
     
}

