#include <iostream>
#include <Eigen/Dense>
//#include <gtest/gtest.h>
#include "../include/Poly.h"
#include <tuple>
//#include "PolyAbs.h"
#include "../include/LinearLS.h"
#include <vector>
#include <Eigen/Core>
#include <opencv2/opencv.hpp>


std::pair<Eigen::Matrix<float,3,4>, Eigen::Matrix<float,3,4>> SetupGeneralCameraConfiguration()
{
    Eigen::Matrix<float,3,4> P0;
    Eigen::Matrix<float,3,4> P1;
    Eigen::Matrix3f K;

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

std::pair<Eigen::Matrix<float,3,4>, Eigen::Matrix<float,3,4>> SetupSecondCameraRotatedLeftConfiguration()
{
   Eigen::Matrix<float,3,4> P0;
    Eigen::Matrix<float,3,4> P1;

    P0 << 0.999701, 0.0174497, -0.017145, 0,
          -0.0171452, 0.999695, 0.0177517, 0,
          0.0174497, -0.0174524, 0.999695, 0;

    P1 << 1, 0, 0, -1000,
          0, 1, 0, 0,
          0, 0, 1, 0;

    Eigen::Matrix3f K;
    K << 7291.67, 0, 639.5,
         0, 7291.67, 511.5,
         0, 0, 1;

    P0 = K * P0;
    P1 = K * P1;

    return std::make_pair(P0, P1);
}


std::pair<Eigen::Matrix<float,3,4>, Eigen::Matrix<float,3,4>> SetupSecondCameraRotatedRightConfiguration()
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
std::pair<Eigen::Matrix<float,3,4>, Eigen::Matrix<float,3,4>> SetupGeneralCameraConfiguration();
std::pair<Eigen::Matrix<float,3,4>, Eigen::Matrix<float,3,4>> SetupSecondCameraRotatedRightConfiguration();
std::pair<Eigen::Matrix<float,3,4>, Eigen::Matrix<float,3,4>> SetupSecondCameraRotatedLeftConfiguration();
void EvaluateResult(const Eigen::Vector3f& result, const Eigen::Vector3f& expected_result);
void TestPoly(std::pair<Eigen::Matrix<float,3,4>, Eigen::Matrix<float,3,4>> result_setup,Eigen::Vector3f point1,Eigen::Vector3f point2, Eigen::Vector3f  expected_result);
bool TriangulateTest(Eigen::Vector3f &x_c1, Eigen::Vector3f &x_c2,
Eigen::Matrix<float,3,4> &Tc1w ,Eigen::Matrix<float,3,4> &Tc2w ,   Eigen::Matrix3f mK_1, Eigen::Matrix3f mK_2,
Eigen::Matrix3f  Rcw1 ,Eigen::Matrix3f Rcw2,Eigen::Vector3f tcw1,
Eigen::Vector3f tcw2 ,Eigen::Vector3f &x3D);

void TestPoly(std::pair<Eigen::Matrix<float,3,4>, Eigen::Matrix<float,3,4>> result_setup,Eigen::Vector3f point1,Eigen::Vector3f point2, Eigen::Vector3f  expected_result)
{

  Eigen::Matrix<float,3,4> P0 =std::get<0>(result_setup);
        Eigen::Matrix<float,3,4> P1 =std::get<1>(result_setup);
        cv::Mat K0_cv,R0_cv, T0_cv, K1_cv,R1_cv, T1_cv;
        cv::Mat cvProjectionMatrix0(3, 4, CV_32F);
        cv::Mat cvProjectionMatrix1(3, 4, CV_32F);
        Eigen::Map<Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>>(
        cvProjectionMatrix0.ptr<float>(), cvProjectionMatrix0.rows, cvProjectionMatrix0.cols) = P0;
        Eigen::Map<Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>>(
        cvProjectionMatrix1.ptr<float>(), cvProjectionMatrix1.rows, cvProjectionMatrix1.cols) = P1;
        cv::decomposeProjectionMatrix(cvProjectionMatrix0,K0_cv,R0_cv,T0_cv);
        cv::decomposeProjectionMatrix(cvProjectionMatrix1,K1_cv,R1_cv,T1_cv);
        Eigen::MatrixXf R0_eigen = Eigen::Map<Eigen::Matrix<float, 3, 3, Eigen::RowMajor>>(R0_cv.ptr<float>());
        Eigen::Vector4f T0_eigen = Eigen::Map<Eigen::Vector4f>(T0_cv.ptr<float>());
         Eigen::Vector4f T1_eigen = Eigen::Map<Eigen::Vector4f>(T1_cv.ptr<float>());
        Eigen::MatrixXf K0_eigen = Eigen::Map<Eigen::Matrix<float, 3, 3, Eigen::RowMajor>>(K0_cv.ptr<float>());
        Eigen::MatrixXf K1_eigen = Eigen::Map<Eigen::Matrix<float, 3, 3, Eigen::RowMajor>>(K1_cv.ptr<float>());
        Eigen::MatrixXf R1_eigen = Eigen::Map<Eigen::Matrix<float, 3, 3, Eigen::RowMajor>>(R1_cv.ptr<float>());
        // Divide each coordinate by the last element
        Eigen::Vector3f t0 = T0_eigen.head<3>().array() / T0_eigen(3);
        Eigen::Vector3f t1 = T1_eigen.head<3>().array() / T1_eigen(3);
         Eigen::Vector3f x3D;
        bool result = TriangulateTest(point1,point2, P0,P1,K0_eigen,K1_eigen,R0_eigen,R1_eigen,t0,t1,x3D);
        EvaluateResult(x3D,expected_result);
}
int main(int argc, char *argv[])
{

     //General
    std::cout << "Testing SetupGeneralCameraConfiguration:" << std::endl;
	auto result_setup = SetupGeneralCameraConfiguration();
    Eigen::Vector3f point1(146, 642.288,1);
    Eigen::Vector3f point2(1137.31, 385.201,1);
    Eigen::Vector3f expected_result(0.0, 100.0, 10000.0);
    TestPoly(result_setup,point1,point2,expected_result);
  //Left
   std::cout << "Testing SetupSecondCameraRotatedLeftConfiguration:" << std::endl;
    auto result_setup_L = SetupSecondCameraRotatedLeftConfiguration();
    Eigen::Vector3f point1_L(878.821, 634.619,1);
    Eigen::Vector3f point2_L(274.917, 511.5,1);
    Eigen::Vector3f expected_result_L(500.0, 0.0, 10000.0);
    TestPoly(result_setup_L,point1_L,point2_L,expected_result_L);
//Right
   std::cout << "Testing SetupSecondCameraRotatedRightConfiguration:" << std::endl;
    auto result_setup_R = SetupSecondCameraRotatedRightConfiguration();
    Eigen::Vector3f point1_R(1004.08, 511.5,1);
    Eigen::Vector3f point2_R(150.068, 634.618,1);
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

bool TriangulateTest(Eigen::Vector3f &x_c1, Eigen::Vector3f &x_c2,
Eigen::Matrix<float,3,4> &Tc1w ,Eigen::Matrix<float,3,4> &Tc2w ,   Eigen::Matrix3f mK_1, Eigen::Matrix3f mK_2,
Eigen::Matrix3f  Rcw1 ,Eigen::Matrix3f Rcw2,Eigen::Vector3f tcw1,
Eigen::Vector3f tcw2 ,Eigen::Vector3f &x3D)
{
  Eigen::VectorXf point1(2);
  point1 << x_c1[0], x_c1[1];
  Eigen::VectorXf point2(2);
 point2 << x_c2[0], x_c2[1];
 
Triangulation::Poly p(Tc1w,Tc2w,mK_1,mK_2,Rcw1,Rcw2,tcw1,tcw2);

    Eigen::Vector3f result = p.triangulate(point1,point2);
     if(result.isZero())
    {
        std::cout << "error in triangulate" << std::endl;
            return false;
    }
    x3D = result;
    return  true;
    
}