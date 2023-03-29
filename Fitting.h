#pragma once

#include "VxApi.h"
#include "Eigen/Core"
#include "Eigen/Dense"

int FittingInit(int format, void* data);
int FittingExit(void);
int Fitting(int idData, void* echo);
int FittingEO(int idData, void* echo);
int FittingTerm(int idData);
//绘制多段线echo函数
int DrawPolyLineEO(int idData);
//绘制多段线
int DrawPolyLine(svxPoint* ptrPntList, int nPntCount, evxColor lineColor = VX_COLOR_BLACK);

//******************************  成员变量  ******************************** //
const int m_nPolyExp = 4;//
Eigen::VectorXd m_vecdCoe;
Eigen::VectorXd m_vecdCoeX;
Eigen::VectorXd m_vecdCoeY;
const int m_iSampleCount = 100;
svxPoint* m_ptrFittingPntList = nullptr;
const double m_dStandardVariance = 10;
const double m_dnLambda = 0.5;

double m_dStartX = 0;
double m_dEndX = 0;
Eigen::MatrixXd m_matdCoe;


//******************************
//       插值型拟合方法     //
//************************

//多项式基函数插值拟合
//void Interpolation_PolynomialBaseFunction(svxPoint* ptrPntList, int nPntCount);
void Interpolation_PolynomialBaseFunction(svxPoint* ptrPntList, int nPntCount, svxPoint* ptrResultPntList = nullptr);
//高斯基函数插值拟合
void Interpolation_GaussBaseFunction(svxPoint* ptrPntList, int nPntCount);

//Lagrange插值多项式

//Newton插值多项式

//计算多项式拟合表达式结果
//void CalculatePolynominalExpress(double dStart, double dEnd);
void CalculatePolynominalExpress(double dStart, double dEnd, svxPoint* resultPntList = m_ptrFittingPntList);

//计算Gauss拟合曲线上采样点的值
void CalculateGaussExpress(double dStart, double dEnd, svxPoint* ptrPntList);

//********************************
//       逼近型拟合方法         //
//******************************

//最小二乘逼近拟合
bool Approximation_LeastSquare(svxPoint* ptrPntList, int nPntCount);

//岭回归逼近拟合
bool Approximation_RidgeRegression(svxPoint* ptrPntList, int nPntCount);

//******************************  二维平面曲线拟合  ******************************** //


//均匀参数化
bool UniformParameterization(svxPoint* ptrPntList, int nPntCount);

//计算多项式拟合表达式结果
void Calculate2DPolynominalExpress();


