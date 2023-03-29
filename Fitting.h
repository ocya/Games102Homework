#pragma once

#include "VxApi.h"
#include "Eigen/Core"
#include "Eigen/Dense"

int FittingInit(int format, void* data);
int FittingExit(void);
int Fitting(int idData, void* echo);
int FittingEO(int idData, void* echo);
int FittingTerm(int idData);
//���ƶ����echo����
int DrawPolyLineEO(int idData);
//���ƶ����
int DrawPolyLine(svxPoint* ptrPntList, int nPntCount, evxColor lineColor = VX_COLOR_BLACK);

//******************************  ��Ա����  ******************************** //
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
//       ��ֵ����Ϸ���     //
//************************

//����ʽ��������ֵ���
//void Interpolation_PolynomialBaseFunction(svxPoint* ptrPntList, int nPntCount);
void Interpolation_PolynomialBaseFunction(svxPoint* ptrPntList, int nPntCount, svxPoint* ptrResultPntList = nullptr);
//��˹��������ֵ���
void Interpolation_GaussBaseFunction(svxPoint* ptrPntList, int nPntCount);

//Lagrange��ֵ����ʽ

//Newton��ֵ����ʽ

//�������ʽ��ϱ��ʽ���
//void CalculatePolynominalExpress(double dStart, double dEnd);
void CalculatePolynominalExpress(double dStart, double dEnd, svxPoint* resultPntList = m_ptrFittingPntList);

//����Gauss��������ϲ������ֵ
void CalculateGaussExpress(double dStart, double dEnd, svxPoint* ptrPntList);

//********************************
//       �ƽ�����Ϸ���         //
//******************************

//��С���˱ƽ����
bool Approximation_LeastSquare(svxPoint* ptrPntList, int nPntCount);

//��ع�ƽ����
bool Approximation_RidgeRegression(svxPoint* ptrPntList, int nPntCount);

//******************************  ��άƽ���������  ******************************** //


//���Ȳ�����
bool UniformParameterization(svxPoint* ptrPntList, int nPntCount);

//�������ʽ��ϱ��ʽ���
void Calculate2DPolynominalExpress();


