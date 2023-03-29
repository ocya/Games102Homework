#include "Fitting.h"
#include <cmath>
#define	FITTING_UI "Fitting"

int FittingInit(int format, void* data)
{
	cvxCmdFunc("Fitting", (void*)Fitting, VX_CODE_GENERAL);
	cvxCmdFunc("FittingEO", (void*)FittingEO, VX_CODE_GENERAL);
	cvxCmdFunc("FittingTerm", (void*)FittingTerm, VX_CODE_GENERAL);
	return 0;
}


int FittingExit(void)
{
	cvxCmdFuncUnload("Fitting");
	return 0;
}

//int cmp(const void* pntFirst, const void* pntSecond) 
//{
//	return (*((svxPoint*)pntSecond)).x - (*((svxPoint*)pntFirst)).x;
//}

int Fitting(int idData, void *echo)
{
	//使模型空间的视角变为xy平面
	cvxCmdBuffer("^viewTop", 0);
	//acquire user input points list
	int nRet = 0;
	int nPntCount = 0;
	svxPoint* ptrPntList = nullptr;
	nRet = cvxDataGetPnts(idData, 2, &nPntCount, &ptrPntList);
	if (1 == nRet)
	{
		cvxMemFree((void**)&ptrPntList);
		return nRet;
	}

	//qsort(ptrPntList, nPntCount, sizeof(svxPoint), cmp);
	//对用户输入的点集合进行排序	
	//int nListLength = sizeof(ptrPntList) / sizeof(svxPoint);
	
	//std::sort(ptrPntList, ptrPntList + nPntCount, [](svxPoint pntFirst, svxPoint pntSecond) {
	//	return pntFirst.x < pntSecond.x;
	//	});

	//m_dStartX = ptrPntList[0].x;
	//m_dEndX = ptrPntList[nPntCount - 1].x;

	//绘制用户选择点
	nRet = cvxPartPnts(nPntCount, ptrPntList);
	if (1 == nRet)
	{
		cvxMemFree((void**)&ptrPntList);
		return nRet;
	}
	//绘制多项式基函数插值拟合曲线
	//nRet = DrawPolyLine(ptrPntList, nPntCount);
	/*	
	Interpolation_PolynomialBaseFunction(ptrPntList, nPntCount);
	if (nullptr != m_ptrFittingPntList)
	{
		nRet = DrawPolyLine(m_ptrFittingPntList, m_iSampleCount, VX_COLOR_BROWN);
	}
	*/


	//绘制Gauss基函数插值拟合曲线
	/*
	Interpolation_GaussBaseFunction(ptrPntList, nPntCount);
	if (nullptr != m_ptrFittingPntList)
	{
		nRet = DrawPolyLine(m_ptrFittingPntList, m_iSampleCount, VX_COLOR_MED_BLUE);
	}
	*/


	//绘制最小二乘逼近拟合曲线
	/*
	Approximation_LeastSquare(ptrPntList, nPntCount);
	if (nullptr != m_ptrFittingPntList)
	{
		nRet = DrawPolyLine(m_ptrFittingPntList, m_iSampleCount, VX_COLOR_YELLOW);
	}

	//绘制岭回归逼近拟合曲线
	Approximation_RidgeRegression(ptrPntList, nPntCount);
	if (nullptr != m_ptrFittingPntList)
	{
		nRet = DrawPolyLine(m_ptrFittingPntList, m_iSampleCount, VX_COLOR_RED);
	}
	*/
	//绘制均匀参数化二维曲线拟合
	bool bRet = UniformParameterization(ptrPntList, nPntCount);
	if (true == bRet && nullptr != m_ptrFittingPntList)
	{
		DrawPolyLine(m_ptrFittingPntList, m_iSampleCount, VX_COLOR_GREEN);
	}


	//清理内存
	cvxMsgDisp("Fitting successed");
	cvxMemFree((void**)&ptrPntList); //清零用户输入点的内存
	//cvxMemFree((void**)&m_ptrFittingPntList); //清零拟合方程上采样点的内存


	return nRet;
}


int FittingEO(int idData, void* echo)
{
	cvxEchoStart();
	DrawPolyLineEO(idData);
	cvxEchoEnd();
	return 0;
}


int FittingTerm(int idData)
{
	return 0;
}


//绘制多段线echo函数
int DrawPolyLineEO(int idData)
{
	
	int nPolylineIdData = 0;	
	svxData tempData;
	cvxDataInit("CdLnPoly", &nPolylineIdData);
	cvxDataZero(&tempData);
	tempData.isPoint = 1;
	tempData.PntType = VX_PNT3_ABS;
	
	int nRet = 0;
	int nPntCount = 0;
	svxPoint* ptrPntList = nullptr;

	//update acquired points list
	cvxDataGetPnts(idData, 2, &nPntCount, &ptrPntList);

	for (int i = 0; i < nPntCount; i++)
	{
		tempData.Pnt = ptrPntList[i];
		cvxDataSet(nPolylineIdData, 1, &tempData);
	}
	cvxCmdExec(nPolylineIdData);


	return 0;
}

//绘制多段线
int DrawPolyLine(svxPoint* ptrPntList, int nPntCount, evxColor lineColor/* = VX_COLOR_BLACK*/)
{
	cvxDispColorSet(VX_DISP_LINE, lineColor);
	//初始化多段线绘制的ui
	int nRet = 0;
	int nPolylineIdData = 0;
	svxData tempData;
	nRet = cvxDataInit("CdLnPoly", &nPolylineIdData);
	if (1 == nRet)
	{
		return nRet;
	}

	//数据赋值
	for (int i = 0; i < nPntCount; i++)
	{
		tempData.Pnt = ptrPntList[i];
		//对于PointList，持续调用cvxDataSet，即为向List中添加Point
		nRet = cvxDataSet(nPolylineIdData, 1, &tempData);
		if (1 == nRet)
		{
			return nRet;
		}
	}

	//执行指令
	nRet = cvxCmdExec(nPolylineIdData);
	//cvxDataFree();

	return nRet;
}

//多项式拟合
void Interpolation_PolynomialBaseFunction(svxPoint* ptrPntList, int nPntCount, svxPoint* ptrResultPntList/* = nullptr*/)
{
	//根据点构造x矩阵和y向量
	Eigen::MatrixXd matdX(nPntCount, nPntCount);
	Eigen::VectorXd vecdY(nPntCount);
	for (int i = 0; i < nPntCount; i++)
	{
		vecdY(i) = ptrPntList[i].y;
		matdX(i, 0) = 1;
		for (int j = 1; j < nPntCount; j++)
		{
			matdX(i, j) = matdX(i, j - 1) * ptrPntList[i].x;
		}
	}
	//Eigen::VectorXd vecdCoe = matdX.inverse() * vecdY;
	m_vecdCoe.resize(nPntCount);// n行 1列
	m_vecdCoe = matdX.inverse() * vecdY; //inverse()转置矩阵
	
	
	//对拟合方程进行采样，获得绘图时所需要的数据点
	//double dStart = ptrPntList[0].x;
	//double dEnd = ptrPntList[nPntCount - 1].x;
	CalculatePolynominalExpress(m_dStartX, m_dEndX, ptrResultPntList);
}


void Interpolation_GaussBaseFunction(svxPoint* ptrPntList, int nPntCount)
{
	//根据点构造x矩阵和y向量
	Eigen::MatrixXd matdX(nPntCount+1, nPntCount + 1);// n+1行 n+1列
	Eigen::VectorXd vecdY(nPntCount+1);// 1行 n+1列

	//由于未知数个数为n+1，方程个数为n，增加一个约束点，即令拟合曲线通过X均值点，Y均值点。
	double dConstraintX = 0;
	double dConstraintY = 0;
	for (int i = 0; i < nPntCount; i++)
	{
		vecdY(i) = ptrPntList[i].y;
		dConstraintY += ptrPntList[i].y / nPntCount;//防止加数过大导致溢出
		dConstraintX += ptrPntList[i].x / nPntCount;
		matdX(i, 0) = 1;
		double dX = ptrPntList[i].x;
		double dXi = 0;
		for (int j = 1; j < nPntCount + 1; j++)
		{
			dXi = ptrPntList[j - 1].x;
			matdX(i, j) = exp((dX - dXi) * (dX - dXi) / (-2 * m_dStandardVariance * m_dStandardVariance));
		}
	}
	vecdY(nPntCount) = dConstraintY;
	matdX(nPntCount, 0) = 1;
	for (int i = 1; i < nPntCount + 1; i++)
	{
		matdX(nPntCount, i) = exp((dConstraintX - ptrPntList[i - 1].x) * (dConstraintX - ptrPntList[i - 1].x) / (-2 * m_dStandardVariance * m_dStandardVariance));
	}

	m_vecdCoe.resize(nPntCount+1);//系数矩阵1行 n+1列
	m_vecdCoe = matdX.inverse() * vecdY;

	//对拟合方程进行采样，获得绘图时所需要的数据点
	CalculateGaussExpress(m_dStartX, m_dEndX, ptrPntList);

}

//void CalculatePolynominalExpress(double dStart, double dEnd)
//{
//	
//	double x = 0;
//	double dStepLength = (dEnd - dStart) / (m_iSampleCount - 1);
//	int n = m_vecdCoe.size();//多项式的幂指数
//	m_ptrFittingPntList = new svxPoint[m_iSampleCount];
//	
//	for (int i = 0; i < m_iSampleCount; i++)
//	{
//		x = dStart + i * dStepLength;
//		m_ptrFittingPntList[i].x = x;
//		double dXExpTemp = 1;
//		double y = 0;
//		for (int j = 0; j < n; j++)
//		{
//			if(0 != j)
//			{
//				dXExpTemp *= x;//对x求幂指数
//			}		
//			y += m_vecdCoe[j] * dXExpTemp;
//		}
//		m_ptrFittingPntList[i].y = y;
//		m_ptrFittingPntList[i].z = 0;
//	}
//}

void CalculatePolynominalExpress(double dStart, double dEnd, svxPoint* resultPntList/* = m_ptrFittingPntList*/)
{
	double x = 0;
	double dStepLength = (dEnd - dStart) / (m_iSampleCount - 1);
	int n = m_vecdCoe.size();//多项式的幂指数
	/*resultPntList = new svxPoint[m_iSampleCount];*/

	for (int i = 0; i < m_iSampleCount; i++)
	{
		x = dStart + i * dStepLength;
		resultPntList[i].x = x;
		double dXExpTemp = 1;
		double y = 0;
		for (int j = 0; j < n; j++)
		{
			if (0 != j)
			{
				dXExpTemp *= x;//对x求幂指数
			}
			y += m_vecdCoe[j] * dXExpTemp;
		}
		resultPntList[i].y = y;
		resultPntList[i].z = 0;
	}
}

void CalculateGaussExpress(double dStart, double dEnd, svxPoint* ptrPntList)
{
	double dX = 0;
	double dStepLength = (dEnd - dStart) / (m_iSampleCount - 1);
	int n = m_vecdCoe.size();//多项式的幂指数
	m_ptrFittingPntList = new svxPoint[m_iSampleCount];

	for (int i = 0; i < m_iSampleCount; i++)
	{
		dX = dStart + i * dStepLength;
		m_ptrFittingPntList[i].x = dX;

		double dGXi = 1;
		double y = 0;
		double dXi = 0;
		for (int j = 0; j < n; j++)
		{
			if (0 != j)
			{
				dXi = ptrPntList[j - 1].x;
				dGXi = exp((dX - dXi) * (dX - dXi) / (-2 * m_dStandardVariance * m_dStandardVariance));
			}
			y += m_vecdCoe[j] * dGXi;
		}
		m_ptrFittingPntList[i].y = y;
		m_ptrFittingPntList[i].z = 0;
	}
}

bool Approximation_LeastSquare(svxPoint* ptrPntList, int nPntCount)
{
	if (m_nPolyExp >= nPntCount)
	{
		return false;
	}
	//根据点构造 x矩阵 和 y向量
	Eigen::MatrixXd matdX(nPntCount, m_nPolyExp);
	Eigen::VectorXd vecdY(nPntCount);
	for (int i = 0; i < nPntCount; i++)
	{
		vecdY(i) = ptrPntList[i].y;
		matdX(i, 0) = 1;
		for (int j = 1; j < m_nPolyExp; j++)
		{
			matdX(i, j) = matdX(i, j - 1) * ptrPntList[i].x;
		}
	}

	m_vecdCoe.resize(m_nPolyExp);
	m_vecdCoe = (matdX.transpose() * matdX).inverse() * (matdX.transpose() * vecdY); //最小二乘解法

	//对拟合方程进行采样，获得绘图时所需要的数据点
	CalculatePolynominalExpress(m_dStartX, m_dEndX);

	return true;
}

bool Approximation_RidgeRegression(svxPoint* ptrPntList, int nPntCount)
{
	if (m_nPolyExp >= nPntCount)
	{
		return false;
	}
	Eigen::MatrixXd matdX(nPntCount, m_nPolyExp);
	Eigen::VectorXd vecdY(nPntCount);
	for (int i = 0; i < nPntCount; i++)
	{
		matdX(i, 0) = 1;
		vecdY(i) = ptrPntList[i].y;
		double dTempX = ptrPntList[i].x;
		for (int j = 1; j < m_nPolyExp; j++)
		{
			matdX(i, j) = matdX(i, j - 1) * dTempX;
		}
	}
	m_vecdCoe.resize(m_nPolyExp);
	m_vecdCoe = (matdX.transpose() * matdX + m_dnLambda * Eigen::MatrixXd::Identity(m_nPolyExp, m_nPolyExp)).inverse() * matdX.transpose() * vecdY;

	CalculatePolynominalExpress(m_dStartX, m_dEndX);
	return true;

}

bool UniformParameterization(svxPoint* ptrPntList, int nPntCount)
{
	svxPoint* ptrTXPntList = (svxPoint*)malloc(sizeof(svxPoint) * nPntCount);// X = X(t)
	svxPoint* ptrTYPntList = (svxPoint*)malloc(sizeof(svxPoint) * nPntCount); // Y = Y(t)
	for (int i = 0; i < nPntCount; i++)
	{
		ptrTXPntList[i].x = (double)i;
		ptrTXPntList[i].y = ptrPntList[i].x;
		ptrTYPntList[i].x = (double)i;
		ptrTYPntList[i].y = ptrPntList[i].y;
 	}
	m_dStartX = 0;
	m_dEndX = nPntCount - 1;
	svxPoint* ptrPntListX = new svxPoint[m_iSampleCount];;
	svxPoint* ptrPntListY = new svxPoint[m_iSampleCount];;
	Interpolation_PolynomialBaseFunction(ptrTXPntList, nPntCount, ptrPntListX);
	Interpolation_PolynomialBaseFunction(ptrTYPntList, nPntCount, ptrPntListY);

	m_ptrFittingPntList = new svxPoint[m_iSampleCount];
	for (int i = 0; i < m_iSampleCount; i++)
	{
		m_ptrFittingPntList[i].x = ptrPntListX[i].y;
		m_ptrFittingPntList[i].y = ptrPntListY[i].y;
		m_ptrFittingPntList[i].z = 0;
	}
	delete[] ptrPntListX;
	ptrPntListX = nullptr;
	delete[] ptrPntListY;
	ptrPntListY = nullptr;

	return true;
}

void Calculate2DPolynominalExpress()
{
	
	double dParameter = 0;
	double dStepLength = 1.0 / (m_iSampleCount - 1);

	m_ptrFittingPntList = new svxPoint[m_iSampleCount];

	for (int i = 0; i < m_iSampleCount; i++)
	{
		double x = 0;
		double y = 0;
		double dTExpTemp = 1.0;
		for (int j = 0; j < m_nPolyExp; j++)
		{
			if (0 != j)
			{
				dTExpTemp *= dParameter;//求幂指数
			}
			x += m_matdCoe(j, 0) * dTExpTemp;
			y += m_matdCoe(j, 1) * dTExpTemp;
		}
		dParameter += i * dStepLength;
		m_ptrFittingPntList[i].x = x;
		m_ptrFittingPntList[i].y = y;
		m_ptrFittingPntList[i].z = 0;
	}
}

//void memoryClear()
//{
//	cvxMemFree((void**)&m_ptrFittingPntList);
//}
