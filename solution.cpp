#include <iostream> 
#include <complex>  
#include <vector> 
#define _USE_MATH_DEFINES
#include <cmath>
#include <math.h>
using namespace std;

double x_0 = -2; double y_0 = -2;
int n_x = 400; double h_x = 0.01; double x_f = x_0 + h_x * n_x;
int n_y = 400; double h_y = 0.01; double y_f = y_0 + h_y * n_y;
int n_t = 3500; double h_t = 0.000005; 
double pi = M_PI; complex<double> i(0, 1); complex<double> zero(0, 0);
vector<vector<complex<double>>> F(n_x);vector<vector<complex<double>>> A(n_x);
vector<vector<complex<double>>> B(n_x);vector<vector<complex<double>>> C(n_x);

double V(const double& x1, const double& y1)
{
	if (-0.1 <= x1 && x1 <= 0)
	{
		if (-2 <= y1 && y1 <= -0.15)
		{
			return 1000;
		}
		if (-0.075 <= y1 && y1 <= 0.075)
		{
			return 1000;
		}
		if (0.15 <= y1 && y1 <= 2)
		{
			return 1000;
		}
	}
	return 0;
}

complex<double> laplacian(const complex<double>& y0, const complex<double>& y1, const complex<double>& y2)
{
	complex<double> y; 
	y = (y2 - (double)2*y1 + y0)/pow(h_x,2);
	return y;
}

complex<double> psi(const double& x1, const double& y1)
{
	double A = 2.3; float sigma = 0.2; float k = 60;
	double x2 = 0.3; double y2 = 0;
	double r = pow(x1-x2, 2) + pow(y1-y2, 2);
	complex<double> y = A*exp(-r/(2*pow(sigma,2)))*(sin(x1*k)+i*cos(x1*k));
	return y;
}

complex<float> Schrodinger(vector<vector<complex<double>>>& F1, vector<vector<complex<double>>>& A1, vector<vector<complex<double>>>& B1, vector<vector<complex<double>>>& C1)
{
	vector<complex<double>> boundary(n_y);
	for (int j = 0; j < n_y; j++)
	{
		boundary[j] = 0;
	}

	for (int n = 0; n < n_t; n++)
	{
		vector<vector<complex<double>>> F2(F1.size()); vector<vector<complex<double>>> A2(A1.size());
		vector<vector<complex<double>>> B2(B1.size()); vector<vector<complex<double>>> C2(C1.size());

		F2[0] = boundary;
		F2[n_x - 1] = boundary;
		for (int k = 1; k < n_x - 1; k++)
		{
			double x_1 = x_0 + k*h_x;
			vector<complex<double>> row(n_y);
			row[0] = 0;
			row[n_y - 1] = 0;
			for (int j = 1; j < n_y - 1; j++)
			{
				double y_1 = y_0 + j*h_y;
				complex<double> k_1, k_2, k_3, k_4;
				complex<double> LX, LY, AX, AY, BX, BY, CX, CY;
				LX = laplacian(F1[k - 1][j], F1[k][j], F1[k + 1][j]);
				LY = laplacian(F1[k][j - 1], F1[k][j], F1[k][j + 1]);

				AX = laplacian(A1[k - 1][j], A1[k][j], A1[k + 1][j]);
				AY = laplacian(A1[k][j - 1], A1[k][j], A1[k][j + 1]);
				
				BX = laplacian(B1[k - 1][j], B1[k][j], B1[k + 1][j]);
				BY = laplacian(B1[k][j - 1], B1[k][j], B1[k][j + 1]);
				
				CX = laplacian(C1[k - 1][j], C1[k][j], C1[k + 1][j]);
				CY = laplacian(C1[k][j - 1], C1[k][j], C1[k][j + 1]);
				
				k_1 = (i/(double)2)*((LX + LY) - V(x_1, y_1)*F1[k][j]);
				k_2 = (i/(double)2)*((LX + LY) + (h_t/2)*(AX + AY) - V(x_1, y_1)*(F1[k][j] + (h_t/2)*A1[k][j]));
				k_3 = (i/(double)2)*((LX + LY) + (h_t/2)*(BX + BY) - V(x_1, y_1)*(F1[k][j] + (h_t/2)*B1[k][j]));
				k_4 = (i/(double)2)*((LX + LY) + h_t*(CX + CY) - V(x_1, y_1)*(F1[k][j] + h_t*C1[k][j]));
				row[j] = F1[k][j] + (h_t/6)*(k_1 + (double)2*k_2 + (double)2*k_3 + k_4);
			}
			F2[k] = row;
		}
		F1 = F2;

		A2[0] = boundary;
		A2[n_x - 1] = boundary;
		for (int k = 1; k < n_x - 1; k++)
		{
			double x_1 = x_0 + k * h_x;
			vector<complex<double>> row(n_y);
			row[0] = 0;
			row[n_y - 1] = 0;
			for (int j = 1; j < n_y - 1; j++)
			{
				double y_1 = y_0 + j * h_y;
				complex<double> LX, LY;
				LX = laplacian(F1[k - 1][j], F1[k][j], F1[k + 1][j]);
				LY = laplacian(F1[k][j - 1], F1[k][j], F1[k][j + 1]);
				row[j] = (i/(double)2)*((LX + LY) - V(x_1, y_1)*F1[k][j]);
			}
			A2[k] = row;
		}
		A1 = A2;

		B2[0] = boundary;
		B2[n_x - 1] = boundary;
		for (int k = 1; k < n_x - 1; k++)
		{
			double x_1 = x_0 + k * h_x;
			vector<complex<double>> row(n_y);
			row[0] = 0;
			row[n_y - 1] = 0;
			for (int j = 1; j < n_y - 1; j++)
			{
				double y_1 = y_0 + j * h_y;
				complex<double> LX, LY, AX, AY;
				LX = laplacian(F1[k - 1][j], F1[k][j], F1[k + 1][j]);
				LY = laplacian(F1[k][j - 1], F1[k][j], F1[k][j + 1]);
				AX = laplacian(A1[k - 1][j], A1[k][j], A1[k + 1][j]);
				AY = laplacian(A1[k][j - 1], A1[k][j], A1[k][j + 1]);
				row[j] = (i/(double)2)*((LX + LY) + (h_t/2)*(AX + AY) - V(x_1, y_1)*(F1[k][j] + (h_t/2)*A1[k][j]));
			}
			B2[k] = row;
		}
		B1 = B2;

		C2[0] = boundary;
		C2[n_x - 1] = boundary;
		for (int k = 1; k < n_x - 1; k++)
		{
			double x_1 = x_0 + k * h_x;
			vector<complex<double>> row(n_y);
			row[0] = 0;
			row[n_y - 1] = 0;
			for (int j = 1; j < n_y - 1; j++)
			{
				double y_1 = y_0 + j * h_y;
				complex<double> LX, LY, LA, BX, BY;
				LX = laplacian(F1[k - 1][j], F1[k][j], F1[k + 1][j]);
				LY = laplacian(F1[k][j - 1], F1[k][j], F1[k][j + 1]);
				BX = laplacian(B1[k - 1][j], B1[k][j], B1[k + 1][j]);
				BY = laplacian(B1[k][j - 1], B1[k][j], B1[k][j + 1]);
				row[j] = (i/(double)2)*((LX + LY) + (h_t/2)*(BX + BY) - V(x_1, y_1)*(F1[k][j] + (h_t/2)*B1[k][j]));
			}
			C2[k] = row;
		}
		C1 = C2;
	}
}

int main() 
{	
	for (int k = 0; k < n_x; k++)
	{
		vector<complex<double>> row(n_y);
		double x_1 = x_0 + k*h_x;
		for (int j = 0; j < n_y; j++)
		{
			double y_1 = y_0 + j*h_y;
			row[j] = psi(x_1, y_1);
		}
		F[k] = row;
	}

	
	vector<complex<double>> boundary(n_y);
	for (int j = 0; j < n_y; j++)
	{
		boundary[j] = 0;
	}

	A[0] = boundary;
	A[n_x - 1] = boundary;
	for (int k = 1; k < n_x - 1; k++)
	{
		double x_1 = x_0 + k*h_x;
		vector<complex<double>> row(n_y);
		row[0] = 0;
		row[n_y - 1] = 0;
		for (int j = 1; j < n_y - 1; j++)
		{
			double y_1 = y_0 + j*h_y;
			complex<double> LX, LY;
			LX = laplacian(F[k-1][j], F[k][j], F[k+1][j]);
			LY = laplacian(F[k][j-1], F[k][j], F[k][j+1]);
			row[j] = (i/(double)2)*((LX + LY) - V(x_1, y_1)*F[k][j]);
		}
		A[k] = row;
	}

	B[0] = boundary;
	B[n_x - 1] = boundary;
	for (int k = 1; k < n_x - 1; k++)
	{
		double x_1 = x_0 + k * h_x;
		vector<complex<double>> row(n_y);
		row[0] = 0;
		row[n_y - 1] = 0;
		for (int j = 1; j < n_y - 1; j++)
		{
			double y_1 = y_0 + j * h_y;
			complex<double> LX, LY, LA, AX, AY;
			LX = laplacian(F[k - 1][j], F[k][j], F[k + 1][j]);
			LY = laplacian(F[k][j - 1], F[k][j], F[k][j + 1]);
			AX = laplacian(A[k - 1][j], A[k][j], A[k + 1][j]);
			AY = laplacian(A[k][j - 1], A[k][j], A[k][j + 1]);
			row[j] = (i/(double)2)*((LX + LY) + (h_t/2)*(AX + AY) - V(x_1, y_1)*(F[k][j] + (h_t/2)*A[k][j]));
		}
		B[k] = row;
	}

	C[0] = boundary;
	C[n_x - 1] = boundary;
	for (int k = 1; k < n_x - 1; k++)
	{
		double x_1 = x_0 + k * h_x;
		vector<complex<double>> row(n_y);
		row[0] = 0;
		row[n_y - 1] = 0;
		for (int j = 1; j < n_y - 1; j++)
		{
			double y_1 = y_0 + j * h_y;
			complex<double> LX, LY, LA, BX, BY;
			LX = laplacian(F[k - 1][j], F[k][j], F[k + 1][j]);
			LY = laplacian(F[k][j - 1], F[k][j], F[k][j + 1]);
			BX = laplacian(B[k - 1][j], B[k][j], B[k + 1][j]);
			BY = laplacian(B[k][j - 1], B[k][j], B[k][j + 1]);
			row[j] = (i/(double)2)*((LX + LY) + (h_t/2)*(BX + BY) - V(x_1, y_1)*(F[k][j] + (h_t/2)*B[k][j]));
		}
		C[k] = row;
	}
	Schrodinger(F, A, B, C);
}
