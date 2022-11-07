// October 28, 2022
// by SM



//import <iostream>;
//import <array>;
//import <vector>;
//import <cmath>;

#include <iostream>
#include <array>
#include <vector>
#include <format>
#include <cmath>

const double epsilon{ 1e-6 };

struct Point
{
	double x{};
	double y{};
	double z{};
};

double VecNormalize(Point& pt)
{
	// normalize vector
	// return magnitude of the original vector

	double magnitude{ std::hypot(pt.x, pt.y, pt.z) };
	pt.x /= magnitude;
	pt.y /= magnitude;
	pt.z /= magnitude;

	return magnitude;

}

void PointToLines(Point& S, Point& T, Point& Pt, Point& O)
{
	// project Pt on the line passing S and pointing T
	// S: starting point
	// T: unit direction vector
	// O: Projected Point

	double dotProduct{ S.x * T.x + S.y * T.y + S.z * T.z };
	O.x = S.x + (Pt.x - S.x) * T.x * T.x;
	O.y = S.y + (Pt.y - S.y) * T.y * T.y;
	O.z = S.z + (Pt.z - S.z) * T.z * T.z;
}

void VecCrossProd(const Point& a, const Point& b, Point& c)
{
	c.x = a.y * b.z - a.z * b.y;
	c.y = a.z * b.x - a.x * b.z;
	c.z = a.x * b.y - a.y * b.x;
}

bool Intersect3DLines(const Point& P0, const Point& T0, const Point& P2, const Point& T2, Point& Pt)
{
	/*
    Calculate the line segment PaPb that is the shortest route between
    two lines P1P2 and P3P4. Calculate also the values of mua and mub where
        Pa = P1 + mua (P2 - P1)
        Pb = P3 + mub (P4 - P3)
    Return true if Pa == Pb

	Input: P0, T0, P2, T2
	Output: Pt => intersection point if return value is true
	*/

	Point p13, p43, p21;

	p13.x = P0.x - P2.x;
	p13.y = P0.y - P2.y;
	p13.z = P0.z - P2.z;

	p43.x = T2.x - P2.x;
	p43.y = T2.y - P2.y;
	p43.z = T2.z - P2.z;

	if (std::abs(p43.x) < epsilon && std::abs(p43.y) < epsilon && std::abs(p43.z) < epsilon)
	{
		return false;
	}

	p21.x = T0.x - P0.x;
	p21.y = T0.y - P0.y;
	p21.z = T0.z - P0.z;

	if (std::abs(p21.x) < epsilon && std::abs(p21.y) < epsilon && std::abs(p21.z) < epsilon)
	{
		return false;
	}
	

	double d1343 { p13.x * p43.x + p13.y * p43.y + p13.z * p43.z};
	double d4321 { p43.x * p21.x + p43.y * p21.y + p43.z * p21.z};
	double d1321 { p13.x * p21.x + p13.y * p21.y + p13.z * p21.z};
	double d4343 { p43.x * p43.x + p43.y * p43.y + p43.z * p43.z};
	double d2121{ p21.x * p21.x + p21.y * p21.y + p21.z * p21.z };

	double denom{ d2121 * d4343 - d4321 * d4321 };

	if (std::abs(denom) < epsilon)
	{
		return false;
	}

	double numer{ d1343 * d4321 - d1321 * d4343 };

	double mua{ numer / denom };
	double mub{ (d1343 + d4321 * mua) / d4343 };

	Point pa, pb;
	pa.x = P0.x + mua * p21.x;
	pa.y = P0.y + mua * p21.y;
	pa.z = P0.z + mua * p21.z;

	pb.x = P2.x + mub * p43.x;
	pb.y = P2.y + mub * p43.y;
	pb.z = P2.z + mub * p43.z;

	if (std::abs(pa.x - pb.x) < epsilon && std::abs(pa.y - pb.y) < epsilon && std::abs(pa.z - pb.z) < epsilon)
	{
		Pt = pa;
		return true;
	}

	return false;
}


void MakeRevolvedSurface(Point& S, Point& T, double theta, int m, std::vector<Point>& Pj, std::vector<double>& wj,
	int n, std::vector<double>& U, std::vector<std::vector<Point>>& Pij, std::vector<std::vector<double>>& wij)
{
	// Algorithm A8.1
	// Create NURBS surface of revolution
	// Input: S, T, theta, m, Pj[], wj[]
	// Output: n, U, Pij[][], wij[][]

	int narcs{};
	double dtheta{};

	if (theta <= 90.0)
	{
		narcs = 1;
	}
	else if (theta <= 180.0)
	{
		narcs = 2; U[3] = U[4] = 0.5;
	}
	else if (theta <= 270.0)
	{
		narcs = 3;
		U[3] = U[4] = 1.0 / 3.0;
		U[5] = U[6] = 2.0 / 3.0;
	}
	else
	{
		narcs = 4;
		U[3] = U[4] = 0.25;
		U[5] = U[6] = 0.5;
		U[7] = U[8] = 0.75;
	}

	dtheta = theta / narcs;
	int J{ 3 + 2 * (narcs - 1) }; // load end knots
	for (int i{}; i < 3; ++J, ++i)
	{
		U[i] = 0.0;
		U[J] = 1.0;
	}
	n = 2 * narcs;
	double wm{ cos(dtheta / 2.0) }; // dtheta / 2.0 is base angle
	double angle{ 0.0 }; // compute sine and cosine only once
	std::vector<double> cosines, sines;

	for (int i{ 1 }; i <= narcs; ++i)
	{
		angle += dtheta;
		cosines.push_back(cos(angle));
		sines.push_back(sin(angle));
	}

	Point X, Y, O, P0, P2, T0, T2;
	double r{};
	int index{};

	for (int j{ 1 }; j <= m; ++j)
	{
		// loop and compute each u row of ctrl pts and weights
		PointToLines(S, T, Pj[j], O);
		X.x = Pj[j].x - O.x;
		X.y = Pj[j].y - O.y;
		X.z = Pj[j].z - O.z;
		r = VecNormalize(X);
		VecCrossProd(T, X, Y);
		Pij[0][j] = P0 = Pj[j]; // Initialize first
		wij[0][j] = wj[j]; // ctrl pt and weight
		T0 = Y;
		index = 0;
		angle = 0.0;
		for (int i{ 1 }; i <= narcs; ++i) // compute u row
		{
			P2.x = O.x + r * cosines[i] * X.x + r * sines[i] * Y.x;
			P2.y = O.y + r * cosines[i] * X.y + r * sines[i] * Y.y;
			P2.z = O.z + r * cosines[i] * X.z + r * sines[i] * Y.z;
			Pij[index + 2][j] = P2;
			wij[index + 2][i] = wj[j];
			T2.x = -sines[i] * X.x + cosines[i] * Y.x;
			T2.y = -sines[i] * X.y + cosines[i] * Y.y;
			T2.z = -sines[i] * X.z + cosines[i] * Y.z;
			Intersect3DLines(P0, T0, P2, T2, Pij[index + 1][j]);
			wij[index + 1][j] = wm * wj[j];
			index += 2;
			if (i < narcs)
			{
				P0 = P2;
				T0 = T2;
			}
		}
	}
}

void test_intersect()
{
	Point P0{-1, 2, 0};
	Point P1{ 1, 2, 0 };
	Point P2{ -1, 1, 0 };
	Point P3{ 1, 1, 0 };
	Point Pt;
	if (Intersect3DLines(P0, P1, P2, P3, Pt))
	{
		std::cout << std::format("Intersection Point: ({}, {}, {})\n", Pt.x, Pt.y, Pt.z);
	}

	P0 = {-1, 1, 0};
	P1 = { 1, -1, 0 };
	P2 = { -1, -1, 0 };
	P3 = { 1, 1, 0 };
	if (Intersect3DLines(P0, P1, P2, P3, Pt))
	{
		std::cout << std::format("Intersection Point: ({}, {}, {})\n", Pt.x, Pt.y, Pt.z);
	}

	P0 = { -1, 1, 0 };
	P1 = { 0, 0, 0 };
	P2 = { -1, 0, 0 };
	P3 = { 0, 1, 0 };
	if (Intersect3DLines(P0, P1, P2, P3, Pt))
	{
		std::cout << std::format("Intersection Point: ({}, {}, {})\n", Pt.x, Pt.y, Pt.z);
	}

	P0 = { -1, 1, 0 };
	P1 = { 0, 0, 0 };
	P2 = { -1, 1, 0 };
	P3 = { 0, 0, 0 };
	if (Intersect3DLines(P0, P1, P2, P3, Pt))
	{
		std::cout << std::format("Intersection Point: ({}, {}, {})\n", Pt.x, Pt.y, Pt.z);
	}
}

int main()
{
	test_intersect();
}