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
#include <numbers>
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
	Original code from
	http://paulbourke.net/geometry/pointlineplane/
    Calculate the line segment PaPb that is the shortest route between
    two lines P1P2 and P3P4. Calculate also the values of mua and mub where
        Pa = P1 + mua (P2 - P1) (1)
        Pb = P3 + mub (P4 - P3) (2)
    Return true if Pa == Pb

	In this algorithm
	Input, (the point, tangent) pairs: (P0, T0), (P2, T2)
	Output: Pt => intersection point if return value is true

	Therefore, (1) and (2) become
	Pa = P0 + mua * T0
	Pb = P2 + mub * T2

	Pa == Pb => inteersection point
	*/

	Point p13, p43, p21;

	p13.x = P0.x - P2.x; // p13.x = P1.x - P3.x;
	p13.y = P0.y - P2.y; // p13.y = P1.y - P3.y;
	p13.z = P0.z - P2.z; // p13.z = P1.z - P3.z;

	// p43 is T2
	p43 = T2;

	//p43.x = P4.x - P3.x;
	//p43.y = P4.y - P3.y;
	//p43.z = P4.z - P3.z;

	if (std::abs(p43.x) < epsilon && std::abs(p43.y) < epsilon && std::abs(p43.z) < epsilon)
	{
		return false;
	}

	// p21 is T0
	p21 = T0;

	//p21.x = P2.x - P1.x;
	//p21.y = P2.y - P1.y;
	//p21.z = P2.z - P1.z;

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
	pa.x = P0.x + mua * p21.x; //pa.x = P1.x + mua * p21.x;
	pa.y = P0.y + mua * p21.y; //pa.y = P1.y + mua * p21.y;
	pa.z = P0.z + mua * p21.z; //pa.z = P1.z + mua * p21.z;

	pb.x = P2.x + mub * p43.x; // pb.x = P3.x + mub * p43.x;
	pb.y = P2.y + mub * p43.y; // pb.y = P3.y + mub * p43.y;
	pb.z = P2.z + mub * p43.z; // pb.z = P3.z + mub * p43.z;

	if (std::abs(pa.x - pb.x) < epsilon && std::abs(pa.y - pb.y) < epsilon && std::abs(pa.z - pb.z) < epsilon)
	{
		Pt = pa;
		return true;
	}

	return false;
}


void MakeRevolvedSurface(Point& S, Point& T, double theta, int m, std::vector<Point>& Pj, std::vector<double>& wj,
	int& n, std::vector<double>& U, std::vector<std::vector<Point>>& Pij, std::vector<std::vector<double>>& wij)
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
		U.resize(6);
	}
	else if (theta <= 180.0)
	{
		narcs = 2;
		U.resize(8);
		U[3] = U[4] = 0.5;
	}
	else if (theta <= 270.0)
	{
		narcs = 3;
		U.resize(10);
		U[3] = U[4] = 1.0 / 3.0;
		U[5] = U[6] = 2.0 / 3.0;
	}
	else
	{
		narcs = 4;
		U.resize(12);
		U[3] = U[4] = 0.25;
		U[5] = U[6] = 0.5;
		U[7] = U[8] = 0.75;
	}

	dtheta = theta / narcs;
	int J{ 3 + 2 * (narcs - 1) };
	// load end knots
	for (int i{}; i < 3; ++J, ++i)
	{
		U[i] = 0.0;
		U[J] = 1.0;
	}

	n = 2 * narcs;

	double wm{ cos(dtheta * std::numbers::pi / 360.0) }; // dtheta / 2.0 is base angle
	double angle{ 0.0 }; // compute sine and cosine only once
	std::vector<double> cosines(narcs + 1), sines(narcs + 1);

	for (int i{ 1 }; i <= narcs; ++i)
	{
		angle += dtheta;
		double radAngle{ angle * std::numbers::pi / 180.0 };
		cosines[i] = std::cos(radAngle);
		sines[i] = std::sin(radAngle);
	}

	Point X, Y, O, P0, P2, T0, T2;
	double r{};
	int index{};

	Pij.resize(n + 1);
	wij.resize(n + 1);
	for (int i{ 0 }; i <= n; ++i)
	{
		Pij[i].resize(m + 1);
		wij[i].resize(m + 1);
	}

	for (int j{ 1 }; j <= m; ++j)
	{
		// loop and compute each u row of ctrl pts and weights
		PointToLines(S, T, Pj[j], O);
		X.x = Pj[j].x - O.x;
		X.y = Pj[j].y - O.y;
		X.z = Pj[j].z - O.z;
		r = VecNormalize(X);
		VecCrossProd(T, X, Y);

		// Initialize first ctrl pt and weight
		P0 = Pj[j];
		Pij[0][j] = P0;
		wij[0][j] = wj[j];
		T0 = Y;
		index = 0;
		angle = 0.0;

		for (int i{ 1 }; i <= narcs; ++i) // compute u row
		{
			P2.x = O.x + r * cosines[i] * X.x + r * sines[i] * Y.x;
			P2.y = O.y + r * cosines[i] * X.y + r * sines[i] * Y.y;
			P2.z = O.z + r * cosines[i] * X.z + r * sines[i] * Y.z;

			Pij[index + 2][j] = P2;
			wij[index + 2][j] = wj[j];

			T2.x = -sines[i] * X.x + cosines[i] * Y.x;
			T2.y = -sines[i] * X.y + cosines[i] * Y.y;
			T2.z = -sines[i] * X.z + cosines[i] * Y.z;

			//if (!Intersect3DLines(P0, P2, T0, T2, Pij[index + 1][j]))
			if (!Intersect3DLines(P0, T0, P2, T2, Pij[index + 1][j]))
			{
				std::cout << "Intersect3DLines failed\n";
			}
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



int main()
{
	//void MakeRevolvedSurface(Point& S, Point& T, double theta, int m, std::vector<Point>& Pj, std::vector<double>& wj,
	//int& n, std::vector<double>& U, std::vector<std::vector<Point>>& Pij, std::vector<std::vector<double>>& wij)
	//S, T, theta, m, Pj[], wj[]
	// Output: n, U, Pij[][], wij[][]

	Point startPt{ 0,0,0 };
	Point directionVec{ 0,0,1 };
	double angle{ 360.0 };
	int m{8};
	std::vector<Point> controlPts{ {1, 0, 2}, {1, 0, 3}, {2, 0, 3}, {3, 0, 3}, {3, 0, 2}, {3, 0, 1}, {2, 0, 1}, {1, 0, 1}, {1, 0, 2} };
	std::vector<double> wj{1, 0.7071, 1, 0.7071, 1, 0.7071, 1, 0.7071, 1};

	int n{};
	std::vector<double> U;
	std::vector<std::vector<Point>> Pij;
	std::vector<std::vector<double>> wij;

	MakeRevolvedSurface(startPt, directionVec, angle, m, controlPts, wj, n, U, Pij, wij);
}