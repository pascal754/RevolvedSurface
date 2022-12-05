module RevolvedSurface;

import std;
import Point3D;

constexpr double epsilon{ 1e-6 };

double VecNormalize(Point3D& pt)
{
    // normalize vector
    // return magnitude of the original vector

    double magnitude{ std::hypot(pt.x, pt.y, pt.z) };
    pt.x /= magnitude;
    pt.y /= magnitude;
    pt.z /= magnitude;

    return magnitude;

}

void PointToLines(const Point3D& S, const Point3D& T, const Point3D& Pt, Point3D& O)
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

void VecCrossProd(const Point3D& a, const Point3D& b, Point3D& c)
{
    c.x = a.y * b.z - a.z * b.y;
    c.y = a.z * b.x - a.x * b.z;
    c.z = a.x * b.y - a.y * b.x;
}

bool Intersect3DLines(const Point3D& P0, const Point3D& T0, const Point3D& P2, const Point3D& T2, Point3D& Pt)
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

    Point3D p13, p43, p21;

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

    double d1343{ p13.x * p43.x + p13.y * p43.y + p13.z * p43.z };
    double d4321{ p43.x * p21.x + p43.y * p21.y + p43.z * p21.z };
    double d1321{ p13.x * p21.x + p13.y * p21.y + p13.z * p21.z };
    double d4343{ p43.x * p43.x + p43.y * p43.y + p43.z * p43.z };
    double d2121{ p21.x * p21.x + p21.y * p21.y + p21.z * p21.z };

    double denom{ d2121 * d4343 - d4321 * d4321 };

    if (std::abs(denom) < epsilon)
    {
        return false;
    }

    double numer{ d1343 * d4321 - d1321 * d4343 };

    double mua{ numer / denom };
    double mub{ (d1343 + d4321 * mua) / d4343 };

    Point3D pa, pb;
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


void MakeRevolvedSurface(const Point3D& S, const Point3D& T, double theta, int m, const std::vector<Point3D>& Pj, const std::vector<double>& wj,
    int& n, std::vector<double>& U, std::vector<std::vector<Point3D>>& Pij, std::vector<std::vector<double>>& wij)
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

    double wm{ std::cos(dtheta * std::numbers::pi / 360.0) }; // dtheta / 2.0 is base angle
    double angle{ 0.0 }; // compute sine and cosine only once
    std::vector<double> cosines(narcs + 1), sines(narcs + 1);

    for (int i{ 1 }; i <= narcs; ++i)
    {
        angle += dtheta;
        double radAngle{ angle * std::numbers::pi / 180.0 };
        cosines[i] = std::cos(radAngle);
        sines[i] = std::sin(radAngle);
    }

    Point3D X, Y, O, P0, P2, T0, T2;
    double r{};
    int index{};

    Pij.resize(n + 1);
    wij.resize(n + 1);
    for (int i{ 0 }; i <= n; ++i)
    {
        Pij[i].resize(m + 1);
        wij[i].resize(m + 1);
    }

    for (int j{ 0 }; j <= m; ++j)
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

void make_nurbs_circle(const Point3D& O, const Point3D& X, const Point3D& Y,
    double r, double ths, double the, int& n, std::vector<double>& U, std::vector<Point3D>& Pw)
{
    // Algorithm A7.1 pp. 308
    // Create arbitrary NURBS circular arc
    // Input: O, X, Y, r, ths, the
    // O: center of circle (origin of local coordinate system)
    // X: unit length vector lying in the plane of definition of the circle
    // Y: unit length vector in the plane of definition of the circle, and orthogonal to X
    // r: radius
    // ths, the: start and end angles, measured with respect to X
    // 
    // Output: n, U, Pw

    if (the < ths)
    {
        the += 360.0;
    }

    double theta{ the - ths };
    int narcs{};

    if (theta <= 90.0)
    {
        narcs = 1; // get number of arcs
    }
    else if (theta <= 180.0)
    {
        narcs = 2;
    }
    else if (theta <= 270.0)
    {
        narcs = 3;
    }
    else
    {
        narcs = 4;
    }

    double dtheta{ theta / narcs };
    n = 2 * narcs; // n + 1: control points
    double w1{ std::cos(std::numbers::pi * dtheta / 360.0) }; // dtheta / 2 is base angle

    Pw.resize(n + 1);

    Point3D P0;

    double ths_rad{ ths * std::numbers::pi / 180.0 };
    P0.x = O.x + r * std::cos(ths_rad) * X.x + r * std::sin(ths_rad) * Y.x;
    P0.y = O.y + r * std::cos(ths_rad) * X.y + r * std::sin(ths_rad) * Y.y;
    P0.z = O.z + r * std::cos(ths_rad) * X.z + r * std::sin(ths_rad) * Y.z;

    Point3D T0;

    T0.x = -std::sin(ths_rad) * X.x + std::cos(ths_rad) * Y.x; // Initialize start values
    T0.y = -std::sin(ths_rad) * X.y + std::cos(ths_rad) * Y.y;
    T0.z = -std::sin(ths_rad) * X.z + std::cos(ths_rad) * Y.z;

    Pw[0] = P0;

    int index{};
    double angle{ ths };
    Point3D P1;
    Point3D P2;
    Point3D T2;

    for (int i{ 1 }; i <= narcs; ++i) // create narcs segments
    {
        angle += dtheta;
        double angle_rad{ angle * std::numbers::pi / 180.0 };
        P2.x = O.x + r * std::cos(angle_rad) * X.x + r * std::sin(angle_rad) * Y.x;
        P2.y = O.y + r * std::cos(angle_rad) * X.y + r * std::sin(angle_rad) * Y.y;
        P2.z = O.z + r * std::cos(angle_rad) * X.z + r * std::sin(angle_rad) * Y.z;

        Pw[index + 2] = P2;

        T2.x = -std::sin(angle_rad) * X.x + std::cos(angle_rad) * Y.x;
        T2.y = -std::sin(angle_rad) * X.y + std::cos(angle_rad) * Y.y;
        T2.z = -std::sin(angle_rad) * X.z + std::cos(angle_rad) * Y.z;

        if (!Intersect3DLines(P0, T0, P2, T2, P1))
        {
            std::cerr << "Intersect3DLines failed\n";
        }

        Pw[index + 1].x = w1 * P1.x;
        Pw[index + 1].y = w1 * P1.y;
        Pw[index + 1].z = w1 * P1.z;

        index += 2;

        if (i < narcs)
        {
            P0 = P2;
            T0 = T2;
        }
    }

    switch (narcs)
    {
    case 1:
        U.resize(6);
        break;
    case 2:
        U.resize(8);
        U[3] = U[4] = 0.5;
        break;
    case 3:
        U.resize(10);
        U[3] = U[4] = 1.0 / 3.0;
        U[5] = U[6] = 2.0 / 3.0;
        break;
    case 4:
        U.resize(12);
        U[3] = U[4] = 0.25;
        U[5] = U[6] = 0.5;
        U[7] = U[8] = 0.75;
        break;
    default:
        std::cerr << "error\n";
        break;
    }

    int J{ 2 * narcs + 1 }; // load the knot vector
    for (int i{}; i < 3; ++i)
    {
        U[i] = 0.0;
        U[i + J] = 1.0;
    }
}
