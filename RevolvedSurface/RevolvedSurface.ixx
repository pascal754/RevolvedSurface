export module RevolvedSurface;

import std;
import Point3D;

export double VecNormalize(Point3D& pt);
export void PointToLines(const Point3D& S, const Point3D& T, const Point3D& Pt, Point3D& O);
export void VecCrossProd(const Point3D& a, const Point3D& b, Point3D& c);
export bool Intersect3DLines(const Point3D& P0, const Point3D& T0, const Point3D& P2, const Point3D& T2, Point3D& Pt);
export void MakeRevolvedSurface(const Point3D& S, const Point3D& T, double theta, int m, const std::vector<Point3D>& Pj, const std::vector<double>& wj,
    int& n, std::vector<double>& U, std::vector<std::vector<Point3D>>& Pij, std::vector<std::vector<double>>& wij);
export void make_nurbs_circle(const Point3D& O, const Point3D& X, const Point3D& Y,
    double r, double ths, double the, int& n, std::vector<double>& U, std::vector<Point3D>& Pw);
