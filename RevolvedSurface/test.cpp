void test_intersect()
{
	Point P0{ -1, 2, 0 };
	Point P1{ 1, 2, 0 };
	Point P2{ -1, 1, 0 };
	Point P3{ 1, 1, 0 };
	Point Pt;
	if (Intersect3DLines(P0, P1, P2, P3, Pt))
	{
		std::cout << std::format("Intersection Point: ({}, {}, {})\n", Pt.x, Pt.y, Pt.z);
	}

	P0 = { -1, 1, 0 };
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

	P0 = { -1, 1, 1 };
	P1 = { 0, 0, 1 };
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

void test_vecnormalize()
{
	Point p{ 3, 4, 5 };
	std::cout << std::format("magnitude: {}, normalized: ({}, {}, {})\n", VecNormalize(p), p.x, p.y, p.z);
}

void test_pointtolines()
{
	Point P0{ 0,0,0 };
	Point v{ 0, 1, 0 };
	Point Pt{ 3, 2, 0 };
	Point O;
	PointToLines(P0, v, Pt, O);
	std::cout << std::format("Projected point: ({}, {}, {})\n", O.x, O.y, O.z);

	v = { 0, -1, 0 };
	PointToLines(P0, v, Pt, O);
	std::cout << std::format("Projected point: ({}, {}, {})\n", O.x, O.y, O.z);

	v = { 1, 0, 0 };
	PointToLines(P0, v, Pt, O);
	std::cout << std::format("Projected point: ({}, {}, {})\n", O.x, O.y, O.z);

	v = { -1, 0, 0 };
	PointToLines(P0, v, Pt, O);
	std::cout << std::format("Projected point: ({}, {}, {})\n", O.x, O.y, O.z);

	v = { 0, 0, 1 };
	PointToLines(P0, v, Pt, O);
	std::cout << std::format("Projected point: ({}, {}, {})\n", O.x, O.y, O.z);
}

void test_veccrossprod()
{
	Point v0{ 1, 0, 0 };
	Point v1{ 0, 1, 0 };
	Point v2;
	VecCrossProd(v0, v1, v2);
	std::cout << std::format("cross vector: ({}, {}, {})\n", v2.x, v2.y, v2.z);

	VecCrossProd(v1, v0, v2);
	std::cout << std::format("cross vector: ({}, {}, {})\n", v2.x, v2.y, v2.z);
}