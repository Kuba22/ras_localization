#pragma once

#include <stdio.h>
#include <stdlib.h>
#include <armadillo>

using namespace arma;
using namespace std;

template<typename T>
T mod(T a, double n)
{
	return a - floor(a / n)*n;
}

mat readLines(string map_file)
{
	mat lines(1, 4);

	ifstream map_fs;
	map_fs.open(map_file.c_str());
	if (!map_fs.is_open()) {
	}

	string line;
	int wall_id = 0;
	while (getline(map_fs, line)) {
		if (line[0] == '#') {
			// comment -> skip
			continue;
		}

		double max_num = std::numeric_limits<double>::max();
		double x1 = max_num,
			x2 = max_num,
			y1 = max_num,
			y2 = max_num;

		std::istringstream line_stream(line);

		line_stream >> x1 >> y1 >> x2 >> y2;
		rowvec linevec({ x1, y1, x2, y2 });
		lines = join_cols(lines, linevec);
	}

	lines.shed_row(0);
	return lines;
}

vec crossProduct(vec a, vec b) {
	vec ret(3);
	if (a.size() == 2)
	{
		ret[0] = a[1] - b[1];
		ret[1] = b[0] - a[0];
		ret[2] = a[0] * b[1] - a[1] * b[0];
	}
	else
	{
		ret[0] = a[1] * b[2] - a[2] * b[1];
		ret[1] = a[2] * b[0] - a[0] * b[2];
		ret[2] = a[0] * b[1] - a[1] * b[0];
	}
	return ret;
}

double roundDouble(double val) {
	double x = pow(10, 6);
	return round(val * x) / x;
}

bool notProperX(rowvec line, vec X, vec P, vec R) {
	double x = roundDouble(X[0] / X[2]);
	double y = roundDouble(X[1] / X[2]);
	double px = roundDouble(P[0]);
	double py = roundDouble(P[1]);
	double rx = roundDouble(R[0]);
	double ry = roundDouble(R[1]);
	double l1x = roundDouble(line[0]);
	double l1y = roundDouble(line[1]);
	double l2x = roundDouble(line[2]);
	double l2y = roundDouble(line[3]);
	if (x > l1x && x > l2x
		|| x < l1x && x < l2x
		|| y > l1y && y > l2y
		|| y < l1y && y < l2y
		|| x > rx && x > px
		|| x < rx && x < px
		|| y > ry && y > py
		|| y < ry && y < py)
		return true;
	return false;
}

double getRange(mat lines, vec R, double phi) {
	double minR = 100.0;
	double theta = R[2];
	vec P({ 100.0*cos(phi + theta) + R(0), 100.0*sin(phi + theta) + R(1) });
	for (int i = 0; i < lines.n_rows; i++) {
		vec ll = crossProduct(vec({ lines(i, 0), lines(i, 1) }), vec({ lines(i, 2), lines(i, 3) }));
		vec RP = crossProduct(R(span(0, 1)), P);
		vec x = crossProduct(ll, RP);
		if (notProperX(lines.row(i), x, P, R))
			continue;
		double r = sqrt((R[0] - x[0] / x[2])*(R[0] - x[0] / x[2]) + (R[1] - x[1] / x[2])*(R[1] - x[1] / x[2]));
		if (r < minR) {
			minR = r;
		}
	}
	return minR;
}