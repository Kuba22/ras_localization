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
		throw std::runtime_error("map file not open");
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
		join_cols(lines, linevec);
	}

	lines.shed_row(lines.n_rows - 1);
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

vec perpendicularThroughPoint(vec line, vec point){
	vec ret(3);
	ret[0] = -line[1];
	ret[1] = line[0];
	ret[2] = line[1] * point[0] - line[0] * point[1];
	return ret;
}

vec measurement(vec line, vec robot)
{
	vec l1({ line[0], line[1] });
	vec l2({ line[2], line[3] });
	vec ll = crossProduct(l1, l2);
	vec pl = perpendicularThroughPoint(ll, robot);
	vec r = crossProduct(ll, pl);
	double rho = sqrt((r[0] / r[2] - robot[0])*(r[0] / r[2] - robot[0]) + (r[1] / r[2] - robot[1])*(r[1] / r[2] - robot[1]));
	double phi = atan2(r[1] / r[2] - robot[1], r[0] / r[2] - robot[0]);
	double theta = phi - robot[2];
	return vec({ rho, theta });
}
