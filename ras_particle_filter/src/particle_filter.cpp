#include <ros/ros.h>
#include <armadillo>

int main(int argc, char *argv[])
{
	ros::init(argc, argv, "particle_filter");
	ROS_INFO("Particle Filter Node");

	arma::mat m(2,2);
	m(0,0) = 1;
	m(0,1) = 2;
	m(1,0) = 3;
	m(1,1) = 4;
	arma::vec v(2);
	v(0) = 10; v(1) = 10;
	arma::vec mv = m*v;
	ROS_INFO("(%f, %f)", mv(0), mv(1));
}
