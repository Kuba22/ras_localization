#include <ros/ros.h>
#include <armadillo>
#include "pf.h"
#include <geometry_msgs/Twist.h>
#include <ras_line_detector/LineSegmentList.h>
#include <ras_particle_filter/Particles.h>

class ParticleFilterNode
{
public:
    ros::NodeHandle nh;
    ros::Subscriber twist_sub, lines_sub;
    ros::Publisher particles_pub;
    ParticleFilter pf;

    geometry_msgs::Twist last_twist;
    ras_line_detector::LineSegmentList line_segments;

    void TwistCallback(const geometry_msgs::Twist::ConstPtr& msg){
        last_twist = *msg;
        // take into account the time intervals between the messages
    }

    void LinesCallback(const ras_line_detector::LineSegmentList::ConstPtr& msg){
        line_segments = *msg;
    }

    ParticleFilterNode()
    {
        lines_sub = nh.subscribe("/lines", 1, &ParticleFilterNode::LinesCallback, this);
        twist_sub = nh.subscribe("/localization/odometry_twist", 20, &ParticleFilterNode::TwistCallback, this);
        particles_pub = nh.advertise<ras_particle_filter::Particles>("/particles", 1);

    }
};

int main(int argc, char *argv[])
{
	ros::init(argc, argv, "particle_filter");
	ROS_INFO("Particle Filter Node");

	ParticleFilterNode node;

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
