#include <ros/ros.h>
#include <phidgets/motor_encoder.h>
#include <geometry_msgs/PoseStamped.h>

const double r = 0.036;
const double b = 0.242;
const int ticks_per_rev = 360;
phidgets::motor_encoder encoder_left, encoder_right;

void EncoderLeftCallback(const phidgets::motor_encoder::ConstPtr &msg)
{
    encoder_left = *msg;
}

void EncoderRightCallback(const phidgets::motor_encoder::ConstPtr &msg)
{
    encoder_right = *msg;
}

int main(int argc, char *argv[])
{
    ros::init(argc, argv, "odometry");
    ros::NodeHandle nh;

    ros::Subscriber encoder_left_sub = nh.subscribe("/motorcontrol/encoder_left", 1, &EncoderLeftCallback);
    ros::Subscriber encoder_right_sub = nh.subscribe("/motorcontrol/encoder_right", 1, &EncoderRightCallback);

    ros::Publisher pose_pub = nh.advertise<geometry_msgs::PoseStamped>("localization/odometry_pose", 1);

    int frequency = 30; // this could be a parameter
    ros::Rate loop_rate(frequency);
    double theta, x, y; // this will have to be initialized with the use of the laser scan
    while(nh.ok())
    {
	double w = (encoder_right.count_change - encoder_left.count_change)/b*2*3.1415*ticks_per_rev;
	double v = (encoder_right.count_change + encoder_left.count_change)*r*0.5;
	theta += w;
	x += v*cos(theta);
	y += v*sin(theta);

	geometry_msgs::PoseStamped pose_stamped;
	pose_stamped.header.stamp = ros::Time::now();
	pose_stamped.header.frame_id = "odom";
	pose_stamped.pose.position.x = x;
	pose_stamped.pose.position.y = y;
	pose_stamped.pose.orientation.z = sin(theta*0.5);
	pose_stamped.pose.orientation.w = cos(theta*0.5);
	
	pose_pub.publish(pose_stamped);

	ros::spinOnce();
	loop_rate.sleep();
    }
}
