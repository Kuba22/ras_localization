#include <ros/ros.h>
#include <phidgets/motor_encoder.h>
#include <geometry_msgs/PoseStamped.h>

const double r = 0.036;
const double b = 0.242;
double ticks_per_rev;
phidgets::motor_encoder encoder_left, encoder_right;

int accumulated_left = 0, accumulated_right = 0;

void EncoderLeftCallback(const phidgets::motor_encoder::ConstPtr &msg)
{
    encoder_left = *msg;
    accumulated_left += encoder_left.count_change;
}

void EncoderRightCallback(const phidgets::motor_encoder::ConstPtr &msg)
{
    encoder_right = *msg;
    accumulated_right += encoder_right.count_change;
}

int main(int argc, char *argv[])
{
    ros::init(argc, argv, "odometry");
    ros::NodeHandle nh;

    int frequency;
    nh.param("frequency", frequency, 30);
    nh.param("encoder_res", ticks_per_rev, 900.0);

    ros::Subscriber encoder_left_sub = nh.subscribe("/motorcontrol/encoder_left", 100, &EncoderLeftCallback);
    ros::Subscriber encoder_right_sub = nh.subscribe("/motorcontrol/encoder_right", 100, &EncoderRightCallback);

    ros::Publisher pose_pub = nh.advertise<geometry_msgs::PoseStamped>("localization/odometry_pose", 1);

    ros::Rate loop_rate(frequency);
    double theta = 0.0, x = 0.0, y = 0.0; // this will have to be initialized with the use of the laser scan
    while(nh.ok())
    {
        double w_right = ((double)accumulated_right) * 2 * 3.14159 / ticks_per_rev * frequency;
        double w_left = ((double)accumulated_left) * 2 * 3.14159 / ticks_per_rev * frequency;
        double w = (w_right - w_left)*r/b;
        double v = (w_right + w_left)*r*0.5;
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

        accumulated_left = accumulated_right = 0;

        ros::spinOnce();
        loop_rate.sleep();
    }
}
