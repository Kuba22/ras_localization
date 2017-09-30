#include <string>
#include <ros/ros.h>
#include <sensor_msgs/JointState.h>
#include <tf/transform_broadcaster.h>

int main(int argc, char** argv) {
    ros::init(argc, argv, "state_publisher");
    ros::NodeHandle n;
    //ros::Publisher joint_pub = n.advertise<sensor_msgs::JointState>("joint_states", 1);
    tf::TransformBroadcaster broadcaster;
    ros::Rate loop_rate(30);

    const double degree = M_PI/180;

    // robot state
    double tilt = 0, tinc = degree, swivel=0, angle=0, height=0, hinc=0.005;

    // message declarations
    geometry_msgs::TransformStamped odom_trans;
    //sensor_msgs::JointState joint_state;
    odom_trans.header.frame_id = "odom";
    odom_trans.child_frame_id = "base_link";

    while (ros::ok()) {
        //update joint_state
        //joint_state.header.stamp = ros::Time::now();


        // update transform
        // (moving in a circle with radius=2)
        odom_trans.header.stamp = ros::Time::now();
        odom_trans.transform.translation.x = 2;
        odom_trans.transform.translation.y = 2;
        odom_trans.transform.translation.z = 0.2;
        odom_trans.transform.rotation = tf::createQuaternionMsgFromYaw(30.0);

        //send the joint state and transform
        //joint_pub.publish(joint_state);
        broadcaster.sendTransform(odom_trans);


        // This will adjust as needed per iteration
        loop_rate.sleep();
    }


    return 0;
}
