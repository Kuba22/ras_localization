#include <ros/ros.h>
#include <armadillo>
#include "pf.h"
#include <geometry_msgs/Twist.h>
#include <sensor_msgs/LaserScan.h>
#include <geometry_msgs/PoseArray.h>
#include <geometry_msgs/PoseStamped.h>
#include <tf/transform_broadcaster.h>

#include <iostream>

bool received_scan = false;

class ParticleFilterNode
{
public:
    ros::NodeHandle nh;
    ros::Subscriber twist_sub, scan_sub;
    ros::Publisher particles_pub, pose_estimate_pub;
	tf::TransformBroadcaster broadcaster;
    ParticleFilter pf;

    geometry_msgs::Twist twist;
    sensor_msgs::LaserScan scan;

    double Lambda_psi;
    std::vector<double> start_pose;
    std::vector<double> bound = {0.0, 1.0, 0.0, 1.0};
    std::vector<double> R_vec = {0.01, 0.01, 0.01};
    std::vector<double> Q_vec = {0.01, 0.01};
    double part_bound;
    int M;
    int update_freq, predict_freq;
    int npp, n_phi;
    string map_file;
    double init_particle_spread;

    mat S, R, Q, S_bar, z, W;
    cube Psi;
    rowvec outlier;

    void TwistCallback(const geometry_msgs::Twist::ConstPtr& msg)
    {
        twist = *msg;
    }

    void ScanCallback(const sensor_msgs::LaserScan::ConstPtr& msg)
    {
        scan = *msg;
        received_scan = true;
    }

    ParticleFilterNode()
    {
        scan_sub = nh.subscribe("/scan", 1, &ParticleFilterNode::ScanCallback, this);
        twist_sub = nh.subscribe("/localization/odometry_twist", 1, &ParticleFilterNode::TwistCallback, this);
        particles_pub = nh.advertise<geometry_msgs::PoseArray>("/particles", 1);
        pose_estimate_pub = nh.advertise<geometry_msgs::PoseStamped>("/localization/pose", 1);
        LoadParams();
        InitializePf();
    }

   void LoadParams()
    {
        nh.param("update_freq", update_freq, 5);
        nh.param("predict_freq", predict_freq, 10);
        nh.param("Lambda_psi", Lambda_psi, 0.0001);
        nh.getParam("start_pose", start_pose);
        nh.getParam("bound", bound);
        nh.param("part_bound", part_bound, 20.0);
        nh.param("M", M, 1000);
        nh.getParam("R", R_vec);
        nh.getParam("Q", Q_vec);
        nh.param("n_particles_to_pub", npp, 20);
        npp = npp > M ? M : npp;
        nh.getParam("map_file", map_file);
        nh.param("init_particle_spread", init_particle_spread, -1.0);
        nh.param("n_phi", n_phi, 10);
    }

    void InitializePf()
    {
        R = diagmat(vec(R_vec));
        Q = diagmat(vec(Q_vec));
        pf.init(vec(bound), part_bound, vec(start_pose), S, R, Q, M, init_particle_spread);
        W = PfHelper::readLines(map_file);
    }

    void PublishParticles()
    {
        geometry_msgs::PoseArray particles_msg;
        particles_msg.header.frame_id = "odom";
        for(int i = 0; i < npp; i++){
            geometry_msgs::Pose pose;
            pose.position.x = S(0, i);
            pose.position.y = S(1, i);
            pose.orientation.z = sin(0.5*S(2, i));
            pose.orientation.w = cos(0.5*S(2, i));
            particles_msg.poses.push_back(pose);
        }
        particles_pub.publish(particles_msg);
    }

    void PublishPoseEstimate()
    {
        mat s = S.rows(0, 2).cols(0, npp-1);
        vec m = mean(s, 1);
        geometry_msgs::PoseStamped pose;
        pose.header.frame_id = "odom";
        pose.header.stamp = ros::Time::now();
        pose.pose.position.x = m(0);
        pose.pose.position.y = m(1);
        pose.pose.orientation.z = sin(0.5*m(2));
        pose.pose.orientation.w = cos(0.5*m(2));
        pose_estimate_pub.publish(pose);

		geometry_msgs::TransformStamped odom_trans;
	    odom_trans.header.frame_id = "odom";
	    odom_trans.child_frame_id = "base_link";
		odom_trans.header.stamp = ros::Time::now();
	    odom_trans.transform.translation.x = m(0);
	    odom_trans.transform.translation.y = m(1);
	    odom_trans.transform.translation.z = 0;
	    odom_trans.transform.rotation = tf::createQuaternionMsgFromYaw(m(2));
	    broadcaster.sendTransform(odom_trans);
    }

    void ObservationsFromScan(){
        int ntheta = std::round((scan.angle_max - scan.angle_min)/scan.angle_increment+1);
        int d = ntheta/n_phi;
        z = mat(2, 1);
        for(int i=0; i<n_phi; i++){
            if(std::isinf(scan.ranges[i*d]))
                continue;
            z = join_rows(z, vec({scan.ranges[i*d], i*d*scan.angle_increment}));
        }
        z.shed_col(0);
    }

    void MCL(long& ct, double& t, int freq_ratio)
    {
        double dt = ros::Time::now().toSec() - t;
        double v = twist.linear.x;
        double w = twist.angular.z;
        S_bar = pf.predict(S, v, w, dt);
        t = ros::Time::now().toSec();
        if(ct++%freq_ratio==0)
        {
            if(!received_scan){
                ROS_WARN("No scan messages received yet.");
                S = S_bar;
                return;
            }
            ObservationsFromScan();
            pf.associate(S_bar, z, W, Lambda_psi, Q, outlier, Psi);
            int outliers = (int)double(arma::as_scalar(arma::sum(outlier)));
            if (outliers == outlier.n_elem) {
                S = S_bar;
                return;
            }
            S_bar = pf.weight(S_bar, Psi, outlier);
            S = pf.resample(S_bar);
        }
        else{
            S = S_bar;
        }
    }
};

int main(int argc, char *argv[])
{
	ros::init(argc, argv, "particle_filter");
        ROS_INFO("Particle Filter Node");
        ParticleFilterNode pfn;
		
        ros::Rate rate(pfn.predict_freq);
        int freq_ratio = (double)pfn.predict_freq / (double)pfn.update_freq;
        long ct = 0;
        double t = ros::Time::now().toSec();
        while(ros::ok())
        {
            ros::spinOnce();

            pfn.MCL(ct, t, freq_ratio);

            pfn.PublishParticles();
            pfn.PublishPoseEstimate();

            rate.sleep();
        }
}
