#include <ros/ros.h>
#include <armadillo>
#include "pf.h"
#include <geometry_msgs/Twist.h>
#include <ras_line_detector/LineSegmentList.h>
#include <ras_line_detector/LineSegment.h>
#include <ras_particle_filter/Particles.h>

#include <iostream>

class ParticleFilterNode
{
public:
    ros::NodeHandle nh;
    ros::Subscriber twist_sub, lines_sub;
    ros::Publisher particles_pub, pose_estimate_pub;
    ParticleFilter pf;

    geometry_msgs::Twist twist;
    ras_line_detector::LineSegmentList line_segments;

    double Lambda_psi;
    std::vector<double> start_pose;
    std::vector<double> bound = {0.0, 1.0, 0.0, 1.0};
    std::vector<double> R_vec = {0.01, 0.01, 0.01};
    std::vector<double> Q_vec = {0.01, 0.01};
    double part_bound;
    int M;
    int update_freq, predict_freq;
    int npp;
    string map_file;
    double init_particle_spread;

    mat S, R, Q, S_bar, z, W;
    cube Psi;
    rowvec outlier;

    void TwistCallback(const geometry_msgs::Twist::ConstPtr& msg)
    {
        twist = *msg;
    }

    void LinesCallback(const ras_line_detector::LineSegmentList::ConstPtr& msg)
    {
        line_segments = *msg;
        int n;
        if((n = line_segments.line_segments.size()) > 0){
            z = mat(2, n);
            for(int i = 0; i < n; i++){
                z.col(i) = vec({line_segments.line_segments[n].radius,
                               line_segments.line_segments[n].angle*datum::pi/180.0});
            }
        }
    }

    ParticleFilterNode()
    {
        lines_sub = nh.subscribe("/lines", 1, &ParticleFilterNode::LinesCallback, this);
        twist_sub = nh.subscribe("/localization/odometry_twist", 1, &ParticleFilterNode::TwistCallback, this);
        particles_pub = nh.advertise<ras_particle_filter::Particles>("/particles", 1);
        pose_estimate_pub = nh.advertise<geometry_msgs::Pose2D>("/localization/pose", 1);
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
    }

    void InitializePf()
    {
        R = diagmat(vec(R_vec));
        Q = diagmat(vec(Q_vec));
        pf.init(vec(bound), part_bound, vec(start_pose), S, M, init_particle_spread);
    }

    void PublishParticles()
    {
        ras_particle_filter::Particles particles_msg;
        for(int i = 0; i < npp; i++){
            geometry_msgs::Pose2D pose;
            pose.x = S(0, i);
            pose.y = S(1, i);
            pose.theta = S(2, i);
            particles_msg.poses.push_back(pose);
        }
        particles_pub.publish(particles_msg);
    }

    void PublishPoseEstimate()
    {
        mat s = S.rows(0, 2).cols(0, npp-1);
        vec m = mean(s, 1);
        geometry_msgs::Pose2D pose;
        pose.x = m(0);
        pose.y = m(1);
        pose.theta = m(2);
        pose_estimate_pub.publish(pose);
    }

    void MCL(long& ct, double& t, int freq_ratio)
    {
        double dt = ros::Time::now().toSec() - t;
        double v = twist.linear.x;
        double w = twist.angular.z;
        S_bar = pf.predict(S, v, w, R, dt);
        t = ros::Time::now().toSec();
        if(ct++%freq_ratio==0)
        {
            W = readLines(map_file);
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

        arma::mat m(2,2);
        m(0,0) = 1;
        m(0,1) = 2;
        m(1,0) = 3;
        m(1,1) = 4;
        arma::vec v(2);
        v(0) = 10; v(1) = 10;
        arma::vec mv = m*v;
        ROS_INFO("(%f, %f)", mv(0), mv(1));

        ros::Rate rate(pfn.predict_freq);
        int freq_ratio = pfn.predict_freq / pfn.update_freq;
        long ct = 0;
        double t = ros::Time::now().toSec()-0.1;
        while(ros::ok())
        {
            ros::spinOnce();

            pfn.MCL(ct, t, freq_ratio);

            pfn.PublishParticles();
            pfn.PublishPoseEstimate();

            rate.sleep();
        }
}
