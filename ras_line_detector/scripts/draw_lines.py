#!/usr/bin/env python
import rospy
from ras_line_detector.msg import LineSegmentList
from ras_line_detector.msg import LineSegment
import pygame
from math import atan2, sqrt, cos, sin, asin
from geometry_msgs.msg import PoseStamped

pi = 3.14159

scale = 100
dx = 200
dy = 600 - 200
	
def drawLine(screen, line, p):
	(rho, theta) = line
	x = rho*cos(theta) + p[0]
	y = rho*sin(theta) + p[1]
	p1 = [x, y]
	pygame.draw.lines(screen,(255,0,0), False, [(p[0]*scale+dx,-p[1]*scale+dy),(p1[0]*scale+dx,-p1[1]*scale+dy)],1)
	pygame.display.update()

class Drawer:
	def __init__(self):
		rospy.init_node("drawer", anonymous = True)
		self.scan_sub = rospy.Subscriber("/lines", LineSegmentList, self.linesCallback)
		self.scan_sub = rospy.Subscriber("/localization/odometry_pose", PoseStamped, self.poseCallback)
		pygame.init()
		self.screen = pygame.display.set_mode((600,600))

	def poseCallback(self, msg):
		self.pose = msg.pose
		self.x = msg.pose.position.x
		self.y = msg.pose.position.y
		self.theta = 2*asin(msg.pose.orientation.z)

	def linesCallback(self, msg):
		self.screen.fill((255,255,255))
		lines = msg.line_segments
		for i in range(len(lines)):
			rho = lines[i].radius
			theta = lines[i].angle*pi/180 + self.theta
			drawLine(self.screen, (rho, theta), (self.x, self.y))
		pygame.draw.circle(self.screen, (0, 0, 255), (int(self.x*scale)+dx, -int(self.y*scale)+dy), 4, 0)
		pygame.display.update()

	
if __name__=='__main__':
	try:
		node = Drawer()
		rate = rospy.Rate(10)
		while not rospy.is_shutdown():
			rate.sleep()
			for event in pygame.event.get():
				if event.type == pygame.QUIT:
					pygame.quit()
	except rospy.ROSInterruptException:
		pass	

