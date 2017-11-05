#!/usr/bin/env python

import rospy
from sensor_msgs.msg import LaserScan
from ras_line_detector.msg import LineSegmentList
from ras_line_detector.msg import LineSegment
import numpy as np
from scipy.signal import *

def getPrevious(a, n):
	if np.isinf(a[n-1]):
		return getPrevious(a, n-1)
	else:
		return a[n-1]

def	moveUntilDifferent(a, n, direction):
	# dir +-1
	curr_val = a[n]	
	if n + direction < 0:
		n += len(a)
	if n + direction >= len(a):
		n -= len(a)
	if a[n+direction] == curr_val:
		return moveUntilDifferent(a, n+direction, direction)
	else:
		return n+direction
	
def getMin(a):
	mins = []
	for i in range(len(a)):
		j = moveUntilDifferent(a, i, -1)
		k = moveUntilDifferent(a, i, 1)
		if a[k] > a[i] and a[j] > a[i]:
			if k > j:
				m = int(round((k + j) * 0.5))
				mins.append(m if not m < 0 else m + len(a))
			else:
				m = int(round((k + j + len(a)) * 0.5))
				mins.append(m if not m >= len(a) else m - len(a))
	return set(mins)


class LineDetector:

	def __init__(self):
		self.n = 0
		self.ranges = list()
		self.bearings = list()
		rospy.init_node("line_detector", anonymous = True)
		self.scan_sub = rospy.Subscriber("/scan", LaserScan, self.scanCallback)
		self.line_pub = rospy.Publisher("/lines", LineSegmentList, queue_size = 1)
		self.scan_msg = None

	def publish(self):
		if self.scan_msg is None:
			return
		scan = self.scan_msg
		ranges = self.ranges
		bearings = self.bearings
		noinf = np.array([None]*len(ranges))
		for i in range(len(ranges)):
			noinf[i] = ranges[i] if not np.isinf(ranges[i]) else getPrevious(ranges, i)
		N = 5
		noinf_Nlonger = np.append(noinf, noinf[:N])
		noinf_Nlonger = np.insert(noinf_Nlonger, 0, noinf[-N:])
		gaussian_smoothed = np.convolve(general_gaussian(2*N+1, p=1, sig=16), noinf_Nlonger, mode='valid')
		grad = np.gradient(noinf)
		incontinuities = np.where(np.absolute(grad) > 0.1)
		minima = getMin(gaussian_smoothed)
		lines = list()
		for m in minima:
			line = LineSegment()
			line.radius = noinf[m]
			line.angle = bearings[m]
			lines.append(line)
		msg = LineSegmentList()
		msg.line_segments = lines
		self.line_pub.publish(msg)

	def scanCallback(self, msg):
		self.scan_msg = msg
		self.ranges = msg.ranges
		self.bearings = [i for i in range(int(round((self.scan_msg.angle_max-self.scan_msg.angle_min)/self.scan_msg.angle_increment))+1)]


if __name__=='__main__':
	try:
		node = LineDetector()
		rate = rospy.Rate(10)
		while not rospy.is_shutdown():
			node.publish()
			rate.sleep()
	except rospy.ROSInterruptException:
		pass
