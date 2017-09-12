#!/usr/bin/env python

# Copyright (C) 2017 Electric Movement Inc.
#
# This file is part of Robotic Arm: Pick and Place project for Udacity
# Robotics nano-degree program
#
# All Rights Reserved.

# Author: Harsh Pandya

# import modules
import rospy
import tf
from kuka_arm.srv import *
from trajectory_msgs.msg import JointTrajectory, JointTrajectoryPoint
from geometry_msgs.msg import Pose
from mpmath import *
from sympy import *
from math import copysign


def Get_TFMartix(d, a, q, alpha):
    TF_Matrix = Matrix([
        [ cos(q),               -sin(q),            0,              a    ],
        [sin(q)*cos(alpha), cos(q)*cos(alpha), -sin(alpha), -sin(alpha)*d],
        [sin(q)*sin(alpha), cos(q)*sin(alpha),  cos(alpha),  cos(alpha)*d],
        [       0,                  0,              0,              1    ]
        ])
    return TF_Matrix


def handle_calculate_IK(req):
    rospy.loginfo("Received %s eef-poses from the plan" % len(req.poses))
    if len(req.poses) < 1:
        print "No valid poses received"
        return -1
    else:
		
        ### Your FK code here
        # Create symbols

        d1, d2, d3, d4, d5, d6, d7 = symbols('d1:8')
        a0, a1, a2, a3, a4, a5, a6 = symbols('a0:7')
        q1, q2, q3, q4, q5, q6, q7 = symbols('q1:8')
        alp0, alp1, alp2, alp3, alp4, alp5, alp6 = symbols('alp0:7')
	#
	#   
	# Create Modified DH parameters

        DH_table = { 
                alp0:       0,  a0:      0,  d1:  0.75, q1:  q1,
                alp1: -pi/2.0,  a1:   0.35,  d2:     0, q2:  -pi/2.0 + q2,
                alp2:       0,  a2:   1.25,  d3:     0, q3:  q3,
                alp3: -pi/2.0,  a3: -0.054,  d4:   1.5, q4:  q4,
                alp4:  pi/2.0,  a4:      0,  d5:     0, q5:  q5,
                alp5: -pi/2.0,  a5:      0,  d6:     0, q6:  q6,
                alp6: 0,        a6:      0,  d7: 0.303, q7:  0
                }
	#
	#            
	# Define Modified DH Transformation matrix

        T0_1 = Get_TFMartix(d1, a0, q1, alp0).subs(DH_table)
        T1_2 = Get_TFMartix(d2, a1, q2, alp1).subs(DH_table)
        T2_3 = Get_TFMartix(d3, a2, q3, alp2).subs(DH_table)
        T3_4 = Get_TFMartix(d4, a3, q4, alp3).subs(DH_table)
        T4_5 = Get_TFMartix(d5, a4, q5, alp4).subs(DH_table)
        T5_6 = Get_TFMartix(d6, a5, q6, alp5).subs(DH_table)
        T6_EE = Get_TFMartix(d7, a6, q7, alp6).subs(DH_table)

        T0_EE = T0_1 * T1_2 * T2_3 * T3_4 * T4_5 * T5_6 * T6_EE
	#
	#
	# Create individual transformation matrices
	#
	#
	# Extract rotation matrices from the transformation matrices

        r, p, y = symbols('r p y')

        rot_x = Matrix([[1,     0,      0  ],
                        [0, cos(r), -sin(r)],
                        [0, sin(r), cos(r)]])

        rot_y = Matrix([[cos(p), 0 ,sin(p)],
                        [0,      1 ,    0  ],
                        [-sin(p),0 ,cos(p)]])

        rot_z = Matrix([[cos(y), -sin(y), 0],
                        [sin(y), cos(y),  0],
                        [0,      0,       1]])

        rot_err = rot_z.subs(y, radians(180)) * rot_y.subs(p, radians(-90))

        rot_ee = rot_z * rot_y * rot_x * rot_err

	#
	#
        ###

        # Initialize service response
        joint_trajectory_list = []
        for x in xrange(0, len(req.poses)):
            # IK code starts here
            joint_trajectory_point = JointTrajectoryPoint()

	    # Extract end-effector position and orientation from request
	    # px,py,pz = end-effector position
	    # roll, pitch, yaw = end-effector orientation
            px = req.poses[x].position.x
            py = req.poses[x].position.y
            pz = req.poses[x].position.z

            (roll, pitch, yaw) = tf.transformations.euler_from_quaternion(
                [req.poses[x].orientation.x, req.poses[x].orientation.y,
                    req.poses[x].orientation.z, req.poses[x].orientation.w])
     
            ### Your IK code here 
	    # Compensate for rotation discrepancy between DH parameters and Gazebo
	    #
            rot_ee = rot_ee.subs({'r':roll, 'p':pitch, 'y':yaw})
	    #
	    # Calculate joint angles using Geometric IK method
	    #

            ee_vector = Matrix([[px], [py], [pz]])

            #'Calculating coordinates of Wrist Center'

            wc = ee_vector - (0.303) * rot_ee[:,2]

            #calculate geometry for first 3 joints

            #print('Calculating geometry for first 3 joints')

            #triangle sides

            A = 1.501 #from URDF file
            C = 1.25 #from DH table
            #B = sqrt(pow((sqrt(wc[0]*wc[0] + wc[1]*wc[1]) - DH_table[a1]), 2) + pow(wc[2] - DH_table[d1],2))
            B = sqrt(pow((sqrt(wc[0]*wc[0] + wc[1]*wc[1]) - DH_table[a1]), 2) + pow((wc[2] - DH_table[d1]),2))
            #print('Calculating triangle angles')

            #triangle angles
            a = acos((B*B + C*C - A*A) / (2 * B * C))
            b = acos((A*A + C*C - B*B) / (2 * A * C))
            c = acos((A*A + B*B - C*C) / (2 * A * B))

            # to compensate 55mm offset in link4
            #offset_angle = atan2(0.054, 1.501)
            offset_angle = atan2(abs(DH_table[a3]), A)

            #print(offset_angle)

            #print('Calculating first 3 thetas')

            theta1 = atan2(wc[1], wc[0])
            theta2 = pi/2 - a - atan2(wc[2] - DH_table[d1], sqrt(wc[0]*wc[0] + wc[1]*wc[1]) - DH_table[a1]) % (2*pi)
            theta3 = pi/2 - (b + offset_angle) % (2*pi)

            #print('Calculating rotation matrix for joints 1-3')

            R0_3 = T0_1[0:3, 0:3] * T1_2[0:3, 0:3] * T2_3[0:3, 0:3] 
            R0_3 = R0_3.evalf(subs={q1: theta1, q2: theta2, q3: theta3})


            #print('Inverting R0_3 matrix')

            #R0_3_inv = R0_3.inv("LU")

            #print('Calculating rotation matrix for joints 3-6')

            R3_6 = R0_3.inv("LU") * rot_ee

            #print('Calculating thetas 3-6')

            theta4 = atan2(R3_6[2,2], -R3_6[0,2]) % ((2*pi).evalf())
            theta5 = atan2(sqrt(R3_6[0,2] * R3_6[0,2] + R3_6[2,2] * R3_6[2,2]), R3_6[1,2]) % ((2*pi).evalf())
            theta6 = atan2(-R3_6[1,1], R3_6[1,0]) % ((2*pi).evalf)
	    #
            ###
		
            # Populate response for the IK request
            # In the next line replace theta1,theta2...,theta6 by your joint angle variables
            #print([theta1, theta2, theta3, theta4, theta5, theta6])
	    joint_trajectory_point.positions = [theta1, theta2, theta3, theta4, theta5, theta6]
	    joint_trajectory_list.append(joint_trajectory_point)

        rospy.loginfo("length of Joint Trajectory List: %s" % len(joint_trajectory_list))
        return CalculateIKResponse(joint_trajectory_list)


def IK_server():
    # initialize node and declare calculate_ik service
    rospy.init_node('IK_server')
    s = rospy.Service('calculate_ik', CalculateIK, handle_calculate_IK)
    print "Ready to receive an IK request"
    rospy.spin()

if __name__ == "__main__":
    IK_server()
