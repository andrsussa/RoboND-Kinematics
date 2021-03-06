from sympy import *
from time import time
from mpmath import radians
import tf
import numpy as np

'''
Format of test case is [ [[EE position],[EE orientation as quaternions]],[WC location],[joint angles]]
You can generate additional test cases by setting up your kuka project and running `$ roslaunch kuka_arm forward_kinematics.launch`
From here you can adjust the joint angles to find thetas, use the gripper to extract positions and orientation (in quaternion xyzw) and lastly use link 5
to find the position of the wrist center. These newly generated test cases can be added to the test_cases dictionary.
'''

test_cases = {1:[[[2.16135,-1.42635,1.55109],
                  [0.708611,0.186356,-0.157931,0.661967]],
                  [1.89451,-1.44302,1.69366],
                  [-0.65,0.45,-0.36,0.95,0.79,0.49]],
              2:[[[-0.56754,0.93663,3.0038],
                  [0.62073, 0.48318,0.38759,0.480629]],
                  [-0.638,0.64198,2.9988],
                  [-0.79,-0.11,-2.33,1.94,1.14,-3.68]],
              3:[[[-1.3863,0.02074,0.90986],
                  [0.01735,-0.2179,0.9025,0.371016]],
                  [-1.1669,-0.17989,0.85137],
                  [-2.99,-0.12,0.94,4.06,1.29,-4.12]],
              4:[],
              5:[]}


def test_code(test_case):
    ## Set up code
    ## Do not modify!
    x = 0
    class Position:
        def __init__(self,EE_pos):
            self.x = EE_pos[0]
            self.y = EE_pos[1]
            self.z = EE_pos[2]
    class Orientation:
        def __init__(self,EE_ori):
            self.x = EE_ori[0]
            self.y = EE_ori[1]
            self.z = EE_ori[2]
            self.w = EE_ori[3]

    position = Position(test_case[0][0])
    orientation = Orientation(test_case[0][1])

    class Combine:
        def __init__(self,position,orientation):
            self.position = position
            self.orientation = orientation

    comb = Combine(position,orientation)

    class Pose:
        def __init__(self,comb):
            self.poses = [comb]

    req = Pose(comb)
    start_time = time()
    
    ########################################################################################
    ## 

    print('Defining symbols')

    ## Insert IK code here!
    d1, d2, d3, d4, d5, d6, d7 = symbols('d1:8')
    a0, a1, a2, a3, a4, a5, a6 = symbols('a0:7')
    q1, q2, q3, q4, q5, q6, q7 = symbols('q1:8')
    alp0, alp1, alp2, alp3, alp4, alp5, alp6 = symbols('alp0:7')

    print('Reading positions')
    ##End effector coordinates
    px = req.poses[x].position.x
    py = req.poses[x].position.y
    pz = req.poses[x].position.z

    ## Get RPY angles of end effector
    print('Getting RPY angles')

    (roll, pitch, yaw) = tf.transformations.euler_from_quaternion([
        req.poses[x].orientation.x,req.poses[x].orientation.y,req.poses[x].orientation.z,
        req.poses[x].orientation.w])

    #DH table
    print('Defining DH table')

    DH_table = { 
            alp0:       0,  a0:      0,  d1:  0.75, q1:  q1,
            alp1: -pi/2.0,  a1:   0.35,  d2:     0, q2:  -pi/2.0 + q2,
            alp2:       0,  a2:   1.25,  d3:     0, q3:  q3,
            alp3: -pi/2.0,  a3: -0.054,  d4:   1.5, q4:  q4,
            alp4:  pi/2.0,  a4:      0,  d5:     0, q5:  q5,
            alp5: -pi/2.0,  a5:      0,  d6:     0, q6:  q6,
            alp6: 0,        a6:      0,  d7: 0.303, q7:  0
            }

    print('Defining transformation matrices')

    T0_1 = Get_TFMartix(d1, a0, q1, alp0).subs(DH_table)
    T1_2 = Get_TFMartix(d2, a1, q2, alp1).subs(DH_table)
    T2_3 = Get_TFMartix(d3, a2, q3, alp2).subs(DH_table)
    T3_4 = Get_TFMartix(d4, a3, q4, alp3).subs(DH_table)
    T4_5 = Get_TFMartix(d5, a4, q5, alp4).subs(DH_table)
    T5_6 = Get_TFMartix(d6, a5, q6, alp5).subs(DH_table)
    T6_EE = Get_TFMartix(d7, a6, q7, alp6).subs(DH_table)

    print('Calculating total transformation matrix')

    T0_EE = T0_1 * T1_2 * T2_3 * T3_4 * T4_5 * T5_6 * T6_EE


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

    print('Calculating total rotation matrix')

    rot_ee = rot_z * rot_y * rot_x * rot_err

    rot_ee = rot_ee.subs({'r':roll, 'p':pitch, 'y':yaw})
    ee_vector = Matrix([[px], [py], [pz]])

    print('Calculating coordinates of Wrist Center')

    wc = ee_vector - (DH_table[d7]) * rot_ee[:,2]

    #calculate geometry for first 3 joints

    print('Calculating geometry for first 3 joints')

    #triangle sides

    A = 1.501 #from URDF file
    C = 1.25 #from DH table
    B = sqrt(pow((sqrt(wc[0]*wc[0] + wc[1]*wc[1]) - DH_table[a1]), 2) + pow((wc[2] - DH_table[d1]),2))
    print('Calculating triangle ABC angles')

    #triangle angles
    a = acos((B*B + C*C - A*A) / (2 * B * C))
    b = acos((A*A + C*C - B*B) / (2 * A * C))
    c = acos((A*A + B*B - C*C) / (2 * A * B))

    # to compensate 55mm offset in link4
    offset_angle = atan2(abs(DH_table[a3]), A)

    print(offset_angle)

    print('Calculating thetas 1-3')

    theta1 = atan2(wc[1], wc[0]).evalf()
    theta2 = pi/2 - a - atan2(wc[2] - DH_table[d1], sqrt(wc[0]*wc[0] + wc[1]*wc[1]) - DH_table[a1])
    theta3 = pi/2 - (b + offset_angle)

    print('Calculating rotation matrix for joints 1-3')

    R0_3 = T0_1[0:3, 0:3] * T1_2[0:3, 0:3] * T2_3[0:3, 0:3] 
    R0_3 = R0_3.evalf(subs={q1: theta1, q2: theta2, q3: theta3})


    R3_6 = R0_3.inv("LU") * rot_ee

    R3_6 = R3_6.evalf()



    print('Calculating thetas 4-6')
    theta4 = (atan2(R3_6[2,2],-R3_6[0,2])
    theta5 = (atan2(sqrt((R3_6[0,2])**2 + (R3_6[2,2])**2), R3_6[1,2])).evalf()
    theta6 = (atan2(-R3_6[1,1],R3_6[1,0])).evalf()


    ## 
    ########################################################################################
    
    ########################################################################################
    ## For additional debugging add your forward kinematics here. Use your previously calculated thetas
    ## as the input and output the position of your end effector as your_ee = [x,y,z]

    ## (OPTIONAL) YOUR CODE HERE!
    FK = T0_EE.evalf(subs={q1: theta1, q2: theta2, q3: theta3, q4: theta4, q5: theta5, q6: theta6})

    ## End your code input for forward kinematics here!
    ########################################################################################

    ## For error analysis please set the following variables of your WC location and EE location in the format of [x,y,z]
    your_wc = [wc[0],wc[1],wc[2]] # <--- Load your calculated WC values in this array
    your_ee = [FK[0,3],FK[1,3],FK[2,3]] # <--- Load your calculated end effector value from your forward kinematics
    ########################################################################################

    ## Error analysis
    print ("\nTotal run time to calculate joint angles from pose is %04.4f seconds" % (time()-start_time))

    # Find WC error
    if not(sum(your_wc)==3):
        wc_x_e = abs(your_wc[0]-test_case[1][0])
        wc_y_e = abs(your_wc[1]-test_case[1][1])
        wc_z_e = abs(your_wc[2]-test_case[1][2])
        wc_offset = sqrt(wc_x_e**2 + wc_y_e**2 + wc_z_e**2)
        print ("\nWrist error for x position is: %04.8f" % wc_x_e)
        print ("Wrist error for y position is: %04.8f" % wc_y_e)
        print ("Wrist error for z position is: %04.8f" % wc_z_e)
        print ("Overall wrist offset is: %04.8f units" % wc_offset)

    # Find theta errors
    t_1_e = abs(theta1-test_case[2][0])
    t_2_e = abs(theta2-test_case[2][1])
    t_3_e = abs(theta3-test_case[2][2])
    t_4_e = abs(theta4-test_case[2][3])
    t_5_e = abs(theta5-test_case[2][4])
    t_6_e = abs(theta6-test_case[2][5])
    print ("\nTheta 1 error is: %04.8f" % t_1_e)
    print ("Theta 2 error is: %04.8f" % t_2_e)
    print ("Theta 3 error is: %04.8f" % t_3_e)
    print ("Theta 4 error is: %04.8f" % t_4_e)
    print ("Theta 5 error is: %04.8f" % t_5_e)
    print ("Theta 6 error is: %04.8f" % t_6_e)
    print ("\n**These theta errors may not be a correct representation of your code, due to the fact \
           \nthat the arm can have muliple positions. It is best to add your forward kinmeatics to \
           \nconfirm whether your code is working or not**")
    print (" ")

    # Find FK EE error
    if not(sum(your_ee)==3):
        ee_x_e = abs(your_ee[0]-test_case[0][0][0])
        ee_y_e = abs(your_ee[1]-test_case[0][0][1])
        ee_z_e = abs(your_ee[2]-test_case[0][0][2])
        ee_offset = sqrt(ee_x_e**2 + ee_y_e**2 + ee_z_e**2)
        print ("\nEnd effector error for x position is: %04.8f" % ee_x_e)
        print ("End effector error for y position is: %04.8f" % ee_y_e)
        print ("End effector error for z position is: %04.8f" % ee_z_e)
        print ("Overall end effector offset is: %04.8f units \n" % ee_offset)


def Get_TFMartix(d, a, q, alpha):
    TF_Matrix = Matrix([
        [ cos(q),               -sin(q),            0,              a    ],
        [sin(q)*cos(alpha), cos(q)*cos(alpha), -sin(alpha), -sin(alpha)*d],
        [sin(q)*sin(alpha), cos(q)*sin(alpha),  cos(alpha),  cos(alpha)*d],
        [       0,                  0,              0,              1    ]
        ])
    return TF_Matrix


if __name__ == "__main__":
    # Change test case number for different scenarios
    test_case_number = 2

    test_code(test_cases[test_case_number])
