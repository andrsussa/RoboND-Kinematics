## Project: Kinematics Pick & Place
### Writeup Template: You can use this file as a template for your writeup if you want to submit it as a markdown file, but feel free to use some other method and submit a pdf if you prefer.

---


**Steps to complete the project:**  


1. Set up your ROS Workspace.
2. Download or clone the [project repository](https://github.com/udacity/RoboND-Kinematics-Project) into the ***src*** directory of your ROS Workspace.  
3. Experiment with the forward_kinematics environment and get familiar with the robot.
4. Launch in [demo mode](https://classroom.udacity.com/nanodegrees/nd209/parts/7b2fd2d7-e181-401e-977a-6158c77bf816/modules/8855de3f-2897-46c3-a805-628b5ecf045b/lessons/91d017b1-4493-4522-ad52-04a74a01094c/concepts/ae64bb91-e8c4-44c9-adbe-798e8f688193).
5. Perform Kinematic Analysis for the robot following the [project rubric](https://review.udacity.com/#!/rubrics/972/view).
6. Fill in the `IK_server.py` with your Inverse Kinematics code. 


[//]: # (Image References)

[image1]: ./misc_images/misc1.png
[image2]: ./misc_images/misc3.png
[image3]: ./misc_images/misc2.png

[DH_parameters]: ./misc_images/DH_parameters.jpg
[triangle_sides]: ./misc_images/triangle_sides.jpg
[thetas1-3]: ./misc_images/thetas1-3.jpg
[thetas4-6]: ./misc_images/thetas4-6.jpg
[transform_matrices]: ./misc_images/transform_matrices.jpg
[rotation_matrices]: ./misc_images/rotation_matrices.jpg
[correction_matrix]: ./misc_images/correction_matrix.jpg
[kuka_arm0]: ./misc_images/kuka_arm0.jpg
[kuka_arm1]: ./misc_images/kuka_arm.jpg
[kuka_arm_fk]: ./misc_images/kuka_arm_fk.jpg

## [Rubric](https://review.udacity.com/#!/rubrics/972/view) Points
### Here I will consider the rubric points individually and describe how I addressed each point in my implementation.  

---
### Writeup / README

#### 1. Provide a Writeup / README that includes all the rubric points and how you addressed each one.  You can submit your writeup as markdown or pdf.  

You're reading it!

### Kinematic Analysis
#### 1. Run the forward_kinematics demo and evaluate the kr210.urdf.xacro file to perform kinematic analysis of Kuka KR210 robot and derive its DH parameters.

Run roslaunch kuka_arm forward_kinematics.launch

![Kuka arm forward kinematics][kuka_arm_fk]

#### DH parameters calculation


![DH parameters][DH_parameters]

DH-Table

i | alpha_i-1 | a_i-1 | d_i | theta_i
--- | --- | --- | --- | ---
T0_1 | 0 | 0 | 0.75 | theta_1
T1_2 | -90° | 0.35 | 0 | theta_2-90°
T2_3 | 0 | 1.25 | 0 | theta_3
T3_4 | -90° | -0.054 | 1.5 | theta_4
T4_5 | 90° | 0 | 0 | theta_5
T5_6 | -90° | 0 | 0 | theta_6
T6_G | 0 | 0 | 0.303 | 0


#### 2. Using the DH parameter table you derived earlier, create individual transformation matrices about each joint. In addition, also generate a generalized homogeneous transform between base_link and gripper_link using only end-effector(gripper) pose.

##### 2.1 Individual transformation matrices

![Transform matrices][transform_matrices]

##### 2.1 (transformation matrices continued) + Rotation correction matrix

![Correction matrices][correction_matrix]


##### 2.2 Calculating rotation matrices and wrist center position

![Rotation matrices][rotation_matrices]



#### 3. Decouple Inverse Kinematics problem into Inverse Position Kinematics and inverse Orientation Kinematics; doing so derive the equations to calculate all individual joint angles.
 

##### 3.1 Calculating triangle sides

![Calculating ABC triangle sides][triangle_sides]

##### 3.2 Calculating triangle angles and Thetas 1-3

![Calculating theta 1-3][thetas1-3]

##### 3.3 Calculating triangle angles and Thetas 4-6

![Calculcating theta 4-6][thetas4-6]

### Project Implementation

#### 1. Fill in the `IK_server.py` file with properly commented python code for calculating Inverse Kinematics based on previously performed Kinematic Analysis. Your code must guide the robot to successfully complete 8/10 pick and place cycles. Briefly discuss the code you implemented and your results. 


*  After all calculations have been done on paper, programming part was actually pretty straight forward:
1. Implemented IK calculations in IK_debug.py script, and it worked perfetly only for 1-st test case. As per suggestion I had to implement FK calculations as well in order to check previously calculated IK solution.
2. Moved code to IK_server.py and run simulator. Although it was seemed to be moving along the suggested trajectory, values of thetas3-6 apparently were way off optimal - kuka was swinging its wrist and even kicked out target couple of times. The following steps have been introduced to improve quality of theta angles:
    1. Detecting and resolving wrist gimbal lock if theta5 == 0
    2. Correcting theta3-6 to fit into (-pi, pi) range and pick closest path from previous angle to current angle
    3. Chosing solution based on sign of sin(theta5)  (sign shows in which quarter theta5 is located)

* Things to improve

1. Although speed and accuracy of arm movement was noticeably increased, it was still pretty slow. There is still a lot to improve in for theta3-6 calculation/filtering algorithm
2. IK calculation for every trajectory point takes significant time, which probably could be decreased by using nupmy as a main matrix manipulation library, rather than sympy

Below are couple of screenshots of final results:

![kuka arm final][kuka_arm0]
![kuka arm final][kuka_arm1]


