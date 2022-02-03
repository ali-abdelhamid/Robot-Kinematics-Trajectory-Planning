# Robot-Kinematics-&-Trajectory-Planning
Custom MATLAB Library for Kinematic and Dynamic Modeling of n-axis robotic manipulators

Script Breakdown:
1. Home Position Homogeneous Transformation Matrix
2. S-List (screw axis) based on Denavit Hartenberg Convention
3. Product of Exponentials for Forward Kinematics using skew symmetric representation of S-list
4. Compute the Jacobian (matrix relating the joint velocities to the end effector velocities)
5. Inverse Kinematics calculation for initial and final joint values
6. Plotting Task Space Trajectory Equations using the polynomial method where the order of the polynomial is defined by the total number of constraints in terms of position and speed
7. Iterative, non-linear Inverse Kinematics solver using the inverse of the jacobian for an approximate derivative (speed) of the joint.
8. This iterator is nested in a for loop that runs (in order) over all task space positions, passing each incremental step on to the Inverse Kinematics solver until a full list of joint values is derived
9. This list of joint values can be used for joint-space trajectory mapping and for dynamics based modeling for specific joint torques
