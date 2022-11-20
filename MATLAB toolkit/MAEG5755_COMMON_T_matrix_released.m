% Transform matrix
% written by ZHANG Yuelin

% degree
angle1 = 90;
angle2 = 180;
angle3 = -90;

% radius 
theta1 = (angle1) * pi / 180;  % Z rotation
theta2 = (angle2) * pi / 180;  % X rotation
theta3 = (angle3) * pi / 180;  % X rotation

R1 = [cos(theta1) -sin(theta1) 0;
      sin(theta1) cos(theta1) 0;
      0 0 1]  % Z rotation
R2 = [1 0 0;
      0 cos(theta2) -sin(theta2);
      0 sin(theta2) cos(theta2)]  % X rotation
R3 = [1 0 0;
      0 cos(theta3) -sin(theta3);
      0 sin(theta3) cos(theta3)]  % X rotation

t1 = [0;0;0];  % Translation 1
t2 = [0;0;0];  % Translation 2
t3 = [0;0;0];  % Translation 3

T01 = [R1,t1;0 0 0 1]  % Transform 01
T02 = [R1*R2,t2;0 0 0 1]  % Transform 02
T03 = [R1*R2*R3,t3;0 0 0 1]  % Transform 03
