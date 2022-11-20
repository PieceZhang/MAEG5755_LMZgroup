% DH convention solver
% generate T matrix from DH table
% Written by ZHANG Yuelin, 2022.11.12

% define params and DH table
clear

%%%%% HW1 problem 3 %%%%%
disp('HW1 problem 3')
syms H
syms theta1
syms L1
syms theta2
syms d3
syms theta4
syms d5
syms theta6
syms L2
DH = [0 0 H theta1;
      pi/2 0 L1 pi/2+theta2;
      -pi/2 0 d3 pi/2;
      pi/2 0 0 -pi/2+theta4;
      0 0 d5 0;
      -pi/2 0 L1 theta6]  % DH table


% calculate
T_final = eye(4);
sizeDH = size(DH);
for row = 1:sizeDH(1)
    row
    rowDH = DH(row,:);
    alphai1 = rowDH(1);
    ai1 = rowDH(2);
    di = rowDH(3);
    thetai = rowDH(4);
    T = [cos(thetai) -sin(thetai) 0 ai1;
         sin(thetai)*cos(alphai1) cos(thetai)*cos(alphai1) -sin(alphai1) -sin(alphai1)*di;
         sin(thetai)*sin(alphai1) cos(thetai)*sin(alphai1) cos(alphai1) cos(alphai1)*di;
         0 0 0 1];
    T = simplify(T)
    T_final = T_final*T;
end
T_final = simplify(T_final)

