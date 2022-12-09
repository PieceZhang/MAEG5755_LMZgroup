% DH convention solver
% generate T matrix from DH table
% Written by ZHANG Yuelin, 2022.11.12

% define params and DH table
clear

disp('CUTER')

% notation 1
% syms theta1
% syms theta2
% syms theta3
% DH = [0 0 10.18 theta1;
%       pi/2 0 0 theta2-0.1488;  % pi*8.5264/180=0.1488
%       0 19.6269 0 theta3+0.1488;
%       0 20.2 0 0]

% notation 2
% syms theta1
% syms theta2
% syms theta3
% syms l1
% syms l2
% syms l3
% syms l4
% DH = [0 0 l1 pi/2+theta1;
%       pi/2 0 0 theta2-0.1488;
%       0 sqrt(l2^2+l3^2) 0 theta3+0.1488;
%       0 l4 0 0]

% notation 3
syms theta1
syms theta2
syms theta3
syms l1
syms l2
syms l3
syms beta  % beta=0.1488
DH = [0 0 l1 pi/2+theta1;
      pi/2 0 0 theta2-beta;
      0 l2 0 theta3+beta;
      0 l3 0 0]


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

