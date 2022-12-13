syms theta1
syms theta2
syms theta3
syms l1
syms l2
syms l3
syms beta  % beta=0.1488

fk = [-sin(theta1)*(l3*cos(theta2 + theta3) + l2*cos(beta - theta2))
      cos(theta1)*(l3*cos(theta2 + theta3) + l2*cos(beta - theta2))
      l1 + l3*sin(theta2 + theta3) - l2*sin(beta - theta2)];

% simplify(diff(fk(1), theta1))

J = jacobian(fk, [theta1 theta2 theta3]);
J = simplify(J)