% assignment 1, 4.2 4.3 4.4
% Written by ZHANG Yuelin

%% 4.3

syms a0, syms a1, syms a2
syms a3, syms a4, syms a5

t = 0;
x1t1 = a0 + a1*t + a2*t^2 + a3*t^3 + a4*t^4 + a5*t^5 == -pi/6;
dx1t1 = a1 + 2*a2*t + 3*a3*t^2 + 4*a4*t^3 + 5*a5*t^4== 0;
ddx1t1 = 2*a2 + 6*a3*t + 12*a4*t^2 + 20*a5*t^3== 0;
t = 5;
x1t2 = a0 + a1*t + a2*t^2 + a3*t^3 + a4*t^4 + a5*t^5 == pi/5;
dx1t2 = a1 + 2*a2*t + 3*a3*t^2 + 4*a4*t^3 + 5*a5*t^4== 0;
ddx1t2 = 2*a2 + 6*a3*t + 12*a4*t^2 + 20*a5*t^3== 0;

[a0, a1, a2, a3, a4, a5] = solve([x1t1 dx1t1 ddx1t1 x1t2 dx1t2 ddx1t2], ...
                                            [a0, a1, a2, a3, a4, a5])

polymat3 = double([a0, a1, a2, a3, a4, a5])


%% 4.4 draw
% for 3
xt = polyval(fliplr(polymat3), 0:0.01:5);
plot(0:0.01:5, xt)
hold on

% for 2
xt = [];
for t = 0:0.01:5
    xt = [xt, blendfunc2(t)];
end
plot(0:0.01:5, xt)

%% function for 4.2
function [x] = blendfunc2(t)
    % function for 2
    tb = 1.0564;
    tf = 5;
    if t>=0 && t<tb
        x = 0.1383*t^2 -pi/6;
    elseif t>=tb && t<=(tf-tb)
        x = 0.2921*(t-1.0564) + 0.1543 - pi/6;
    elseif t>(tf-tb) && t<=tf
        x = -0.1383*(t-5)^2 + pi/5;
    end
end
