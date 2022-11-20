% assignment 1, 4.1
% Written by ZHANG Yuelin


clear
tlist = [1 2 3 4 5 6];
polymat = [];

% for x
% xt = [2.3 3.156 2.829 1.771 1.444 2.3];
% dxt = [1.1309 0.3495 -0.9149 -0.9149 0.3495 1.1309];

% for y
xt = [2.6 1.7708 0.4292 0.4292 1.7708 2.6];
dxt = [0 -1.4341 -0.8863 0.8863 1.4341 0];

for i = 2:length(tlist)-1
    syms a0, syms a1, syms a2, syms a3
    syms b0, syms b1, syms b2, syms b3
    t1 = tlist(i-1);
    t2 = tlist(i);
    t3 = tlist(i+1);
    
    t = tlist(t1);
    x1t1 = a0 + a1*t + a2*t^2 + a3*t^3 == xt(t1);
    dx1t1 = a1 + 2*a2*t + 3*a3*t^2 == dxt(t1);
    t = tlist(t2);
    x1t2 = a0 + a1*t + a2*t^2 + a3*t^3 == xt(t2);
    dx1t2 = a1 + 2*a2*t + 3*a3*t^2 == dxt(t2);
    t = tlist(t2);
    x2t2 = b0 + b1*t + b2*t^2 + b3*t^3 == xt(t2);
    dx2t2 = b1 + 2*b2*t + 3*b3*t^2 == dxt(t2);
    t = tlist(t3);
    x2t3 = b0 + b1*t + b2*t^2 + b3*t^3 == xt(t3);
    dx2t3 = b1 + 2*b2*t + 3*b3*t^2 == dxt(t3);
    
    [a0, a1, a2, a3, b0, b1, b2, b3] = solve([x1t1 dx1t1 x1t2 dx1t2 x2t2 dx2t2 x2t3 dx2t3], ...
                                            [a0 a1 a2 a3 b0 b1 b2 b3]);
    
    v = vpa([a0, a1, a2, a3, b0, b1, b2, b3])

    polymat = [polymat; v(1:4)];
end 
polymat = [polymat; v(5:8)]
polymat = double(polymat);


% plot
polymat = fliplr(polymat)
t = tlist(1):0.01:tlist(end);
xtlist = zeros(size(t));

for i = 1:length(t)-1
    xtlist(i) = polyval(polymat(fix(t(i)),:), t(i));
end
xtlist(end) = polyval(polymat(end,:), t(end));

plot(t, xtlist)
