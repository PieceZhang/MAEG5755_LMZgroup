% A* algorithm
% LMZ group
clc; clear; close all;

%% location
scaler = 0.2;
loc_bigbox = [23 6 -13];loc_bigbox=formatcoor(loc_bigbox,scaler);  % (x y z R)
R_bigbox = 7*scaler;
loc_boxes = [-8 0 -24];loc_boxes=formatcoor(loc_boxes,scaler);  % (x y z R H)
R_boxes = 8*scaler;
H_boxes = 30*scaler;
loc_Penguin_right = [-6 6 -35];loc_Penguin_right=formatcoor(loc_Penguin_right,scaler);
loc_Penguin_left = [24 5 -26];loc_Penguin_left=formatcoor(loc_Penguin_left,scaler);
% loc_pigeon = [7 30 -55];loc_pigeon=formatcoor(loc_pigeon,scaler);
loc_start = [-8 30 -8];loc_start=formatcoor(loc_start,scaler);
% loc_end = [6 25 -32];loc_end=formatcoor(loc_end,scaler);

%% params
obstacleMatrix = loc_bigbox;
RobstacleMatrix = R_bigbox;
cylinderMatrix = loc_boxes(1:2);
cylinderRMatrix = R_boxes;
cylinderHMatrix = H_boxes;
start = loc_start;          % start point
goal = loc_Penguin_right;   % object
[numberOfSphere, ~] = size(obstacleMatrix);
[numberOfCylinder, ~] = size(cylinderMatrix);
Alldirec = [[1,0,0];[0,1,0];[0,0,1];[-1,0,0];[0,-1,0];[0,0,-1];...
            [1,1,0];[1,0,1];[0,1,1];[-1,-1,0];[-1,0,-1];[0,-1,-1];...
            [1,-1,0];[-1,1,0];[1,0,-1];[-1,0,1];[0,1,-1];[0,-1,1];...
            [1,1,1];[-1,-1,-1];[1,-1,-1];[-1,1,-1];[-1,-1,1];[1,1,-1];...
            [1,-1,1];[-1,1,1]];
threshold = 0.7;
stop = threshold*1.5;
g = [start, 0; goal, inf];
Path = [];
Parent = [];
Open = [start, g(findIndex(g,start),4) + getDist(start,goal)];
%% obstacle
[n,~] = size(obstacleMatrix);
for i = 1:n
    [x,y,z] = sphere();
    surfc(RobstacleMatrix(i)*x+obstacleMatrix(i,1),...
        RobstacleMatrix(i)*y+obstacleMatrix(i,2),...
        RobstacleMatrix(i)*z+obstacleMatrix(i,3));
    hold on;
end

[n,~] = size(cylinderMatrix);
for i = 1:n
    [x,y,z] = cylinder(cylinderRMatrix(i));
    z(2,:) = cylinderHMatrix(i);
    surfc(x + cylinderMatrix(i,1),y + cylinderMatrix(i,2),...
        z,'FaceColor','interp');
    hold on;
end

bar1 = scatter3(start(1),start(2),start(3),80,"cyan",'filled','o');hold on
bar2 = scatter3(goal(1),goal(2),goal(3),80,"magenta",'filled',"o");
axis equal
%% find
while ~isempty(Open)
    [xi, index] = findMin(Open);
    Open(index,:) = [];
    if getDist(xi, goal) < stop
        break;
    end
    children = getChildren(xi, Alldirec, threshold, obstacleMatrix, RobstacleMatrix,...
                           cylinderMatrix, cylinderRMatrix, cylinderHMatrix);
    scatter3(children(:,1),children(:,2),children(:,3),10,'filled','o');
    drawnow;
    [n,~] = size(children);
    for i = 1:n
        child = children(i,:);
        if findIndex(g, child) == 0 
            g = [g; child, inf];
        end
        a = g(findIndex(g, xi),4) + getDist(xi,child);
        if a < g(findIndex(g, child),4)
            g(findIndex(g, child),4) = a;
            Parent = setParent(Parent, child,xi);
            Open = setOpen(Open, child, a, goal);
        end
    end  
end
lastPoint = xi;
%% path
x = lastPoint;
Path = x;
[n,~] = size(Parent);
while any(x ~= start)
    for i = 1:n
        if Parent(i,1:3) == x
            Path = [Parent(i,4:6); Path];
            break;
        end
    end
    x = Parent(i,4:6);
end
plot3([Path(:,1);goal(1)],[Path(:,2);goal(2)],[Path(:,3);goal(3)],'LineWidth',3,'color','r');

%% convert back to simulator coor
for i=1:length(Path)
    Path(i,:) = deformatcoor(Path(i,:),scaler);
end
Path
%% func
function children = getChildren(pos, Alldirec, step,circleCenter,circleR, cylinderCenter,cylinderR, cylinderH)
allchild = [];
[n,~] = size(Alldirec);
for i = 1:n
    direc = Alldirec(i,:);
    child = pos + direc * step;
    if ~checkCol(child, circleCenter,circleR, cylinderCenter,cylinderR, cylinderH)
        continue;
    end
    allchild = [allchild; child];
end
children = allchild;
end

function flag = checkCol(pos, circleCenter,circleR, cylinderCenter,cylinderR, cylinderH)
[numberOfSphere, ~] = size(circleCenter);
[numberOfCylinder, ~] = size(cylinderCenter);
flag = true;
for i = 1:numberOfSphere
    if getDist(pos, circleCenter(i,:)) <= circleR(i)
        flag = false;
        break;
    end
end
for i = 1:numberOfCylinder
    if getDist(pos(1:2), cylinderCenter(i,:)) <= cylinderR(i) && pos(3) <= cylinderH(i)
        flag = false;
        break;
    end
end
if pos(3) <= 0 flag = false; end
end

function Par = setParent(Parent, xj, xi)
[n,~] = size(Parent);
if n == 0
    Par = [xj, xi];
else
    for i = 1:n
        if Parent(i,1:3) == xj
            Parent(i,4:6) = xi;
            Par = Parent;
            break;
        end
        if i == n
            Par = [Parent; xj, xi];
        end
    end
end
end

function Ope = setOpen(Open, child, a, goal)
[n,~] = size(Open);
if n == 0
    Ope = [child, a + getDist(child, goal)];
else
    for i = 1:n
        if Open(i,1:3) == child
            Open(i,4) = a + getDist(child, goal);
            Ope = Open;
        end
        if i == n
            Ope = [Open; child, a + getDist(child, goal)];
        end
    end
end
end

function h = heuristic(pos, goal)
h = max([abs(goal(1) - pos(1)),abs(goal(2) - pos(2)),abs(goal(3) - pos(3))]);
end

function index = findIndex(g, pos)
[n,~] = size(g);
index = 0; 
for i = 1:n
    if g(i,1:3) == pos
        index = i;
        break;
    end
end
end

function d = getDist(x,y)
d = sqrt(sum((x - y).^2));
end

function [pos, index] = findMin(Open)
[~,index] = min(Open(:,4));
pos = Open(index,1:3);
end

function [x] = formatcoor(x,scaler)
    if length(x) == 3
        x = [-x(1) -x(3) x(2)]*scaler;
    elseif length(x) > 3
        x = [-x(1)*scaler -x(3)*scaler x(2)*scaler x(4:length(x))];
    end
end

function [x] = deformatcoor(x,scaler)
    x = [-x(1) x(3) -x(2)]*(1/scaler);
end
