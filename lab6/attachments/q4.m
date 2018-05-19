%
clearvars
close all

% x and y axises limit from 0 to x_max and 0 to y_max respectively.
x_max = 1000;
y_max = 7;
phi_max = pi;
phi_min = -pi;

% time of moving at every step
EPST = 0.1;

% maximum iterations
numNodes = 6000;

% attributions of starting point
q_start.coord = [0 2 0]';
q_start.cost = 0;
q_start.parent = 0; %.parent means the index of parent node

% initialize the tree
nodes(1) = q_start;

% plot the safe area
figure(1)
x=[0 300 300 500 500 1000 1000 600 600 200 200 0]; %x coordinates of all the vertices
y=[0 0 4 4 0 0 3 3 7 7 3 3];  %y coordinates of all the vertices
X=[x,x(1)];   %????????????????????
Y=[y,y(1)];   %??
plot(X,Y,'k')  %?????
fill(x,y,'r')  % fill the safe zone with color
hold on

% plot the goal area
figure(1)
axis([0 x_max 0 y_max])
goal_area = rectangle('Position',[900,1,50,0.5],'FaceColor',[0 .5 .5]);
xlabel('x')
ylabel('y')
hold on

%% define the vehicle parameters
vx = 30;
L = 3;

%% grow the tree

for i = 1:1:numNodes
    
    %q_rand = [floor(rand(1)*x_max) floor(rand(1)*y_max*2)-y_max floor(rand(1)*2*phi_max)-phi_max]';
    %plot(q_rand(1), q_rand(2), 'x', 'Color',  [0 0.4470 0.7410])
    
    pan = 0;
    % generate the random points in the given safe area and plot the points
    while ~pan
        q_rand = [rand*x_max rand*y_max floor(rand(1)*2*phi_max)-phi_max];
        pan = inpolygon(q_rand(1),q_rand(2),X,Y);
    end
    plot(q_rand(1), q_rand(2), 'x', 'Color',  [0 0.4470 0.7410])
    
    % Find the nearest point existing on the tree to the random point
    ndist = [];
    for j = 1:1:length(nodes)
        n = nodes(j);
        tmp = dist(n.coord(1:2), q_rand);
        ndist = [ndist tmp];
    end
    [mini_distance, idx] = min(ndist);
    q_nearest = nodes(idx);
        
    %brute force to check all the possible steering angles, and assign the
    %closet to q_new
   
    k = 1;
    tempdist = [];
    q_newPossible = [];
    for delta = -20:2:20
        deltaRad = delta*pi/180;
        dxdt = @(t,x) kinematicsModel(x, deltaRad, vx, L);
        [tsol, xsol] = ode45(dxdt,[0,EPST],q_nearest.coord);
        
        InorOn = all(inpolygon(xsol(:,1),xsol(:,2),X,Y)) && all(xsol(:,3)<=pi) && all(xsol(:,3)>=-pi);
        
        k=k+1;
        if InorOn == 1
            q_newPossible = [q_newPossible; xsol(end,:)];
            tempdist = [tempdist dist(q_newPossible(end,1:2), q_rand)];
        end
        
    end
    [mini_distance2, idx2] = min(tempdist);
    
    if isempty(q_newPossible(idx2,:))
        continue;
    end
    
    q_new.coord = q_newPossible(idx2,:);
    line([q_nearest.coord(1), q_new.coord(1)], [q_nearest.coord(2), q_new.coord(2)],...
        'Color', 'k', 'LineWidth', 2);
    drawnow
    hold on
    q_new.cost = dist(q_new.coord, q_nearest.coord) + q_nearest.cost;
    q_new.parent = idx;
 
    InorOn = inpolygon(q_new.coord(1),q_new.coord(2),X,Y);
    % Append to nodes
    if InorOn == 1
        nodes = [nodes q_new];
    end
    % Break if the link from second to last node to last node intersects any of
    % the four edges of the goal area
    if ~noCollision(q_nearest.coord, q_new.coord(1:2), [900,1,50,0.5]) && q_new.coord(3)>=-pi/6 && q_new.coord(3)<=pi/6 
        break
    end
end

q_end = q_new;
num_node_path = 1;

while q_end.parent ~= 0
    start = q_end.parent;
    line([q_end.coord(1), nodes(start).coord(1)], [q_end.coord(2), nodes(start).coord(2)],...
        'Color', 'r', 'LineWidth', 2);
    hold on
    q_end = nodes(start);
    num_node_path = num_node_path+1;
end

%% total number of node in the tree
num_node_tree = length(nodes)

%% number of nodes in the sequence that reaches goal area
num_node_path

