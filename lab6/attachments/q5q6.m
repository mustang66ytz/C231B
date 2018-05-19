%% Q5
% There are many limitations and draw back with this brute force rrt
% First, the q_rand is selected randomly in the safe zone, which does not
% contain any implication of the direction of the goal region.
% Second, everytime all the possible steering angles are checked to find
% out the input steering angle to the vehicle to move to the q_rand, which
% is computaionally expensive.
% Third, the resulting path can actually be optimized by triming the
% detouring, which results from the way we find q_nearest at each
% iteration, the q_nearest is only the closest to the q_rand, but a better 
% path may not pass the node necessarily.

% There are some modifications that could be make to improve the
% perforamance of the naive rrt algorithm: 
% First, modify the q_rand generating algorithm by including heuristics
% concerns. So the generated random points are more efficient to connect
% the nodes into the goal area.
% Second, analyzing the kinematics model of the vehicle to generate the
% optimal steering angle at each iteration, rather than try out all the
% possible steering angles as input.
% Third, implementing RRT* algorithm rather than RRT to remove the sub-opt
% imal intermiate points.