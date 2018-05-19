clc
clear
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Francesco Borrelli ME C231A 2015
% Kinematic Navigation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
N=20;
sampling=10;
%Var Defintions
z = sdpvar(2,N);
%Initial and terminal condition
z0 = [0;0];
zT = [10;10];
dzmin=-[2;2];
dzmax=[2;2];
zmin=-[20;20];
zmax=[20;20];
%Obstacle list
i=1;
obs{i}.center=[4;4];
obs{i}.LW=[2;2];
obs{i}.theta=30*pi/180; %(in radiants)
i=i+1;
obs{i}.center=[5.5;7];
obs{i}.LW=[1;1];
obs{i}.theta=0*pi/180; %(in radiants)
% i=i+1;
% obs{i}.center=[9;4];
% obs{i}.LW=[2;2.5];
% obs{i}.theta=-30*pi/180; %(in radiants)
% some obstacle postprocessing
for j=1:length(obs)
t=obs{j}.theta;
% generate T matrix for each obstacle
obs{j}.T=[cos(t), -sin(t);sin(t) cos(t)]*diag(obs{j}.LW/2);
% polyehdral representaion
obs{j}.poly=obs{j}.T*unitbox(2)+obs{j}.center;
end
%try to remove/add this one
%z_obs{4}=[3;7];
%d_obs{4}=8;
%Qobs{4}=diag([1,10]);
%Constraints
%Setup Optimization Problem
cost = 0;
constr = [z(:,1)==z0;z(:,N)==zT];
Q=eye(2);
%constr = [zmin<=z(:,N)<= zmax, z(:,1)==z0,z(:,N)==zT];
for t = 2:N
cost=cost+(z(:,t)-z(:,t-1))'*Q*(z(:,t)-z(:,t-1));
constr = constr +[dzmin<= z(:,t)-z(:,t-1)<=dzmax];
constr = constr +[zmin<= z(:,t)<=zmax];
for k = 0:sampling-1
for j=1:length(obs)
zs=z(:,t-1)+k/sampling*(z(:,t)-z(:,t-1));
[H,K]=double(obs{j}.poly);
constr = constr +[ H(1,:)*(zs)>=K(1) ...
| H(2,:)*(zs)>=K(2) ...
| H(3,:)*(zs)>=K(3) ...
| H(4,:)*(zs)>=K(4) ];
end
end
end
options = sdpsettings('solver','gurobi');
%options.ipopt=ipoptset('linear_solver','MUMPS');
solvesdp(constr,cost,options);
z_vec = double(z);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plotting Functions % to add title and labels
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% plot routine
figure
plot(z_vec(1,:),z_vec(2,:),'o')
hold on
for j=1:length(obs)
plot(obs{j}.T*unitbox(2)+obs{j}.center);
end