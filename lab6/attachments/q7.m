clc
clear
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Francesco Borrelli ME C231A 2015
% Kinematic Navigation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
N=50;
sampling=10;
%Var Defintions
z = sdpvar(2,N);

%Initial and terminal condition
z0 = [0;1];
zT = [850;1];
dzmin=-[20;2];
dzmax=[20;2];
zmin = [0; 0];
zmax = [1000; 7];

%Obstacle list
i=1;
obs{i}.center=[400;1];
obs{i}.LW=[200;2];
obs{i}.theta=0; %(in radiants)
i=i+1;
obs{i}.center=[800;5];
obs{i}.LW=[400;4];
obs{i}.theta=0; %(in radiants)

% some obtacle postprocessing 
for j=1:length(obs)
    t=obs{j}.theta;
    % generate T matrix for each obstacle
    obs{j}.T=[cos(t), -sin(t);sin(t) cos(t)]*diag(obs{j}.LW/2);
    % polyehdral representaion
    obs{j}.poly=obs{j}.T*unitbox(2)+obs{j}.center;
end

%try to remove/add this one


%Constraints
%Setup Optimization Problem
cost = 0;
%constr = [z(:,1)==z0;z(:,N)==zT;z(1,:)>=0;z(1,:)<=1000;z(2,:)>=0;z(1,:)<=7];
Q=eye(2);
constr = [z(:,1)==z0,z(:,N)==zT];
for t = 2:N
      cost=cost+(z(:,t)-z(:,t-1))'*Q*(z(:,t)-z(:,t-1));
      constr = constr +[dzmin<= z(:,t)-z(:,t-1)<=dzmax]; 
      constr = constr + [zmin<= z(:,t) <=zmax];
      for k = 0:sampling-1
          for j=1:length(obs) 
              xs=z(:,t-1)+k/sampling*(z(:,t)-z(:,t-1));
              %constr = constr + [norm(inv(obs{j}.T)*(xs - obs{j}.center),inf)>=1];
              constr = constr +[(xs-obs{j}.center)'*inv(obs{j}.T)'*inv(obs{j}.T)*(xs-obs{j}.center)>=2]; 
          end
      end
end
options = sdpsettings('solver','ipopt');
%options.ipopt=ipoptset('linear_solver','MUMPS');
solvesdp(constr,cost,options);
z_vec = double(z);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plotting Functions % to add title and labels
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
th = 0:pi/50:2*pi;
for j=1:length(obs)
    for l=1:length(th)
        z=[cos(th(l));sin(th(l))]*sqrt(2);
        y=obs{j}.T*z+obs{j}.center;
        xobs{j}(l) = y(1);
        yobs{j}(l) = y(2);
    end
end

%% plot routine
figure

plot(z_vec(1,:),z_vec(2,:),'o')
hold on
for j=1:length(obs)
plot(xobs{j}, yobs{j},'b');
plot(obs{j}.T*unitbox(2)+obs{j}.center);
end
xlim([0 1000]);
ylim([0 7]);
