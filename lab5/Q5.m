%% Extended Kalman Filter

%m = 2237; % kg
%Jz = 5112; % kg/m^2
%a = 1.46; % m
%b = 1.55; % m

%A = -6.8357;
%B = 0.0325;
%C = 238.9874;
load('VehicleData.mat');

m0 = [3; 0; 0; 0; 0; 0];
Sx0 = diag([0.1, 0.1, 0.1, 0.1, 0.1, 0.1]);
W = diag([1, 1, 0.1, 1, 1, 0.1, 1, 0.1, 0.05]);
xkk1 = m0;
Sxkk1 = Sx0;
dt = 0.01; % sampling time

N = size(ay_meas,1);
xHist = zeros(size(m0, 1), N);
yHist = zeros(3, N);
i = 1;

E = [eye(6), zeros([6, 3])];
F = [zeros([3, 6]), eye(3)];

dfdx = @(x, u) dfdxsys(x, u);
dhdx = @(x, u) dhdxsys(x, u);
fhandle = @(x, u) discreteModel(x, u);
hhandle = @(x, u) nonlinearVehicleModelMeasure(x,u);

%%
while i <= N
    uk = [Fx_in(i); delta_in(i)];
    yk = [vx_meas(i); ay_meas(i); yawRate_meas(i)];
    [xk1k,Sxk1k,xkk,~,~] = ...
        EKF231B(xkk1,Sxkk1,dfdx,dhdx,E,F,fhandle,hhandle,W,yk,uk);
    
    xkk1 = xk1k;
    Sxkk1 = Sxk1k;
   
    xHist(:,i) = xkk;
    yHist(:,i) = nonlinearVehicleModelMeasure(xkk,uk);
    i = i + 1;
end


%% Plot

vYSeq = xHist(2,:);
vXSeq = xHist(1,:);

betaSeq = atan(vYSeq./vXSeq);

figure(1)
plot(time, betaSeq);
hold on
plot(time, beta_meas);
legend('sideslip estimation', 'sideslip measurement')
title('Sideslip Angle')

figure(2)
plot(time, yHist(2,:));
hold on
plot(time, ay_meas);
legend('lateral acceleration estimation', 'lateral acceleration measurement')
title('Lateral Acceleration')

figure(3)
plot(time, xHist(3,:));
hold on
plot(time, yawRate_meas);
legend('yaw rate estimation', 'yaw rate measurement')
title('Yaw Rate')


%% Helper Function
function result = dfdxsys(zVal, uVal)
syms vx vy wz X Y psi_ Fx delta;
z = [vx; vy; wz; X; Y; psi_];
u = [Fx; delta];

result = jacobian(discreteModel(z, u),z);
result = subs(result,{vx, vy, wz, X, Y, psi_, Fx, delta},...
    {zVal(1), zVal(2), zVal(3), zVal(4), zVal(5), zVal(6), uVal(1), uVal(2)});
result = double(result);
end

function result = dhdxsys(zVal, uVal)
syms vx vy wz X Y psi_ Fx delta;
x = [vx vy wz X Y psi_];
u = [Fx delta];

result = jacobian(nonlinearVehicleModelMeasure(x, u),x);
result = subs(result,{vx, vy, wz, X, Y, psi_, Fx, delta},...
    {zVal(1), zVal(2), zVal(3), zVal(4), zVal(5), zVal(6), uVal(1), uVal(2)});
result = double(result);
end