%% Vehicle Parameter

m = 2237; % kg
Jz = 5112; % kg/m^2
a = 1.46; % m
b = 1.55; % m

A = -6.8357;
B = 0.0325;
C = 238.9874;

load('VehicleData.mat');

%% IC
m0 = [3; 0; 0; 0; 0; 0];
Sigma0 = diag([0.1, 0.1, 0.1, 0.1, 0.1, 0.1]);
W = diag([1, 1, 0.1, 1, 1, 0.1, 1, 0.1, 0.05]);
xkk1 = m0;
Sxkk1 = Sigma0;
dT = 0.01; % sampling time

numOfMeas = size(ay_meas,1);
xHist = zeros(size(m0, 1), numOfMeas);
yHist = zeros(3, numOfMeas);
i = 1;

KF.E = [eye(6), zeros([6, 3])];

KF.F = [zeros([3, 6]), eye(3)];

dfdx = @(x, u) dfdxSymb(x, u);
% dfdu = @(x, u) dfduSymb(x, u);
dhdx = @(x, u) dhdxSymb(x, u);
% dhdu = @(x, u) dhduSymb(x, u);
fhandle = @(x, u) nonlinearVehicleModel(x, u);
hhandle = @(x, u) nonlinearVehicleModelMeasure(x,u);



while i <= numOfMeas
    currentU = [Fx_in(i); delta_in(i)];
    currentY = [vx_meas(i); ay_meas(i); yawRate_meas(i)];
    [xk1k,Sxk1k,xkk,~,~] = ...
        EKF231B(xkk1,Sxkk1,dfdx,dhdx,KF.E,KF.F,fhandle,hhandle,W,currentY,currentU);
    % [xEst, sigmaXEst] = KalmanFilter(cuttentU, currentY, xPrev, signmaXPrev, W, sys, KF);
    % Update states
    xkk1 = xk1k;
    Sxkk1 = Sxk1k;
    % Store Results
    xHist(:,i) = xkk;
    yHist(:,i) = nonlinearVehicleModelMeasure(xkk,currentU);
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
legend('result from KF', 'exp measurement')
title('Sideslip Angle')

figure(2)
plot(time, yHist(2,:));
hold on
plot(time, ay_meas);
legend('result from KF', 'exp measurement')
title('Lateral Acceleration')

figure(3)
plot(time, xHist(3,:));
hold on
plot(time, yawRate_meas);
legend('result from KF', 'exp measurement')
title('Yaw Rate')


%% Helper Function

function result = dfdxSymb(xNum, uNum)

%% Vehicle Parameter

m = 2237; % kg
Jz = 5112; % kg/m^2
a = 1.46; % m
b = 1.55; % m

A = -6.8357;
B = 0.0325;
C = 238.9874;

syms vx vy wz X Y psi_ Fx delta;
x = [
   vx
   vy
   wz
   X
   Y
   psi_
   ];

u = [
    Fx
    delta
    ];

result = jacobian(nonlinearVehicleModel(x, u),x);

result = subs(result,{vx, vy, wz, X, Y, psi_, Fx, delta},...
    {xNum(1), xNum(2), xNum(3), xNum(4), xNum(5), xNum(6), uNum(1), uNum(2)});
result = double(result);
end

function result = dhdxSymb(xNum, uNum)
m = 2237; % kg
Jz = 5112; % kg/m^2
a = 1.46; % m
b = 1.55; % m

A = -6.8357;
B = 0.0325;
C = 238.9874;

syms vx vy wz X Y psi_ Fx delta;
x = [
   vx
   vy
   wz
   X
   Y
   psi_
   ];

u = [
    Fx
    delta
    ];

result = jacobian(nonlinearVehicleModelMeasure(x, u),x);


result = subs(result,{vx, vy, wz, X, Y, psi_, Fx, delta},...
    {xNum(1), xNum(2), xNum(3), xNum(4), xNum(5), xNum(6), uNum(1), uNum(2)});
result = double(result);
end