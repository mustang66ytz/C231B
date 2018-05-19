%% Load data
load('sysLinearDiscrete.mat')
load('VehicleData.mat')
%% Parameters 
m = 2237;
Jz = 5112;
a = 1.46;
b = 1.55;
A = -6.8357;
B = 0.0325;
C = 238.9874;
NIter = size(time,2);
%% Get linearized dynamics model
load('SYS.mat')
Alin = SYS.A;
Blin = SYS.B;
Clin = SYS.C;
Dlin = SYS.D;
%% Simulate the system/KF one step at a time

load('VehicleData.mat')
% define all the dimensions
N = 9120;
nX = 6;
nU = 2;
nW = 9;
nY = 3;
W = eye(9);
W(3,3) = 0.1;
W(6,6) = 0.1;
W(8,8) = 0.1;
W(9,9) = 0.05;
arraySW = repmat(W, [1,1,N]);
Sx0 = 0.1*eye(nX, nX);
m0 = zeros(nX,1);
m0(1) = 3;
% construct the measurement data
yMeasure = zeros(nY, N);
yMeasure(1, :) = vx_meas';
yMeasure(2, :) = ay_meas' + (vx_meas.*yawRate_meas)';
yMeasure(3, :) = yawRate_meas';
uInput = [Fx_in'; delta_in'];
% define the system
arrayA = repmat(Alin, [1,1,N]);  % nX-by-nX-by-N
arrayB = repmat(Blin, [1,1,N]);  % nX-by-nU-by-N
arrayE = repmat([eye(nX) zeros(nX, nW-nX)], [1,1,N]);  % nX-by-nW-by-N
arrayC = repmat(Clin, [1,1,N]);  % nY-by-nX-by-N
arrayD = repmat(Dlin, [1,1,N]);  % nU-by-nY-by-N
arrayF = repmat([zeros(nY, nX) eye(nW-nX)], [1,1,N]);   % nY-by-nW-by-N
% initialize the kalman filter
Sxii1 = Sx0;
xii1 = m0;
xiiSeq = zeros(nX,N);
yiiSeq = zeros(nY,N);
% run the kalman filter
for i=0:N-1
   iMatlab = i+1;
   Ai = arrayA(:,:,iMatlab);
   Bi = arrayB(:,:,iMatlab);
   Ei = arrayE(:,:,iMatlab);
   Ci = arrayC(:,:,iMatlab);
   Di = arrayD(:,:,iMatlab);
   Fi = arrayF(:,:,iMatlab);
   Swi = arraySW(:,:,iMatlab);
   %wi = wSeq(:,iMatlab);
   %xi = xSeq(:,iMatlab);
   % Get r(i), value of reference input
   %ri = arrayR(:,iMatlab);
   % Get u(i) from controller output map
   %ui = outCHan(ni,ri);
   % Get y(i) from system model, using x(i) and w(i)
   %y(:,iMatlab) = Ci*xi + Fi*wi;
   % Get estimates of x(i+1|i) using u(i) and y(i) from KF
   [xi1i,Sxi1i,xii,Sxii,Syii1,Lk,Hk,Gk,wkk] = ...
      KF231B(xii1,Sxii1,Ai,Bi,Ci,Ei,Fi,Swi,uInput(:,iMatlab),yMeasure(:,iMatlab));
   % Save the state estimate xHat_{i|i}
   xiiSeq(:,iMatlab) = xii;
   % Get x(i+1) from system model, using x(i), u(i) and w(i)
   %xSeq(:,iMatlab+1) = Ai*xi + Bi*uInput(:,iMatlab) + Ei*wi;
   % Save any signals or variances that will be used in later calculations.
   %
   % Shift the error-variance estimate so that when loop-index i advances,
   % the initial condition for this variance is correct.
   Sxii1 = Sxi1i;
   xii1 = xi1i;
   yiiSeq(:,iMatlab)= Ci*xii + Di*uInput(:,iMatlab) + Fi * wkk;
end

% plot the result
% sideslip angle
beta_hat = atan(xiiSeq(2,:)./xiiSeq(1,:));
figure;
plot(0:N-1, beta_hat);
hold on;
plot(0:N-1, beta_meas');
legend("estimated sideslip angle","measured sideslip angle");

% lateral acceleration
figure;
plot(0:N-1, yiiSeq(2,:)-yiiSeq(1,:).*yiiSeq(3,:), 'r', 0:N-1, ay_meas', 'b');
legend("estimated lateral acceleration", "measured lateral acceleration");

% yaw rate
figure;
plot(0:N-1, yiiSeq(3,:)', 'r', 0:N-1, yawRate_meas, 'b');
legend("estimated yaw rate", "measured yaw rate");
