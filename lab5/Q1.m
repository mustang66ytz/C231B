%% Q1
%% constant velocity skidpad
u1 = [2500; 10*pi/180];
dzdt = @(t, z) nonlinearVehicleModel(z, u1);
[tsol1, zsol1] = ode45(dzdt,[0 20],[20; 0; 0; 0; 0; 0]);
% plot X Y over time
figure;
plot(zsol1(:,4), zsol1(:,5));
xlabel('X')
ylabel('Y')
% plot the lateral and longitudinal velocity over time
figure;
plot(tsol1,zsol1(:,1),'-o',tsol1,zsol1(:,2),'-o');
xlabel('time')
ylabel('velocity')
legend('longitudinal velocity', 'lateral velocity')
%% constant acceleration skidpad
xinit2 = [5;0;0;0;0;0];
u2 = [8000; 10 / 180 * pi];
dzdt = @(t, z) nonlinearVehicleModel(z, u2);
[tsol2, zsol2] = ode45(dzdt, [0 20], xinit2);
% plot X and Y over time
figure
plot(zsol2(:,4), zsol2(:,5));
xlabel('X')
ylabel('Y')
% plot the lateral and longitudinal velocity over time
figure;
plot(tsol2, zsol2(:,1), '-o', tsol2, zsol2(:,2), '-o');
xlabel('time')
ylabel('velocity')
legend('longitudinal velocity', 'lateral velocity')
%% lane change
xinit3 = [20;0;0;0;0;0];
u3 = @(t) (0 * (t > 0 & t <= 5) + 5 * (t > 5 & t <= 7) + 0 * (t > 7 & t <= 12) - 5 * (t > 12 & t <= 14) + 0 * (t > 14 & t <= 20)) / 180 * pi;
dzdt = @(t,z) nonlinearVehicleModel(z, [0; u3(t)]);
[tsol3, zsol3] = ode45(dzdt, [0 20], xinit3);
% plot X and Y over time
figure
plot(zsol3(:,4), zsol3(:,5));
xlabel('X')
ylabel('Y')
% plot the lateral and longitudinal velocity over time
figure;
plot(tsol3, zsol3(:,1), '-o', tsol3, zsol3(:,2), '-o');
xlabel('time')
ylabel('velocity')
legend('longitudinal velocity', 'lateral velocity')

%% Q2
% forward euler approximation
zsold = zeros(6, 2001);
zsold(:,1) = [20;0;0;0;0;0];
Ts = 0.01;
tspand = 0:0.01:20;
ud =[repmat([0; 0], 1, 500),repmat([0; 5 / 180 * pi], 1, 200),...
    repmat([0; 0], 1, 500), repmat([0; -5 / 180 * pi], 1, 200),... 
    repmat([0; 0], 1, 601)];
for k = 1: 1: 2000
    zsold(:,k+1) = zsold(:,k) + Ts * nonlinearVehicleModel(zsold(:,k), ud(:,k)); 
end
% plot X and Y over time
figure
plot(zsold(4,:), zsold(5,:));
xlabel('X')
ylabel('Y')
% plot the lateral and longitudinal velocity over time
figure;
plot(tspand, zsold(1,:), '-o', tspand, zsold(2,:), '-o');
xlabel('time')
ylabel('velocity')
legend('longitudinal velocity', 'lateral velocity')

%% Q3
% symbolic linearization
% define the equilibrium point first
equil = [10;0;0;0;0;0;0;0];
% define all the symbolic variables
syms vx vy wz X Y phi Fx delta;
zsyms = [vx; vy; wz; X; Y; phi];
usyms = [Fx; delta];

dzdt = nonlinearVehicleModel(zsyms, usyms);
y = nonlinearVehicleModelMeasure(zsyms, usyms);
% find the jacobians at the linearization point (equil)
dzdtJac = jacobian(dzdt, [zsyms; usyms]);
dzdtLinear = subs(dzdtJac, [zsyms; usyms], equil);
yJac = jacobian(y, [zsyms; usyms]);
yLinear = subs(yJac, [zsyms; usyms], equil);
% convert the symbolic expression back to the numeric answer
ABmatrices = double(dzdtLinear);
A = ABmatrices(:, 1:6);
B = ABmatrices(:, 7:8);
CDmatrices = double(yLinear);
C = CDmatrices(:, 1:6);
D = CDmatrices(:, 7:8);

sysLinearDiscrete = ss(A, B, C, D);
save('sysLinearDiscrete')
sysLinearDiscrete

%% Q4
% load the data
load('sysLinearDiscrete.mat')
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
arrayA = repmat(A, [1,1,N]);  % nX-by-nX-by-N
arrayB = repmat(B, [1,1,N]);  % nX-by-nU-by-N
arrayE = repmat([eye(nX) zeros(nX, nW-nX)], [1,1,N]);  % nX-by-nW-by-N
arrayC = repmat(C, [1,1,N]);  % nY-by-nX-by-N
arrayD = repmat(D, [1,1,N]);  % nU-by-nY-by-N
arrayF = repmat([zeros(nY, nX) eye(nW-nX)], [1,1,N]);   % nY-by-nW-by-N
% initialize the kalman filter
Sxii1 = Sx0;
xii1 = m0;
xiiSeq = zeros(nX,N);
yiiSeq = zeros(nY, N);
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
      KF231B(xii1,Sxii1,Ai*dt+eye(nX),Bi*dt,Ci,Di,Ei,Fi,Swi,uInput(:,iMatlab),yMeasure(:,iMatlab));
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
   yiiSeq(:,iMatlab) = Ci*xii + Di*uInput(:, iMatlab) + Fi*wkk;
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
figure
plot(0:N-1, yiiSeq(3,:), 'r', 0:N-1, yawRate_meas, 'b');
legend("estimated yaw rate", "measured yaw rate");

