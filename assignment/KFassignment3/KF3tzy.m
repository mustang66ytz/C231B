%% Q2
%% (a)
% write out the state space model of the system, W is the variance of w(k)
deltaT = 0.1;
A = [1 deltaT 0; 0 1 deltaT; 0 0 1];
E = [0 0 0; 0 0 0; 1 0 0];
C = [1 0 0; 0 0 1];
F = [0 1 0; 0 0 1];
W = [0.3 0 0; 0 1 0; 0 0 1];
TS = deltaT;
% Find the kalman filter gain L:
[K,L,H,G,nIter,estErrVar] = formSteadyStateKFtoFinish(A,E,C,F,W,TS);
%result = ['The result is: [' num2str(L(:).') ']'];
disp("The converged kalman filter gain is: ")
disp(L)
disp("The steady state error variance is: ")
disp(estErrVar)

%% (b) 
% to study the variance of velocity, we need to extract the entry
% describing the variance of steady-state error variance calculated above
errVelEst = estErrVar(2,2);
sigmaErrVelEst = sqrt(errVelEst);
beta = 0.73/sigmaErrVelEst;
prob = 1/(beta*beta);
disp("The probability for the error of x2 greater than 0.73 is:")
disp(prob)
disp("According to the Chebychev inequality, the probability id smaller than 0.253");

%% (c)
N = 200;
repeat = 600;
% Create System
% The Kalman Filter code below assumes that the state-space arrays (A, E,
% C, F) represent time-varying dynamics, and hence should be defined as 3-d
% arrays.
nX = size(A, 1);
nW = size(W, 1);
nY = size(C, 1);
arrayA = repmat(A, [1,1,N]);  % nX-by-nX-by-N
arrayE = repmat(E, [1,1,N]);  % nX-by-nW-by-N
arrayC = repmat(C, [1,1,N]);   % nY-by-nX-by-N
arrayF = repmat(F, [1,1,N]);  % nY-by-nW-by-N

arrayXseq = zeros(nX,N+1,repeat); %nX-by-N-by-repeat
estErrX2 = zeros(repeat,1);
% Create/Declare variance of disturbance and initial condition (mean too)
% The variance of the "w" sequence is allowed to be time-varying, and
% should be defined as a 3-d array.
arraySW = repmat([0.3 0 0; 0 1 0; 0 0 1], [1,1,N]); % nW-by-nW-by-N
Sx0 = eye(nX, nX);      % nX-by-nX
m0 = zeros(nX,1);     % nX-by-1

% Initialize KF states with appropriate values
Sxii1 = Sx0;
xii1 = m0;

% Create a specific initial condition and noise sequence
% Under ideal circumstances, this should be consistent with the statistical
% assumptions made in the previous code cell.  When studying robustness,
% namely how the filter performance degrades as assumptions are not met, it
% may be useful to create an initial condition and noise sequence which is
% not consistent with the assumptions

% Create a specific reference input
%arrayR =   % nR-by-1-by-N
count = 0;
totalwSeq = zeros(nW,N,repeat);
%% (d) Simulate the system/KF one step at a time
for j=1:repeat
    %wSeq = repmat([0.3 0 0; 0 1 0; 0 0 1]*randn(nW,1), [1,1,N]);    % nW-by-1-by-N
    wSeq = [0.3 0 0; 0 1 0; 0 0 1]*randn(nW,N); %nW-by-N
    totalwSeq(:,:,j) = wSeq;
    x0 = Sx0*randn(nX,1) + m0;     % nX-by-1
    emptyB = [];  % this template file is for the case of no control signal
    emptyu = [];  % this template file is for the case of no control signal
    y = zeros(nY,N);
    xSeq = zeros(nX,N);
    xSeq(:,1) = x0;
    xiiSeq = zeros(nX,N);
    for i=0:N-1
       iMatlab = i+1;
       Ai = arrayA(:,:,iMatlab);
       Ei = arrayE(:,:,iMatlab);
       Ci = arrayC(:,:,iMatlab);
       Fi = arrayF(:,:,iMatlab);
       Swi = arraySW(:,:,iMatlab);
       wi = wSeq(:,iMatlab);
       xi = xSeq(:,iMatlab);
       % Get y(i) from system model, using x(i) and w(i)
       y(:,iMatlab) = Ci*xi + Fi*wi;
       % Get estimates of x(i+1|i) using u(i) and y(i) from KF
       [xi1i,Sxi1i,xii,Sxii,Syii1,Lk,Hk,Gk,wkk] = ...
          KF231B(xii1,Sxii1,Ai,emptyB,Ci,Ei,Fi,Swi,emptyu,y(:,iMatlab));
       % Save the state estimate xHat_{i|i}
       xiiSeq(:,iMatlab) = xii;
       % Other estimates can be saved too.  KF231B returns additional estimates
       % (for example, wkk), and saving those sequences could be useful too,
       % depending on the computational exercise.
       %
       % Get x(i+1) from system model, using x(i), u(i) and w(i)
       xSeq(:,iMatlab+1) = Ai*xi + Ei*wi;
       % Shift the error-variance estimate so that when loop-index i advances,
       % the initial condition for this variance is correct.
       Sxii1 = Sxi1i;
       xii1 = xi1i;
    end
    % find the estimation error for each test among 600 repetitions
    estErrX2(j) = abs(xiiSeq(2,N) - xSeq(2,N));
    if estErrX2(j)>=0.73
        count = count + 1;
    end
    arrayXseq(:,:,j) = xSeq;
end
%% (e)
figure;
histogram(estErrX2);
disp("the number of estimation error greater than 0.73 are, which is smaller than 25% of all 600 simulations")
disp(count/repeat);

%% Q3
%% (a)
nX = 4;
nY = 3;
nU = 2;
H = drss(nX,nY,nU);
if max(abs(eig(H.A)))>0.999
    H.A = 0.99*H.A; 
end
isstable(H)

%% (b)
f = @(Omega) freqresp(H,Omega)*freqresp(H,Omega)';
Int = 1/(2*pi) * integral(f,0,2*pi,'ArrayValued',true);
Int = real(Int);

%% (c)
sqrt(trace(Int));
norm(H,2)

%% (d)
% Decide on the duration (in time) of each input/output sequences.
% The formula has k -> INF, but we can only simulate for a finite
% duration. In general, choose a duration "long" compared to the
% time-constants of the system. Here, we simply pick a duration
% of 200 steps.
Nsteps = 200;

% So each instance (associated with a point in the sample space)
% of the u-sequence is an nU-by-200 array. The convention in
% Matlab?s simulation codes (lsim, Simulink, etc) is to
% represent this as a 200-by-nU array.

% Decide on the dimension of the sample-space. Take 2000 outcomes,
NsampleSpace = 2000;

% Hence, we need 2000, 200-by-nU arrays, as generated by RANDN.
% Create a 200-by-nU-by-2000 array with randn.
U = randn(Nsteps,nU,NsampleSpace);

% Initialize a nY-by-nY matrix to hold the expectation
% of the final value (of each simulation) of y*y^T.
E_y_yT = zeros(nY,nY);

% Our construction (from RANDN) had all outcomes equally likely,
% with probability equal to 1/NsampleSpace. So, create 
% constant value for each individual outcome?s probability
p = 1/NsampleSpace;

% Use a FOR loop to compute the response for all the
% sample-space outcomes
for i=1:NsampleSpace
% Use LSIM to compute each response.
    Y = lsim(H,U(:,:,i));
    y_end  = Y(end,:)';
    E_y_yT = E_y_yT + p*(y_end*y_end');
end
% Compare to limiting expected value to the integral
E_y_yT
Int
%% (e)
covar(H, eye(nU))

%% Q4
% (a)

% (b)
% Please refer to the hand-written derivations

% (c)
% Refer to the ekfIterationScript.m

% (d)
% By running the ekfIterationScript.m many times, it is found that the
% estimation errors for x1 and x2 are decreasring overtime, and almost
% always bounded by the upper bound and lower bound. However, the
% estimation error for x3 increases as the time increasing, and exceeds the
% boundaries. 
    
% (e)
% The initial condition is set using the variance and the mean of the init-
% ial conditions, where variance multiplies a nX-by-1 matrix with
% uniformally distributed elelments within -0.5 and 0.5 plus the mean of.

% By examine the noise sequence w, we see that w1 is much smaller than w2,
% which results from the relative big magnitude difference between SW(1,1) 
% and SW(2,2). Intuitively, it implies that the state noise is much smaller
% than the measurement noise.

% (i)
% The initial condition is set randomly, and distributed uniformly.

% (ii)
% The noise sequence is set randomly based on its standard varianca, and is
% distributed normally as a white noise, which is consisitent with the
% probabilistic model of noise.

% (iii)
% By examining the initialization parameters such as altitude, velocity,
% and scale factor, it is found that the variances of x3 is huge compared
% to that of x1 or x2, thus the EKF is designed to trust x3 lessly. and
% thus, the EKF does not learn enough about x3 at the beginning stage of
% the simulation.

%% Q5

% (a)
% General intuition in assisting the verification of the realization
% formulation: 
% State evolution:
%   xG(k+1) = AG*xG(k) + BG*q;
%   xP(k+1) = AP*xP(k) + BP*u; u = CG*xG(k) + DG*q;
%   xP(k+1) = AP*xP(k) + BP*CG*xG(k) + BP*DG*q;
%   xK(k+1) = AK*xK(k) + BK*u; u = CP*xP(k) + DP*(CG*xG(k) + DG*q);
%   xK(k+1) = AK*xK(k+1) + BK*CP*xP(k) + BK*DP*CG*xG(k) + BK*DP*DG*q;
% Output:
%   e(k) = xP(k) - xhatP(k); xhatP(k) = M*CK*xK(k) + M*DK*u; u = CP*xP(k) + DP*CG*xG(k)
%   e(k) = xP(k) - M*CK*xK(k) - M*DK*CP*xP(k) - M*DK*DP*CG*xG(k);
%   e(k) = (I - M*DK*CP)*xP(k) -M*DK*DP*CG*xG(k) -M*CK*xK(k);






