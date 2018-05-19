%% 1.c(ii)

Amat = repmat(1,[1 1 20]);
nX = size(Amat,1);
Emat = repmat(0,[1 1 20]);
Cmat = repmat(1,[1 1 20]);
Fmat = repmat(1,[1 1 20]);
sX = 1000; sW = 1; m0=2;
for k=1:6
   [LkBatch,VkBatch,eVar] = ...
           batchKF(Amat,Emat,Cmat,Fmat,sX,sW,k);
   LkBatch(end-2*nX+1:end-nX,:);
   VkBatch(end-2*nX+1:end-nX);
end
%(xhat_(k-1|k-1)) = L(x[k]|y[k-1]) = LkBatch*Y[k-1] + VkBatch*m0

%% 1.c(iv)Repeat for the following design
Amat = repmat(1,[1 1 20]);
nX = size(Amat,1);
Emat = repmat(0,[1 1 20]);
Cmat = repmat(1,[1 1 20]);
Fmat = repmat(1,[1 1 20]);
sX = 1000; sW = 5; m0=4;
for k=1:6
   [LkBatch,VkBatch,eVar] = ...
           batchKF(Amat,Emat,Cmat,Fmat,sX,sW,k);
   LkBatch(end-2*nX+1:end-nX,:);
   VkBatch(end-2*nX+1:end-nX);
end

% From the results obtained by changing sX and sW, and therefore, their
% relative ratio, we changed the property of the system, where a bigger sX
% represents a bigger variance of the initial condition, and a  bigger sW
% represents a bigger variance of the noises. The optimal linear estimator
% is designed to base purely on the measurement and the noise. Therefore,
% the estimator should put more weight on the more reliable portion of the
% data, thus, the term with less variance should result in a bigger
% coefficient in the estimator design. From the result, indeed, we have
% that LkBatch in (ii) is bigger than that of (iv) due to the more reliable
% measurement data, and the VkBatch in (ii) is less than that of (iv) due
% to the less reliable initial condition.

%% 3
clc;
clear;

A = [0.6 -0.8; 0.7 0.6];
B = [];
C = [1 1];
E = [0.5; 0.5];
F = 1;

xkk1 = [10; 4];
uk = [];
yk = 15;

Sxkk1 = 200;
Swk = 2;

% save the output after the iterations
saveX = zeros(2, 100);
saveLk = zeros(2, 100);

for i = 1:100
    [xk1k, Sxk1k, xkk, Sxkk, Sykk1, Lk, Hk, Gk, wkk] = ...
        KF231B(xkk1, Sxkk1, A, B, C, E, F, Swk, uk, yk);
    saveLk(:,i) = Lk; 
    saveX(:,i) = xk1k;
    xkk1 = xk1k;
    Sxkk1 = Sxk1k;
    
end

figure;
plot(linspace(1,100), saveX(1,:))
hold on
plot(linspace(1,100), saveX(2,:))
legend("x1", "x2")

figure;
plot(linspace(1,100), saveLk(1,:))
hold on
plot(linspace(1,100), saveLk(2,:))
legend("Lk1", "Lk2");
% Indeed, the estimated values of state converges


%% 4
%% Separating Signal from Noise
% ME C231B, UC Berkeley, Spring 2018

%% Sample Time
% For this discrete-time example, set TS=-1. In Matlab, this
% just means an unspecified sampling time, totally within the
% context of pure discrete-time systems.
TS = -1;

% Create high-pass filter
P1 = 0.4*tf([.5 -.5],[1 0],TS);
[A1,E1,C1,F1] = ssdata(P1);
nX1 = size(A1,1); % will be 1
nW1 = size(E1,2); % will be 1

% Create low-pass filter
P2 = tf(.04,[1 -.96],TS);
[A2,E2,C2,F2] = ssdata(P2);
nX2 = size(A2,1); % will be 1
nW2 = size(E2,2); % will be 1

% Bode plot of both
bOpt = bodeoptions;
bOpt.PhaseVisible = 'off';
bOpt.MagUnits = 'abs';
bOpt.MagScale = 'log';
bOpt.FreqScale = 'linear';
bOpt.Xlim = [0 pi];
bOpt.Ylim = [1e-4 2];
bodeplot(P2,'r',P1,'k',bOpt)

% Form overall system which adds the outputs
A = blkdiag(A1,A2);
E = blkdiag(E1,E2);
C=[C1 C2];
F=[F1 F2];
nX = size(A,1);
nY = size(C,1);
nW = size(E,2);

% Noise variance and initial condition variance
% Keep it simple, and make everything Identity (appropriate
% dimension)
SigW = eye(nW);
Sxkk1 = eye(nX);

% Run several iterations to get the steady-state Kalman Gains
nIter = 40;
for i=1:nIter
    Swk = SigW;
    [~,Sxk1k,~,Sxkk,Sykk1,Lk,Hk,Gk,~] = ...
    KF231B([],Sxkk1,A,[],C,E,F,Swk,[],[]);
    Sxkk1 = Sxk1k;  
    saveLk(:,i) = Lk;
    if(i>2)
        if saveLk(:,i)-saveLk(:,i-1)<0.01
            break;
        end
    end
end

% Approximately 10 steps are needed for the convergence

%% Form Kalman filter with 3 outputs
AKF = A-Lk*C;
BKF = Lk;
CKF = [eye(nX);eye(nX)-Hk*C;-Gk*C];
DKF = [zeros(nX,nY);Hk;Gk];
SSKF = ss(AKF,BKF,CKF,DKF,TS);

%% Form matrix to extract estimate of y2_{k|k}
% We need [0 C2]*xhat_{k|k} + [0 F2]*what_{k|k}. Everything
% is scalar dimension, but we can form this matrix properly
% so that the example would work on other systems too.
M = [zeros(nY,nX) zeros(nY,nX1) C2 zeros(nY,nW1) F2];

%% Bode plot of filter
% It makes sense that the filter will try to "pass" some
% low frequencies, to preserve y2, but will cutoff
% high-frequencies to reject y1. The "pass" region should
% extend over the region where P2 has modest gain. The Bode
% magnitude plot confirms this
bodeplot(P2,'r',P1,'k',M*SSKF,bOpt)
legend('P2','P1','Filter');

%% Single Simulation
% Create a w sequence consistent with variance assumption
wSeq = randn(100,2);
%% Get y1 and y2 (separate simulations) for later comparison
y1 = lsim(P1,wSeq(:,1));
y2 = lsim(P2,wSeq(:,2));
y = y1 + y2;
%%
% Form the cascade (system output goes directly to Kalman
% Filter), and simulate, obtaining the outputs of Kalman
% Filter
Est = lsim(SSKF*ss(A,E,C,F,TS),wSeq);

%% Form Estimate of y2
% Est matrix is 100-by-6, so use transpose correctly to do
% reconstruction as a matrix multiply
y2Est = (M*Est')';

%% Plot
subplot(1,2,1);
plot(0:99,y2,'b+',0:99,y2Est,'ko',0:99,y,'r*');
legend('y2 (actual)','y2 (Est)','y (Measured)');
subplot(1,2,2);
plot(0:99,y2,'b+',0:99,y2Est,'ko');
legend('y2 (actual)','y2 (Est)');


%% time varying kalman filter
T = 20;

arrayA = zeros(size(A,1), size(A,2), T);
arrayE = zeros(size(E,1), size(E,2), T);
arrayC = zeros(size(C,1), size(C,2), T);
arrayF = zeros(size(F,1), size(F,2), T);

for i=1:T
    arrayA(:,:,i) = A;
    arrayE(:,:,i) = E;
    arrayC(:,:,i) = C;
    arrayF(:,:,i) = F;
end

arraySW = repmat(SigW,[1, 1, T]);
Sx0 = eye(nX);
m0 = 4;
Sxii1 = Sx0;
xii1 = m0;
% Create a specific initial condition and noise sequence
% Under ideal circumstances, this should be consistent with the statistical
% assumptions made in the previous code cell.  When studying robustness,
% namely how the filter performance degrades as assumptions are not met, it
% may be useful to create an initial condition and noise sequence which is
% not consistent with the assumptions
wSeq = randn(2,T);
x0 = 0;
% Simulate the system/KF one step at a time
emptyB = [];  % this template is for no control signal
emptyu = [];  % this template is for no control signal
y = zeros(nY,T);
xSeq = zeros(nX,T);
xSeq(:,1) = x0;
xEii1 = zeros(nX,T); xEii1(:,1) = m0;
xEii = zeros(nX,T); yEii = zeros(nY,T); y2Eii = zeros(1,T);
for i=0:T-1
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
   % Get estimates of x(i+1|i) using y(i) from KF
   [xi1i,Sxi1i,xii,Sxii,Syii1,Li,Hi,Gi,wii] = KF231B(xEii1(:,iMatlab),Sxii1,...
      Ai,emptyB,Ci,Ei,Fi,Swi,emptyu,y(:,iMatlab));
   xEii(:,iMatlab) = xii;
   xEii1(:,iMatlab+1) = xi1i;
   % Get x(i+1) from system model, using x(i) and w(i)
   xSeq(:,iMatlab+1) = Ai*xi + Ei*wi;
   % Shift the error-variance estimate so that when loop-index i advances,
   % the initial condition for this variance is correct.
   Sxii1 = Sxi1i;
   yEii(:,iMatlab) = Ci*xii;
   y2TVEii(iMatlab) = C2*xii(2) + F2*wii(2);
end

figure
plot(0:19,y2TVEii,'b+',0:19,y2Est(1:20),'ko');
