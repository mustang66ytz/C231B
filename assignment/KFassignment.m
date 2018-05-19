%% MEC231B Kalman 1: Taozheng Yang
%% Q1
% histogram of data uniformly distributed in the range of -3 to 2
clear all;
clc;

N = 500;
X = linspace(-3,2,N);
figure;
subplot(2,1,1); hist(X);
legend("histogram of unpermutated data");
I = randperm(N);
Y = X(I);
isequal(X,Y)
subplot(2,1,2); hist(Y);
legend("histogram of the permutated data");
hold on;

% histogram of data uniformly distributed in the range of -1 to 1 and its
% cosine distribution
figure;
N = 500;
X = linspace(-1,1,N);
subplot(2,1,1); hist(X);
xlabel('Angle, rads'); ylabel('Counts')
subplot(2,1,2); hist(cos(X));
xlabel('COS(Angle, rads)'); ylabel('Counts')
hold on;

% histogram of data logrithmically uniformly distributed in the range of -1
% to 1
figure;
N = 500;
X = logspace(-1,1,N);
subplot(2,1,1); hist(X);
I = randperm(N);
Y = X(I);
isequal(X,Y)
subplot(2,1,2); hist(Y);
hist(X,100);

%% Q6
n = 3;
M = 4;
A = [1 2 1 0; 2 3 2 1; 3 4 3 2];
colmeanA = mean(A, 1);
rowmeanA = mean(A, 2);
B = A - repmat(mean(A,2),[1 size(A,2)]);
% B is constucted with the items (ith row and jth column), which are the 
% differences between Aij and rowmeanA1j (B_ij = A_ij - rowmeanA_1j)
C = A - repmat(mean(A,1),[size(A,1) 1]);
% C is constucted with the items (ith row and jth column), which are the 
% differences between Aij and colmeanAi1 (C_ij = A_ij - colmeanA_i1)

%% Q7
% eigenvalues of symmetric matrix are real numbers, but are not real
% numbers for unsymmetric matrices.
clear;
clc;

A = rand(7,7);
isequal(A, A')
eig(A)
Asym = A + A';  % make symmetric by adding transpose
isequal(Asym, Asym')
eig(Asym)

%% 
%square, randomly generated matrices are almost always invertible
A = rand(5,5);
shouldBeIdentity = A*inv(A)
norm(shouldBeIdentity-eye(5))

%% 
%if W ? Rn×n is invertible, then W'W ? 0
W = rand(5,5);
eig(W'*W)
% all eigenvalues are greater than zero, implies PD.

%% 
%sqrtm computes the unique positive-semidefinite symmetric matrix square 
%root of a positive-semidefinite matrix
W = rand(5,5);
M = W'*W;
isequal(M, M')
S = sqrtm(M);
isequal(S, S')
eig(S)
S*S - M

%% Q8
%uniformly generated
A = rand(1,500);
figure;
hist(A)
legend("distribution of numbers generated randomly between zero and one");

% normally generated
A = randn(1,500);
figure;
hist(A);
legend("distribution of numbers generated randomly between zero and one differently");


%% 8(C)
n = 3;
M = 4;
Xmat1 = rand(n, M);
Xmat2 = randn(n, M);

% the mean of items in Xmat1 should be around 0.5
% the variance of Xmat1 should be around 1/12
mean(mean(Xmat1))
mean(var(Xmat1))

% the mean of items in Xmat2 should be around 0
% the variance of Xmat2 should be around 1
mean(mean(Xmat2))
mean(var(Xmat2))

Ymat = sqrt(12)*(rand(n, M)-0.5)
% the mean of items in Ymat should be around 0
colmean3 = mean(Ymat, 1);
mean3 = mean(colmean3, 2)
% the variance of Ymat should be around 1
Ymatreshape = reshape(Ymat, [1,12]);
var(Ymatreshape)

%% 8(d)
clear;
clc;

M = 500;
n=1;
Xmat = sqrt(12)*(rand(n, M)-0.5); 
Ymat = randn(n, M);
% Verify the two means and variances are very similiar
meanXmat = mean(Xmat);
meanYmat = mean(Ymat);
varXmat = var(Xmat);
varYmat = var(Ymat);

%% 8(e)
M = 500;
n = 3;
Xmat = sqrt(12)*(rand(n, M)-0.5);
mean(mean(Xmat))
mean(var(Xmat))
%% 8(f)
n = 3;
M = 1000;
X = randn(n, M);
S = [3, -1, 1; -1, 2, 0.5; 1, 0.5, 1];
A = sqrtm(S);
Y = A*X;

%% 8(g)
n = 4;
M = 1000;
Xmat = randn(n, M);
muX = mean(Xmat,2)
D = Xmat - repmat(muX,[1 M]);
D*D'/M
I = randperm(n*M);
Ymat = reshape(Xmat(I),[n M]);
muY = mean(Ymat,2)
D = Ymat - repmat(muY,[1 M]);
D*D'/M
% verify the rand case
n = 4;
M = 1000;
Xmat = rand(n, M);
muX = mean(Xmat,2)
D = Xmat - repmat(muX,[1 M]);
D*D'/M
I = randperm(n*M);
Ymat = reshape(Xmat(I),[n M]);
muY = mean(Ymat,2)
D = Ymat - repmat(muY,[1 M]);
D*D'/M


%% Q9
% Generating random variables for initial conditions and disturbances
M = 5000; % outcomes in probability model
nX = 5;
nD = 3;
duration = 100; % length of disturbance sequence
R = randn(nX+nD*duration, M);

x0vec = randn(nX, 1);
dSeq = randn(nD, duration);

%% 9(a)
A = [0.8,-0.6,0.15;1.5,0.02,1.3;-0.26,-0.16,-1];
E = [1,0;0,0;0,0];
C = [1,-1,0;0,1,-1];
F = [0,1;0,-1];

% Initial Conditions
figure;
% Eight initial Conditions
for i = 1:8
    % Randomly geneate initial conditions for eight times
    X = randn(3, 1);
    % S is the variance of x
    S = [3, -1, 1; -1, 2, 0.5; 1, 0.5, 1];
    % define simulation duration
    k = 30;
    % the expectionof initial state is [1;1;1]
    x0 = sqrtm(S)*X + [1;1;1];
    % initialize the state sequence
    x = zeros(3,k);
    x(:,1) = x0;
    % initialize the output sequence
    y = zeros(2,k);
    
    % simulate the system from k=0 to k = 30
    for j = 1:k+1
        % randomly generate disturbance each iteration time
        Sd = [0.5, 0.5; 0.5, 1];
        d = sqrtm(Sd)*randn(2,1);
        % forward pass the system
        y(:,j) = C*x(:,j)+F*d;
        x(:,j+1) = A*x(:,j)+E*d;
    end
duration = linspace(1,k+1,k+1);
% Plot the eight outcomes
plot(duration,y(2,:));
hold on
end
hold off