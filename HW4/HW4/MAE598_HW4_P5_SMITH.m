%% Problem 3
% Definition of variables
clear all
clc


%% P3.a
disp("##########################")
disp("##########5.a#############")
disp("##########################")
%Define given variables
m0=1.0;
c0=0.1;
k0=1.0;
nm=0.1;
nc=1.0;
nk=0.2;

% Define 9-matrix rep
A=[0 1 ;-k0/m0 -c0/m0];

B1 = [0 0 0; -nm -nc/m0 -nk/m0];
B2 = [0;1/m0];

C1=[-k0/m0 -c0/m0; 
    0 c0; 
    k0 0];
C2 = [1 0];
D11 = [-nm -nc/m0 -nk/m0; 0 0 0; 0 0 0];
D12 = [1/m0; 0; 0];
D21 = [0 0 0];
D22 = [0];

% Define 9-matrix rep
% A=[0 1/m0; -k0 -c0];
% B1 = [-nm 0 0; 0 nk nc];
% B2 = [0; 1];
% 
% C1=[0 1/m0; -k0 0; 0 -c0];
% C2 = [1 0];
% 
% D11 = [-nm 0 0; 0 0 0; 0 0 0];
% D12 = [0; 0; 0];
% D21 = [0 0 0];
% D22 = [0];


% Define SS L12.11
A
B=[B1 B2]
C=[C1;C2]
D=[D11 D12; D21 D22]

eps = 0.0001;

X = sdpvar(size(A,1),size(A,2));
%gamma = sdpvar(1,1);
gamma = 1.5; % manually increased and decreased to reach feasability
theta = blkdiag(sdpvar(1),sdpvar(1),sdpvar(1),sdpvar(1));

Mat = [A'*X+X*A X*B; B'*X -theta]+ (1/gamma^2)*[C'; D']*theta*[C D];
gammarange = 1;
Constraints = [
    %theta >= eps*eye(size(theta))
    X >= eps*eye(size(X));
    Mat <= eps*eye(size(Mat))];
options = sdpsettings('solver','mosek','verbose',0);
%Objective = [];
Objective = theta;
optimize(Constraints,Objective,options)
%bis = bisection(Constraints,Objective,options)
%gammaf = value(gamma)
Xf = value(X)
sympref('FloatingPointOutput',true);
Latex_out = latex(sym(Xf))
thetaf = value(theta)
sympref('FloatingPointOutput',true);
Latex_out = latex(sym(thetaf))
%% 3.b
disp("##########################")
disp("##########5.b#############")
disp("##########################")
clear all
%Define given variables
m0=1.0;
c0=0.1;
k0=1.0;
nm=0.1;
nc=1.0;
nk=0.2;

A=[0 1; -k0/m0 -c0/m0];
B=[0 0 0; -nm -nc/m0 -nk/m0]; %B1
M=[0; 1/m0]; %B2

C = [-k0/m0 -c0/m0; 0 c0; k0 0]; %C1
N = [1 0]; % C2

D11 = [-nm -nc/m0 -nk/m0; 0 0 0; 0 0 0];
D12 = [1/m0; 0; 0];
D21 = [0 0 0];
Q= [0]; % D22


eps = 0.0001;

P = sdpvar(size(A,1),size(A,2));
Z = sdpvar(size(B,2),size(B,1));
%theta = sdpvar(size(M,2),size(M,2));
theta = blkdiag(sdpvar(1),sdpvar(1),sdpvar(1),sdpvar(1));

Constraints = P >= eps*eye(size(P));

Mat = [
    A*P+B*Z+P*A'+Z'*B' P*N'+Z'*D12';
    N*P+D12*Z zeros(size(P*N'))];
Mat = Mat + [
    M*theta*M' M*theta*Q'; 
    Q*theta*M' Q*theta*Q'-theta];
Constraints = [Constraints, Mat <= eps*eye(size(Mat))];

sol = optimize(Constraints,[],options);

Z = value(Z);
P = value(P);
K = Z*inv(P)
%% 3.c
disp("##########################")
disp("##########5.c#############")
disp("##########################")
%% 3.d
disp("##########################")
disp("##########5.c#############")
disp("##########################")
P = sdpvar(size(A,1),size(A,2)); 
theta = blkdiag(sdpvar(1),sdpvar(1),sdpvar(1),sdpvar(1));
Z = sdpvar(size(B,2),size(B,1)); 
gamma = sdpvar(1,1);
%theta = size(P*N')
Constraints = P >= eps*eye(size(P));
Constraints = [Constraints, theta >= 0];
Mat = [
    A*P+B*Z+P*A'+Z'*B'+B2*B2'+M*theta*M' (C*P+D22*Z)' P*N'+Z'*D12'; 
    C*P+D22*Z -gamma^2*eye(size(C*P)) zeros(size(P*N')); 
    N*P+D12*Z zeros(size(N*P)) -theta];
Constraints = [Constraints, Mat <= eps*eye(size(Mat))];

sol = optimize(Constraints,gamma);

Z = value(Z);
P = value(P);
K = Z*inv(P)
%% 3.e
disp("##########################")
disp("##########5.c#############")
disp("##########################")

