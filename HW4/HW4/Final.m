%% Problem 3
% Definition of variables
clear all
clc
options = sdpsettings('solver','mosek','verbose',0);


%% P3.a
disp("##########################")
disp("##########5.a#############")
disp("##########################")
%Define given variables
lh=1;
ls=1;
m0=1;
k0=1;

nm=0.1;
nk=0.1;

% Define 9-matrix rep
% A=[0 1 ;-k0/m0 -c0/m0];
% 
% B1 = [0 0 0; -nm -nc/m0 -nk/m0];
% B2 = [0;1/m0];
% 
% C1=[-k0/m0 -c0/m0; 0 c0; k0 0];
% C2 = [1 0];
% D11 = [-nm -nc/m0 -nk/m0; 0 0 0; 0 0 0];
% D12 = [1/m0; 0; 0];
% D21 = [0 0 0];
% D22 = [0];

% Define 9-matrix rep
A=[0 1; (k0*ls+m0*lh) 0];
B1 = [0 0 0 0; 
    -nm nm/(lh^2*m0) nk/(lh^2*m0) nk/(lh^2*m0)
    ];
B2 = [0; -k0*ls/(lh^2*m0)];

C1=[(k0+m0*lh)/(lh^2*m0) 0; lh*m0 0; 0 0; ls*k0 0];
C2 = [1 0];

D11 = [-nm nm/(lh^2*m0) nk/(lh^2*m0) nk/(lh^2*m0); 0 0 0 0; 0 0 0 0; 0 0 0 0];
D12 = [-k0*ls/(lh^2*m0); 0; -ls*k0; 0];
D21 = [0 0 0 0];
D22 = [0];

thetak = sdpvar(1);
thetam= sdpvar(1);
Delta = blkdiag(thetam, thetam,thetak,thetak);

B = [B1 B2];
C = [C1; C2];
D = [D11 D12; D21 D22];

P11=D11;
P12=[C1 D12];
P21=[B1; D21];
P22=[A B2; C2 D22];

A=P11;
M=P12;
N=P21;
Q=P22;

% Find
P = sdpvar(size(A,1),size(A,2)); 
Constraints = [P >= eps*eye(size(P))];
Z = sdpvar(size(B,2),size(B,1)); 

% and
theta = blkdiag(thetam, thetam,thetak,thetak);

%theta = size(P*N')

%Constraints = [Constraints, theta >= 0];
Mat = [
    A*P+B*Z+P*A'+Z'*B' P*N'+Z'+D12';
    N*P+D12*Z zeros(size(P*N'))]+theta*[M*M' M*Q'; Q*M' Q*Q'-eye(size(Q*Q'))];
Constraints = [Constraints, Mat <= eps*eye(size(Mat))];

sol = optimize(Constraints,gamma);

Z = value(Z);
P = value(P);
K = Z*inv(P)
%S_us = P22+P21*inv(eye(size(P11))-Delta*P11)*Delta*P12
