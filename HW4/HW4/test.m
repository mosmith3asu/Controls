%The system xdot = Ax + Bu + Mp + B2 w, p = delta*q, q = Nx + + D_12 u, y = Cx + D_22 u
%is quadratically stable under structured norm-bounded uncertainty if the following LMI in
%A,B,M,B2,N,D_12,C,D22 is feasible. The optimal controller state feedback u = Kx is formed with K = ZP^-1

clear all;
%Define given variables
lh=1;
ls=1;
m0=1;
k0=1;

nm=0.1;
nk=0.1;

% Define 9-matrix rep
A=[0 1; (k0*ls+m0*lh) 0];
B1 = [0 0 0 0; 
    -nm nm/(lh^2*m0) nk/(lh^2*m0) nk/(lh^2*m0)
    ];
B2 = [0; 
    -k0*ls/(lh^2*m0)];

C1=[(k0+m0*lh)/(lh^2*m0) 0; 
    lh*m0 0; 0 0; 
    ls*k0 0];
C2 = [1 0];

D11 = [-nm nm/(lh^2*m0) nk/(lh^2*m0) nk/(lh^2*m0); 0 0 0 0; 0 0 0 0; 0 0 0 0];
D12 = [-k0*ls/(lh^2*m0); 0; -ls*k0; 0];
D21 = [0 0 0 0];
D22 = [0];

B=[B1 B2];
M=[C1 D12 zeros(4,1)];
N=[B1; D21];
Q=[A B2; C2 D22];
A = D11;
% 
% A = [1 4 2;
%     0 1 1; 
%     3 8 5];
% B = [3; 
%     5.5; 
%     9];
% M = [6 0 2;
%     0 1 0; 10 4 3];
% B2 = [4; 0; 2];
% N = [0 0 4; 9 9 1; 3 5 4];
% D12 = [0; 1; 1];
% C = [1 2 0; 4 5 2; 9 8 1];
% D22 = [7; 9; 9];

eps = 0.0001;

P = sdpvar(size(A,1),size(A,2)); theta = sdpvar(size(M,1),size(M,1));
Z = sdpvar(size(B,2),size(B,1)); gamma = sdpvar(1,1);

Constraints = P >= eps*eye(size(P));
Mat = [A*P+B*Z+P*A'+Z'*B'+B2*B2'+M*theta*M' (C*P+D22*Z)' P*N'+Z'*D12'; C*P+D22*Z -gamma^2*eye(size(C*P)) zeros(size(P*N')); N*P+D12*Z zeros(size(N*P)) -theta];
Constraints = [Constraints, Mat <= eps*eye(size(Mat))];

sol = optimize(Constraints,gamma);

if (sol.problem ~= 0)
    'Problem Infeasible: Not Quadratically Stable'
else
    'Problem Feasible! Quadratically Stable!'
end

Z = value(Z);
P = value(P);
K = Z*inv(P)