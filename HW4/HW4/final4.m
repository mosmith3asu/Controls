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


% Define 9-matrix rep%%%%%%%%%%%%%%%%%%%%%
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

B = [B1 B2];
C = [C1; C2];
D = [D11 D12; D21 D22];
% Define 9-matrix rep%%%%%%%%%%%%%%%%%%%%%
% A=[0 1; (k0*ls+m0*lh) 0];
% B1 = [0 0 0 ; 
%     -nm nm/(lh^2*m0) nk/(lh^2*m0)];
% B2 = [0; -k0*ls/(lh^2*m0)];
% 
% C1=[(k0+m0*lh)/(lh^2*m0) 0; lh*m0 0; 0 0];
% C2 = [1 0];
% 
% D11 = [-nm nm/(lh^2*m0) nk/(lh^2*m0); 0 0 0; 0 0 0];
% D12 = [-k0*ls/(lh^2*m0); 0; -ls*k0];
% D21 = [0 0 0];
% D22 = [0];
B = [B1 B2];
C = [C1; C2];
D = [D11 D12; D21 D22];
% P11=D11;
% P12=[C1 D12];
% P21=[B1; D21];
% P22=[A B2; C2 D22];
% %P11=D11;
M=[C1 D12];
N=[B1; D21];
Q=[A B2; C2 D22];


%%

% P = sdpvar(size(A,1),size(A,2)); theta = sdpvar(size(M,1),size(M,1));
% Z = sdpvar(size(B,2),size(B,1)); gamma = sdpvar(1,1);
% 
% Constraints = P >= eps*eye(size(P));
% Mat = [
%     A*P+B*Z+P*A'+Z'*B'+B2*B2'+M*theta*M' (C*P+D22*Z)' P*N'+Z'*D12'; 
%     C*P+D22*Z -gamma^2*eye(size(C*P)) zeros(size(P*N'));
%     N*P+D12*Z zeros(size(N*P)) -theta];
% 
% Constraints = [Constraints, Mat <= eps*eye(size(Mat))];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
P = sdpvar(size(A,1),size(A,2)); 
Z = sdpvar(size(B,2),size(B,1)); 
gamma = sdpvar(1,1);
theta = blkdiag(sdpvar(1),sdpvar(1),sdpvar(1));

Constraints = P >= eps*eye(size(P));
Constraints = [Constraints, theta >= 0];
M11= A*P+B*Z+P*A'+Z'*B'+B2*B2'+M*theta*M';
M12 = (C*P+D22*Z)';
M13 = P*N'+Z'*D12';
M21 = C*P+D22*Z;
M22 = -gamma^2*eye(size(C*P));
M23 = zeros(size(P*N'));
M31 = N*P+D12*Z;
M32 =  zeros(size(N*P));
M33 = -theta;

Mat = [
    M11 M12 M13;
    M21 M22 M23;
    M31 M32 M33];

Constraints = [Constraints, Mat <= eps*eye(size(Mat))];

sol = optimize(Constraints,gamma);

Z = value(Z);
P = value(P);
K = Z*inv(P)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%










