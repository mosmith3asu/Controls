%enforcing quadratic stability using the polytopic framework%% Problem 3
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

B = [B1 B2];
C = [C1; C2];
D = [D11 D12; D21 D22];

% P11=D11;
% P12=[C1 D12];
% P21=[B1; D21];
% P22=[A B2; C2 D22];
% %P11=D11;
% M=[C1 D12];
% N=[B1; D21];
% Q=[A B2; C2 D22];

%A=P11;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% P = sdpvar(size(A,1),size(A,2)); 
% theta = blkdiag(sdpvar(1),sdpvar(1),sdpvar(1),sdpvar(1));
% Z = sdpvar(size(B,2),size(B,1)); 
% gamma = sdpvar(1,1);
% %theta = size(P*N')
% Constraints = P >= eps*eye(size(P));
% Constraints = [Constraints, theta >= 0];
% Mat = [
%     A*P+B*Z+P*A'+Z'*B'+B2*B2'+M*theta*M' (C*P+D22*Z)' P*N'+Z'*D12'; 
%     C*P+D22*Z -gamma^2*eye(size(C*P)) zeros(size(P*N')); 
%     N*P+D12*Z zeros(size(N*P)) -theta];
% Constraints = [Constraints, Mat <= eps*eye(size(Mat))];
% 
% sol = optimize(Constraints,gamma);
% 
% Z = value(Z);
% P = value(P);
% K = Z*inv(P)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% P = sdpvar(size(A,1),size(A,2)); 
% Z = sdpvar(size(B,2),size(B,1)); 
% mu = sdpvar(1)
% 
% Mat = [
%     A*P+B*Z+P*A'+Z'*B' P*N'+Z'+D12';
%     N*P+D12*Z zeros(size(P*N'))]+mu*[M*M' M*Q'; Q*M' Q*Q'-eye(size(Q*Q'))];
% 
% Constraints = [Constraints, mu >= 0];
% Constraints = [Constraints, P >= eps*eye(size(Mat))];
% Constraints = [Constraints, Mat <= eps*eye(size(Mat))];
% 
% sol = optimize(Constraints,gamma);
% 
% Z = value(Z);
% P = value(P);
% K = Z*inv(P)
% %S_us = P22+P21*inv(eye(size(P11))-Delta*P11)*Delta*P12
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% a1_prime = 0.2;
% a1=[1.593,0.398];
% a2=[260.559,51.003];
% a3=[185.488,53.84];
% a4=[1.506,0.421];
% a5=[0.298,0.078];
% 
% % form 2^k aggregate combination of all verticies
%  pre_a_aggr=combvec(a1,a2,a3,a4,a5)';
%  a_aggr=zeros(size(pre_a_aggr));
%  n=5; % number of inputs
%  for i=1:n
%      a_aggr(:,i)=pre_a_aggr(:,n-i+1);
%  end
% a_aggr = a_aggr';
% 
% sympref('FloatingPointOutput',true);
% Latex_out = latex(sym(a_aggr))
% 
% % Extract new verticies
% a1 = a_aggr(1,:);
% a2 = a_aggr(2,:);
% a3 = a_aggr(3,:);
% a4 = a_aggr(4,:);
% a5 = a_aggr(5,:);
% 
% 
% N=size(a1,2);
% for i=1:1:N
% A=[-a4(i) 1 -a5(i);((-a1_prime*a4(i))-a2(i)) (a1_prime-a1(i)) ((-a1_prime*a5(i))-a3(i));0 0 -(1/tz)];
% B1 = (wx/57.3)*[-1 0;-a1_prime JxJy_Jz; 0 0];
% B2 = [0;0;(1/tz)]; 
% C = (1/(57.3*g))*[(57.3*g) 0 0;V*a4(i) 0 V*a5(i)];
% D12 = [0;0]; 
% D11 = (1/(57.3*g))*[0 0; V*b7 0];
% 
% 
% M11 = Y*A' + A*Y + Z'*B2' + B2*Z;
% M21 = B1';
% M31 = C*Y + D12*Z;
% M32 = D11;
% M22 = -gamma*eye(2);
% M33 = -gamma*eye(2);
% M = [M11 M21' M31';
%     M21 M22 M32';
%     M31 M32 M33];
% F = [F; M <= -eps*eye(size(M))];
% end
% 
% 
% F=[F;Y>=eps*eye(size(Y))];
% optimize(F,gamma,options);
% 
% gamma = value(gamma)
% 
% K = value(Z)*(value(Y'))
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
B1 = [0 0 0 ; 
    -nm nm/(lh^2*m0) nk/(lh^2*m0)];
B2 = [0; -k0*ls/(lh^2*m0)];

C1=[(k0+m0*lh)/(lh^2*m0) 0; lh*m0 0; 0 0];
C2 = [1 0];

D11 = [-nm nm/(lh^2*m0) nk/(lh^2*m0); 0 0 0; 0 0 0];
D12 = [-k0*ls/(lh^2*m0); 0; -ls*k0];
D21 = [0 0 0];
D22 = [0];
B = [B1 B2];
C = [C1; C2];
D = [D11 D12; D21 D22];

eps = 0.0001;

X = sdpvar(size(A,1),size(A,2));
%gamma = sdpvar(1,1);
gamma = 0.0001; % manually increased and decreased to reach feasability
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











