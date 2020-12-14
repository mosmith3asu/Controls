% Demo from M. Peet for lecture 10 in MAE 598
% 
% % Model 
clear all
clc

A=[
    1   1   1   0   0   1;
    -1  0   -1  0   0   1;
    1   0   0   -1  1   1;
    -1  1   -1  0   0   0;
    -1  -1  1   1   -1  -1;
    0   -1  0   0   -1  0
    ];
B1 = [
    0   -1  -1;
    0   0   0;
    -1  -1  1;
    -1  0   1;
    0   0   0;
    -1  1   1
    ];
B2 = [
    0   0   0;
    -1  0   1;
    -1  1   0;
    1   -1  0;
    -1  0   -1;
    0   1   1
    ];
C1 = [
    0   1   0   -1  -1  -1;
    0   0   0   -1  0   0;
    1   0   0   0   -1  0
    ];
C2 = zeros(size(C1));
D12 = [
    0   1   1 ;
    0   0   0;
    1   1   1
    ];
D11 = [
    0    0   1;
    -1  0   0;
    0   0   0
    ];
D21 = zeros(size(D11));
D22 = zeros(size(D11));

% measure numbers of inputs and outputs 

eta=.0001;    % degree of strict positivity   
ns=size(A,1);   % number of states
nc=size(B2,2);  % number of actuators
nd=size(B1,2);  % number of external inputs
nr=size(C1,1);  % number of regulated outputs
C2t=eye(ns); D21t=zeros(ns,nd); D22t=zeros(ns,nc);
nm=size(C2,1);  % number of sensors


% H-infinity State Feedback Controller Synthesis


% Declare the variables
gamma=sdpvar(1);               % represents the bound on the H-infinity norm of the CL system.
X=sdpvar(ns);
Z=sdpvar(nc,ns,'full');
W=sdpvar(nr);

% declare constraints
MAT=[A*X+X*A'+B2*Z+Z'*B2'       B1              (C1*X+D12*Z)';
     B1'                        -gamma*eye(nd)   D11';
     C1*X+D12*Z                 D11              -gamma*eye(nr)];
F=[MAT<=0];
F=[F;X>=eta*eye(ns)];

OPTIONS = sdpsettings('solver','mosek','verbose',0);

% Solve the LMI, minimizing gamma
optimize(F,gamma,OPTIONS);
gamman=value(gamma)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% retrieve decision variables
Xn=value(X);
Zn=value(Z);
K=Zn*inv(Xn)

controller=ss(K);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Close the loop with Lower LFT
plant=ss(A,[B1 B2],[C1;C2t],[D11 D12; D21t D22t]);
sys_cl=lft(plant,controller);
%[Acl,Bcl,Ccl,Dcl] = ssdata(sys_cl);
%sys_cl=ss(Acl,Bcl,Ccl,0);
Hinf_Norm = norm(sys_cl,inf)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% compare with Matlab built-in functions
C2tt=ones(1,ns); D21tt=ones(1,nd);D22tt=zeros(1,nc);
sys=ss(A,[B1 B2],[C1;C2tt],[D11 D12; D21tt D22tt]);nm=size(C2tt,1);  % number of sensors
[Knew,CL,GAM,INFO]=hinfsyn(sys,nm,nc);
Matlab_Norm =  INFO.GAMFI % This assumes controller depnds on disturbances too, which is apples to oranges

%% 1B
disp("####################################################################")
disp("1B")
disp("####################################################################")

clear all
A=[
    1   1   1   0   0   1;
    -1  0   -1  0   0   1;
    1   0   0   -1  1   1;
    -1  1   -1  0   0   0;
    -1  -1  1   1   -1  -1;
    0   -1  0   0   -1  0
    ];
B1 = [
    0   -1  -1;
    0   0   0;
    -1  -1  1;
    -1  0   1;
    0   0   0;
    -1  1   1
    ];
B2 = [
    0   0   0;
    -1  0   1;
    -1  1   0;
    1   -1  0;
    -1  0   -1;
    0   1   1
    ];
C1 = [
    0   1   0   -1  -1  -1;
    0   0   0   -1  0   0;
    1   0   0   0   -1  0
    ];
C2 = zeros(size(C1));
D12 = [
    0   1   1 ;
    0   0   0;
    1   1   1
    ];
D11 = [
    0    0   1;
    -1  0   0;
    0   0   0
    ];
D21 = zeros(size(D11));
D22 = zeros(size(D11));

% measure numbers of inputs and outputs 

eta=.0001;    % degree of strict positivity   
ns=size(A,1);   % number of states
nc=size(B2,2);  % number of actuators
nd=size(B1,2);  % number of external inputs
nr=size(C1,1);  % number of regulated outputs
C2t=eye(ns); D21t=zeros(ns,nd); D22t=zeros(ns,nc);
nm=size(C2,1);  % number of sensors


% H-infinity State Feedback Controller Synthesis


% Declare the variables
gamma=sdpvar(1);               % represents the bound on the H-infinity norm of the CL system.
X= blkdiag(sdpvar(2,2),sdpvar(4,4));

Z= blkdiag(sdpvar(1,2),sdpvar(2,4));
%Z=sdpvar(nc,ns,'full');
W=sdpvar(nr);

% declare constraints
MAT=[A*X+X*A'+B2*Z+Z'*B2'       B1              (C1*X+D12*Z)';
     B1'                        -gamma*eye(nd)   D11';
     C1*X+D12*Z                 D11              -gamma*eye(nr)];
F=[MAT<=0];
F=[F;X>=eta*eye(size(ns))];

OPTIONS = sdpsettings('solver','mosek','verbose',0);

% Solve the LMI, minimizing gamma
optimize(F,gamma,OPTIONS);
gamman=value(gamma)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% retrieve decision variables
Xn=value(X); 
Zn=value(Z); 
K=Zn*inv(Xn)
controller=ss(K);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Close the loop with Lower LFT
plant=ss(A,[B1 B2],[C1;C2t],[D11 D12; D21t D22t]);
sys_cl=lft(plant,controller);
%[Acl,Bcl,Ccl,Dcl] = ssdata(sys_cl);
%sys_cl=ss(Acl,Bcl,Ccl,0);
Hinf_norm = norm(sys_cl,inf)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% compare with Matlab built-in functions
C2tt=ones(1,ns); D21tt=ones(1,nd);D22tt=zeros(1,nc);
sys=ss(A,[B1 B2],[C1;C2tt],[D11 D12; D21tt D22tt]);nm=size(C2tt,1);  % number of sensors
[Knew,CL,GAM,INFO]=hinfsyn(sys,nm,nc);
Matlab_norm = INFO.GAMFI % This assumes controller depnds on disturbances too, which is apples to oranges
%INFO.KFI-K % compare full-state gains.


