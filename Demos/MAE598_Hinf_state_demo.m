% Demo from M. Peet for lecture 10 in MAE 598
% 
% % Model 
clear all
t1=.66;t2=2.25;t3=.55;t4=3;t5=.94;t6=.64;
g11=1.32;g22=1;g32=1;g33=1;g35=1;g36=1;g57=1;g64=1;g68=1;
b21=.9;b31=.05;b42=.9;b43=.5;b54=.67;b34=.18;b25=.5;b14=.65;b46=.01;b45=.1;
A=[-1/t1    0       0      b14/t1  0       0;
    b21/t2  -1/t2   0      0       b25/t2  0;
    b31/t3  0       -1/t3  b34/t3  0       0;
    0       b42/t4  b43/t4 -1/t4   b45/t4  b46/t4;
    0       0       0      b54/t5  -1/t5   0;
    0       0       0      0       0       -1/t6];
B2=[g11 0;0 0; 0 0; 0 0; 0 0; 0 g68/t6]
C1=[0 0 0 0 (1-b25-b45) 0; 0 0 0 0 0 0; 0 0 0 0 0 0];
C2=[0 0 1 0 0 0;0 0 0 1 0 0];
B1=[g11 0   0   0   0;
    0   g22 0   0   0 ;
    0   g32 g33 0   0;
    0   0   0   0   0;
    0   0   0   g57 0;
    0   0   0   0   g68;]
D11=zeros(3,5);D12=[0 0;1 0;0 1];D21=zeros(2,5);D22=zeros(2,2);
% clear all 
% delta=.00001;
%  A=[0 0;1 -1]+delta*eye(2);B1=[0;1];B2=[1;0];C1=[1 0;.5 -1];C2=[0 1];D11=[0;0];D12=[1;0];D21=[1];D22=[0];
% 

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

OPTIONS = sdpsettings('solver','sedumi');

% Solve the LMI, minimizing gamma
optimize(F,gamma,OPTIONS);
gamman=value(gamma)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% retrieve decision variables
Xn=value(X); Zn=value(Z); 
K=Zn*inv(Xn);
controller=ss(K);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Close the loop with Lower LFT
plant=ss(A,[B1 B2],[C1;C2t],[D11 D12; D21t D22t]);
sys_cl=lft(plant,controller);
%[Acl,Bcl,Ccl,Dcl] = ssdata(sys_cl);
%sys_cl=ss(Acl,Bcl,Ccl,0);
norm(sys_cl,inf)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% compare with Matlab built-in functions
C2tt=ones(1,ns); D21tt=ones(1,nd);D22tt=zeros(1,nc);
sys=ss(A,[B1 B2],[C1;C2tt],[D11 D12; D21tt D22tt]);nm=size(C2tt,1);  % number of sensors
[Knew,CL,GAM,INFO]=hinfsyn(sys,nm,nc);
INFO.GAMFI % This assumes controller depnds on disturbances too, which is apples to oranges
%INFO.KFI-K % compare full-state gains.
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % compare with Matlab built-in functions
% C2=[0 0 1 0 0 0;0 0 0 1 0 0];
% D21=zeros(2,5);D21(1,1)=1;D21(2,2)=1;D22=zeros(2,2);
% sys=ss(A,[B1 B2],[C1;C2],[D11 D12; D21 D22]);nm=size(C2,1);  % number of sensors
% [Knew,CL,GAM,INFO]=h2syn(sys,nm,nc);
% norm(INFO.GFI)
% INFO.KFI-K % compare full-state gains.
