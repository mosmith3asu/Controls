% Demo from M. Peet for lecture 10 in MAE 598
% 
% % Model 
% clear all
% t1=.66;t2=2.25;t3=.55;t4=3;t5=.94;t6=.64;
% g11=1.32;g22=1;g32=1;g33=1;g35=1;g36=1;g57=1;g64=1;g68=1;
% b21=.9;b31=.05;b42=.9;b43=.5;b54=.67;b34=.18;b25=.5;b14=.65;b46=.01;b45=.1;
% A=[-1/t1    0       0      b14/t1  0       0;
%     b21/t2  -1/t2   0      0       b25/t2  0;
%     b31/t3  0       -1/t3  b34/t3  0       0;
%     0       b42/t4  b43/t4 -1/t4   b45/t4  b46/t4;
%     0       0       0      b54/t5  -1/t5   0;
%     0       0       0      0       0       -1/t6];
% B2=[g11 0;0 0; 0 0; 0 0; 0 0; 0 g68/t6]
% C2=[0 0 1 0 0 0;0 0 0 1 0 0];
% C1=[0 0 0 0 (1-b25-b45) 0; 0 0 0 0 0 0; 0 0 0 0 0 0];
% B1=[g11 0   0   0   0;
%     0   g22 0   0   0 ;
%     0   g32 g33 0   0;
%     0   0   0   0   0;
%     0   0   0   g57 0;
%     0   0   0   0   g68;]
% D11=zeros(3,5);D12=[0 0;1 0;0 1];D21=zeros(2,5);D22=zeros(2,2);
clear all 
delta=.00001
 A=[0 0;1 -1]+delta*eye(2);B1=[0;1];B2=[1;0];C1=[1 0;.5 -1];C2=[0 1];D11=[0;0];D12=[1;0];D21=[1];D22=[0];


% measure numbers of inputs and outputs 

eta=.0001;    % degree of strict positivity   
ns=size(A,1);   % number of states
nc=size(B2,2);  % number of actuators
nm=size(C2,1);  % number of sensors
nd=size(B1,2);  % number of external inputs
no=size(C1,1);  % number of regulated outputs


% H-infinity Dynamic Output Feedback Controller Synthesis


% Declare the variables
gamma=sdpvar(1);               % represents the bound on the H-infinity norm of the CL system.
X1=sdpvar(ns);
Y1=sdpvar(ns);
Z=sdpvar(no);
An=sdpvar(ns,ns,'full')
Cn=sdpvar(nc,ns,'full');
Dn=sdpvar(nc,nm,'full');
Bn=sdpvar(ns,nm,'full');

% declare constraints
%F=[[X1 eye(ns); eye(ns) Y1]>=eps*eye(2*ns)]
%F=[F;Y1>=eps*eye(ns)]
MAT=[A*X1+X1*A'+B2*Cn+Cn'*B2'  (A'+An+(B2*Dn*C2)')'        B1+B2*Dn*D21           ; 
     A'+An+(B2*Dn*C2)'         Y1*A+A'*Y1+Bn*C2+C2'*Bn'    Y1*B1+Bn*D21           ;
     (B1+B2*Dn*D21)'           (Y1*B1+Bn*D21)'             -eye(nd)         ];
MAT2=[  X1                      eye(ns)                    (C1*X1+D12*Cn)';
        eye(ns)                 Y1                          (C1+D12*Dn*C2)';
        C1*X1+D12*Cn            C1+D12*Dn*C2                Z ];

    MAT3=[A*X1+X1*A'+B2*Cn+Cn'*B2'  (A'+An+(B2*Dn*C2)')'; 
     A'+An+(B2*Dn*C2)'         Y1*A+A'*Y1+Bn*C2+C2'*Bn'    ];
 
    F=[MAT<=0];
F=[F;MAT2>=eta*eye(2*ns+no)];
%    F=[F; MAT3<=-.1*eye(2*ns)];
F=[F;(D11+D12*Dn*D21)==0];
F=[F;trace(Z)<=gamma];

OPTIONS = sdpsettings('solver','sedumi');

% Solve the LMI, minimizing gamma
optimize(F,gamma,OPTIONS);
gamman=sqrt(value(gamma))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% retrieve decision variables
X1n=value(Y1); Y1n=value(X1); Ann=value(An);Bnn=value(Bn);Cnn=value(Cn);Dnn=value(Dn);
temp1=[Ann Bnn; Cnn Dnn]-[X1n*A*Y1n zeros(ns,nm); zeros(nc,ns) zeros(nc,nm)]

% Choose X2, Y2, so that X2*Y2=I-X1*Y1;
Y2n=eye(ns);X2n=eye(ns)-X1n*Y1n;

% Reverse variable substitution
temp2=inv([X2n X1n*B2;zeros(nc,ns) eye(nc)])*temp1*inv([Y2n' zeros(ns,nm); C2*Y1n eye(nm)]);
Ak2=temp2(1:ns,1:ns);Bk2=temp2(1:ns,(ns+1):(ns+nm));Ck2=temp2((ns+1):(ns+nc), 1:ns);Dk2=temp2((ns+1):(ns+nc), (ns+1):(ns+nm));
Dk=inv(eye(nc)-Dk2*D22)*Dk2;
Bk=Bk2*(eye(nm)-D22*Dk);
Ck=(eye(nc)-Dk*D22)*Ck2;
Ak=Ak2-Bk*inv(eye(nm)-D22*Dk)*D22*Ck;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Acl2=[A+B2*Dk2*C2 B2*Ck2;Bk2*C2 Ak2];
% Bcl2=[B1+B2*Dk2*D21;Bk2*D21];
% Ccl2=[C1+D12*Dk2*C2 D12*Ck2];
% A=Acl2;B=Bcl2;C=Ccl2;%D=Dcl2;
% Close the loop with Lower LFT
plant=ss(A,[B1 B2],[C1;C2],[D11 D12; D21 D22]);
controller=ss(Ak,Bk,Ck,Dk);
sys_cl=lft(plant,controller);
[Acl,Bcl,Ccl,Dcl] = ssdata(sys_cl);
sys_cl=ss(Acl,Bcl,Ccl,0)
norm(sys_cl,2)
% 
% [Acl,Bcl,Ccl,Dcl] = ssdata(sys_cl);
% sys_cl2=ss(Acl,Bcl,Ccl,0);
% norm(sys_cl2,2)
% A=Acl;B=Bcl;C=Ccl;D=Dcl;

% example from Scherer+Gahinet
A=[0 0 1 0;0 0 0 1;-.101 -.1681 -.04564 -.01075;.06082 -2.1407 -.05578 -.1273];
B=[0 0 0;0 0 0;.1179 .1441 .1478;.1441 1.7057 -.7557];
C=[1 0 0 0;0 1 0 0];
sys_ex=ss(A,B,C,0);
norm(sys_ex)
% measure numbers of inputs and outputs 

eta=.0001;    % degree of strict positivity   
ns=size(A,1);   % number of states
nd=size(B,2);  % number of external inputs
no=size(C,1);  % number of regulated outputs
% Calculate H2 Norm

XX=sdpvar(ns);
gamma2=sdpvar(1);
F2=[XX>=eta*eye(ns)]
F2=[F2;A*XX+XX*A'+B*B'<=0]
F2=[F2;trace(C*XX*C')<=gamma2]
optimize(F2,gamma2,OPTIONS)
sqrt(value(gamma2))

% Calculate H2 Norm
XX=sdpvar(ns);
Q=sdpvar(no);
gamma2=sdpvar(1);
F2=[XX>=eta*eye(ns)]
F2=[F2;[A'*XX+XX*A  XX*B; 
        B'*XX         -eye(nd)]<=0]
F2=[F2;[XX C';C Q]>=0]
F2=[F2;trace(Q)<=gamma2]
optimize(F2,gamma2)
sqrt(value(gamma2))



% compare with Matlab built-in functions
sys=ss(A,[B1 B2],[C1;C2],[D11 D12; D21 D22]);
%[K,CL,GAM,INFO]=h2syn(sys,nm,nc,'METHOD','lmi')
[K,CL,GAM,INFO]=h2syn(sys,nm,nc)

