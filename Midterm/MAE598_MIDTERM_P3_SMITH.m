%% 3A
disp("####################################################################")
disp("3A")
disp("####################################################################")
% % Model 
clear all
clc

A=[
    -1  1   0   1   0   1;
    -1  -2  -1  0   0   1; 
    1   0   -2  -1  1   1;
    -1  -1  1   -3  0   0;
    -1  1   -1  1   -2  -1;
    0   -1  0   0   -1  -2;
    ];
B=[
    0   -1  -1;
    0   0   0;
    -1  1   1;
    -1  0   0;
    0   0   1;
    -1  1   1];
C=[ 
    0   1   0   -1  -1  -1;
    0   0   0   -1  0   0;
    1   0   0   0   -1  0];
D=zeros(3);

disp("9-MATRIX REP###")
I = eye(3);
zB = zeros(size(B));
zC = zeros(size(C));
zD = zeros(size(D));
P = [
    A   zB  B   zB  B;
    C   I   -D  zD   -D;
    zC  zD  zD   zD  I;
    zC  I   zD  zD  zD;
    C   zD  D   I   D
    ];
sympref('FloatingPointOutput',true);
Latex_out = latex(sym(P))
% B1 = [zB B];
% B2 = [zB B];
% C1 = [C; zC];
% C2 = [zC; C];
% 
% D11=[I -D; zD zD];
% D12=[I zD; zD D];
% D21=[zD -D; zD I];
% D22 = [zD zD; I D];


B1 = [zB B zB];
B2 = [B];
C1 = [C; zC];
C2 = [zC; C];

D11=[I   -D  zD;
    zD  zD   zD];
D12=[-D; I];
D21=[I   zD  zD;
    zD  D   I];
D22 = [zD; D];

% Use an LMI or matlab function?
% #########################################################
disp("Hinf and H2 norms of open loop sys")
sys=ss(A,[B1 B2],[C1;C2],[D11 D12; D21 D22]);
disp("The H2 norm is infinite because the system has nonzero feedthrough. ")
disp([D11 D12; D21 D22])
sympref('FloatingPointOutput',true);
Latex_out = latex(sym([D11 D12; D21 D22]))
H2_OpenLoop_Norm = norm(sys,2)
Hinf_OpenLoop_Norm = norm(sys,inf)

%% 3B
disp("####################################################################")
disp("3B")
disp("####################################################################")
% measure numbers of inputs and outputs 

eps=.000001;    % degree of strict positivity   
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
An=sdpvar(ns,ns,'full');
Cn=sdpvar(nc,ns,'full');
Dn=sdpvar(nc,nm,'full');
Bn=sdpvar(ns,nm,'full');

% declare constraints
F=[X1>=eps*eye(ns)]
F=[F;Y1>=eps*eye(ns)]
F=[F;[X1 eye(ns); eye(ns) Y1]>=0];
MAT=[A*Y1+Y1*A'+B2*Cn+Cn'*B2'  (A'+An+(B2*Dn*C2)')'        B1+B2*Dn*D21           (C1*Y1+D12*Cn)'; 
     A'+An+(B2*Dn*C2)'         X1*A+A'*X1+Bn*C2+C2'*Bn'    X1*B1+Bn*D21           (C1+D12*Dn*C2)'  ;
     (B1+B2*Dn*D21)'           (X1*B1+Bn*D21)'             -gamma*eye(nd)          (D11+D12*Dn*D21)'  ;
     C1*Y1+D12*Cn              C1+D12*Dn*C2                D11+D12*Dn*D21         -gamma*eye(no)];

F=[F;MAT<=eps*eye(size(MAT))];
OPTIONS = sdpsettings('solver','mosek','verbose',0);


% Solve the LMI, minimizing gamma
optimize(F,gamma,OPTIONS);
gamman=value(gamma);
Predicted_Hinf_Gain = gamman

disp("####################################################################")
disp("3C Reconstruct Controller with Y2-2*I")
disp("####################################################################")
% retrieve decision variables
X1n=value(X1); 
Y1n=value(Y1);
Ann=value(An);
Bnn=value(Bn);
Cnn=value(Cn);
Dnn=value(Dn);
temp1=[Ann Bnn; Cnn Dnn]-[X1n*A*Y1n zeros(ns,nm); zeros(nc,ns) zeros(nc,nm)];

% Choose X2, Y2, so that X2*Y2=I-X1*Y1;
Y2n=2.0*eye(ns); %############################################ USE Y2=2*I
X2n=(eye(ns)-X1n*Y1n)*2.0;

% Reverse variable substitution
temp2=inv([X2n X1n*B2;zeros(nc,ns) eye(nc)])*temp1*inv([Y2n' zeros(ns,nm); C2*Y1n eye(nm)]);
Ak2=temp2(1:ns,1:ns);Bk2=temp2(1:ns,(ns+1):(ns+nm));Ck2=temp2((ns+1):(ns+nc), 1:ns);Dk2=temp2((ns+1):(ns+nc), (ns+1):(ns+nm));
Dk=inv(eye(nc)-Dk2*D22)*Dk2;
Bk=Bk2*(eye(nm)-D22*Dk);
Ck=(eye(nc)-Dk*D22)*Ck2;
Ak=Ak2-Bk*inv(eye(nm)-D22*Dk)*D22*Ck;
disp("Controller")
sympref('FloatingPointOutput',true);
Latex_out = latex(sym(Ak))
Latex_out = latex(sym(Bk))
Latex_out = latex(sym(Ck))
Latex_out = latex(sym(Dk))
Ak
Bk
Ck
Dk

%% 3D
disp("####################################################################")
disp("3D Ckosed Loop System Matricies")
disp("####################################################################")

C2t=eye(ns); D21t=zeros(ns,nd); D22t=zeros(ns,nc);
controller=ss(Ak,Bk,Ck,Dk);
plant=ss(A,[B1 B2],[C1;C2t],[D11 D12; D21t D22t]);
sys_cl=lft(plant,controller);

%[Acl,Bcl,Ccl,Dcl] = ssdata(sys_cl);
%sys_cl=ss(Acl,Bcl,Ccl,0);
Hinf_norm = norm(sys_cl,inf)
Acl = sys_cl.A
Bcl = sys_cl.B
Ccl = sys_cl.C
Dcl = sys_cl.D
% Acl = A
% Bcl = [B1 B2]
% Ccl = [C1;C2]
% Dcl = [D11 D12; D21 D22]
sympref('FloatingPointOutput',true);
Latex_out = latex(sym(Acl))
Latex_out = latex(sym(Bcl))
Latex_out = latex(sym(Ccl))
Latex_out = latex(sym(Dcl))
% 
% % Close the loop with Lower LFT
% plant=ss(A,[B1 B2],[C1;C2],[D11 D12; D21 D22]);
% controller=ss(Ak,Bk,Ck,Dk);
% sys_cl=lft(plant,controller);
% Controller_Hinf_Gain = norm(sys_cl,Inf)
% 
% % compare with Matlab built-in functions
% sys=ss(A,[B1 B2],[C1;C2],[D11 D12; D21 D22]);
% [K,CL,GAM,INFO]=hinfsyn(sys,nm,nc,'METHOD','lmi');
% [K,CL,GAM,INFO]=hinfsyn(sys,nm,nc);


%% 3E
disp("####################################################################")
disp("3E Hinf Optimal Full State Feeback Controller")
disp("####################################################################")



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

%OPTIONS = sdpsettings('solver','sedumi');

% Solve the LMI, minimizing gamma
optimize(F,gamma,OPTIONS);
gamman=value(gamma);
Xn = value(X)
Zn = value(Z)
K=Zn*inv(Xn)
Latex_out = latex(sym(K))
Predicted_Hinf_Gain = gamman

%% 3F
disp("####################################################################")
disp("3F Hinf Optimal Observer")
disp("####################################################################")

gamma=sdpvar(1);               % represents the bound on the H-infinity norm of the CL system.
P=sdpvar(ns);
G=sdpvar(ns,nr);

MAT = [P*A+A'*P-G*C2-C2'*G' P*B1-G*D21 C1';
        (P*B1-G*D21)' -gamma*eye(9) D11';
        C1 D11 -gamma*eye(6)]
F=[MAT<=0];
F=[F;P>=eta*eye(ns)];

optimize(F,gamma,OPTIONS);

Pn = value(P)
Gn = value(G);
gamman=value(gamma);
Predicted_Hinf_Gain = gamman
L=inv(Pn)*Gn
Latex_out = latex(sym(L))
%% 3E
disp("####################################################################")
disp("Combine the optimal observer and optimal state-feedback controller and construct the closedloop system. What is the resulting H? gain? You may use the Matlab norm command to determine the norm of the closed-loop system.")
disp("####################################################################")
%K = [Ak Bk; Ck Dk]
% At = [ Acl+Bcl*K -Bcl*K
%       zeros(size(Acl)) Acl+L*Ccl ];
% Bt = [ B;
%       zeros(size(Bcl))];
% Ct = [ C    zeros(size(Ccl)) ];
% A=A;
% B=B*K;
% C=-L*Ccl;
% D=Acl+L*Ccl+Bcl*K;
At = [ A+B2*K -B2*K
      zeros(size(A)) A+L*C2 ];
Bt = [ B2;
      zeros(size(B2))];
Ct = [ C2    zeros(size(C2)) ];
sys = ss(At,Bt,Ct,0);
%sys_cl=ss(A,B,C,D);
ObserverStateFeedback_Hinf_Gain = norm(sys_cl,inf)













