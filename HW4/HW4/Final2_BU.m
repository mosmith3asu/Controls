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
g=9.8;
lh=1*g; %meters
ls=0.5;%meters
m0=60; %kg
k0=20; %Nm/rad

nm=0.2;
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

% Define 9-matrix rep%%%%%%%%%%%%%%%%%%%%%
% A=[0 1; (k0*ls+m0*lh) 0];
% B1 = [0 0 0 0; 
%     -nm nm/(lh^2*m0) nk/(lh^2*m0) nk/(lh^2*m0)
%     ];
% B2 = [0; -k0*ls/(lh^2*m0)];
% 
% C1=[(k0+m0*lh)/(lh^2*m0) 0; lh*m0 0; 0 0; ls*k0 0];
% C2 = [1 0];
% 
% D11 = [-nm nm/(lh^2*m0) nk/(lh^2*m0) nk/(lh^2*m0); 0 0 0 0; 0 0 0 0; 0 0 0 0];
% D12 = [-k0*ls/(lh^2*m0); 0; -ls*k0; 0];
% D21 = [0 0 0 0];
% D22 = [0];
% 
% B = [B1 B2];
% C = [C1; C2];
% D = [D11 D12; D21 D22];
% 
% %A=P11;
% eta=.001
% n=2; % number of states
% np=3; % number of uncertain outputs
% nq=3; % number of uncertain inputs
% 
% % Declare Scalings
% th1=sdpvar(1);
% th2=sdpvar(1);
% 
% %Th=[th1 0 0;0 th2 0; 0 0 th3]
% Th=diag([th1;th1;th2;th2],0)
% %Th=blkdiag(blkdiag(th1,th2),th3)
% F=[];
% F=[F;Th>=0];
% 
% % Lyap variable
% tol=.1
% gam_u=100;
% gam_l=0
% X=sdpvar(n);
% F=[F;X>=eta*eye(n)]
% B=B1;C=C
% gam_new=gam_u;
% err=gam_u
% while err>tol
%     MAT=[A'*X+X*A X*B1;B1'*X -Th]+1/gam_new/gam_new*[C1 D11]'*Th*[C1 D11];
%     Ftemp=[F;MAT<=0];
%     DIAG=optimize(Ftemp,[],options);
%     if DIAG.problem==0
%         gam_u=gam_new;
%     else
%         gam_l=gam_new;
%     end
%     gam_new=(gam_u+gam_l)/2;
%     err=gam_u-gam_l;
% end

%% Stabilizing Controller Synthesis
% Define 9-matrix rep%%%%%%%%%%%%%%%%%%%%%
A=[0 1; (k0*ls+m0*lh) 0];
B1 = [0 0 0 0; 
    -nm nm/(lh^2*m0) nk/(lh^2*m0) nk/(lh^2*m0)
    ];
B2 = [0; -k0*ls/(lh^2*m0)];
B3= [0; -k0*ls/(lh^2*m0)];


C1=[(k0+m0*lh)/(lh^2*m0) 0; lh*m0 0; 0 0; ls*k0 0];
C2 = [1 0;
    0 0];
C3 = [1 0];

D11 = [-nm nm/(lh^2*m0) nk/(lh^2*m0) nk/(lh^2*m0); 0 0 0 0; 0 0 0 0; 0 0 0 0];
D12 = [-k0*ls/(lh^2*m0); 0; -ls*k0; 0];
D13=[0;0;0;0];

D21 = [0 0 0 0;0 0 0 0];
D22 = [0;0];
D23= [0;1];

D31=[0 0 0 0];
D32=0;
D33=0;
% Controlled Formulation

eta=.001;
n=2; % number of states
np=4; % number of uncertain outputs
nq=4; % number of uncertain inputs
nc=size(D13,2);

%==================================================
% Declare Scalings
th1=sdpvar(1);
th2=sdpvar(1);
Th=diag([th1;th1;th2;th2],0);
F=[];
F=[F;Th>=0];

X=sdpvar(n);
Z=sdpvar(nc,n,'full');
F=[F;X>=eta*eye(n)];
MAT=[A*X+B3*Z+X*A'+Z'*B3' X*C1'+Z'*D13';C1*X+D13*Z -Th]+[B1; D11]*Th*[B1' D11'];
F=[F;MAT<=0];
DIAG=optimize(F,[],options);
Zn=value(Z);
Xn=value(X);
K=Zn*inv(Xn)

%==================================================
% Part c - Close the loop with stabilizing controller and recalculate mu

Acl=A+B3*K; C1cl=C1+D12*K;
% Declare Scalings
F=[];
F=[F;Th>=0];

% Lyap variable
tol=.01
gam_u=100;
gam_l=0
gam_new=gam_u;
err=gam_u;
while err>tol
    clear X F MAT th1 th2 Th
th1=sdpvar(1);
th2=sdpvar(1);
    Th=diag([th1;th1;th2;th2],0);
    X=sdpvar(n);
    F=[];
    F=[F;X>=eta*eye(n)];
    MAT=[Acl'*X+X*Acl X*B1;B1'*X -Th]+1/gam_new/gam_new*[C1cl D11]'*Th*[C1cl D11];
    Ftemp=[F;MAT<=0];
    DIAG=optimize(Ftemp,[],options);
    if DIAG.problem==0
        gam_u=gam_new;
    else
        gam_l=gam_new;
    end
    gam_new=(gam_u+gam_l)/2;
    err=gam_u-gam_l;
end

% MAT=[A*X+B3*Z+X*A'+Z'*B3'+B1*Th*B1' (C2*X+D23*Z)'    X*C1'+Z'*D13' ;
%     C2*X+D23*Z            -gamma*eye(nro)  zeros(nro,np); 
%     C1*X+D13*Z             zeros(np,nro) -Th];
%     F=[F;MAT<=0];
%     DIAG=optimize(F,gamma);

%==================================================
% Part d- Now do Hinfinity optimal robust control
disp('Hinfinity optimal robust control')
nro=size(C2,1);


%==================================================
% Define Scalings Variables
gamma=sdpvar(1);
th1=sdpvar(1);
th2=sdpvar(1);
Th=diag([th1;th1;th2;th2],0);

%==================================================
% Define Constraints
F=[];
F=[F;Th>=0];
P=sdpvar(n);
Z=sdpvar(nc,n,'full');

F=[F;P>=eta*eye(n)];
MAT=[
    (A*P+B*Z+P*A'+Z'*B'+B2*B2'+M*Th*M') (C*X+D22*Z)'    (P*N'+Z'*D12');
    (C*P+D22*Z)            (-gamma*eye(nro))  (zeros(nro,np));
    (N*P+D12*Z)             (zeros(np,nro)) (-Th)];
F=[F;MAT<=0];

%==================================================
% Run optimization
objective=gamma;
sol=optimize(F,objective,options);

%==================================================
% Return feasable solutions
Zn=value(Z);
Xn=value(P);
K=Zn*inv(Xn);
Thn=value(Th);
gamman=sqrt(value(gamma));



A
B=B3
B2 %maybe
M=B1
C=C2
N=C1
Q=0;
D12 =D13;
D22 = D23

Latex_out = latex(sym(A))
Latex_out = latex(sym(B))
Latex_out = latex(sym(B2))
Latex_out = latex(sym(C))
Latex_out = latex(sym(D12))
Latex_out = latex(sym(D22))
Latex_out = latex(sym(M))
Latex_out = latex(sym(N))
Latex_out = latex(sym(Q))


while err>tol
    clear X F MAT th1 th2 Th
    % =========================================
    % Declare Scalings
    gamma=sdpvar(1);
    th1=sdpvar(1);
    th2=sdpvar(1);
    Th=diag([th1;th1;th2;th2],0);
    
    % =========================================
    % Declare sdp Vars
    X=sdpvar(n);
    
    % =========================================
    % Constraints
    F=[];
    F =[F;X>=eta*eye(n)];
    MAT=[Acl'*X+X*Acl X*Bcl;Bcl'*X -Th]+1/gam_new/gam_new*[C1cl Dcl]'*Th*[C1cl Dcl];
    Ftemp=[F;MAT<=0];
    sol=optimize(Ftemp);
    if sol.problem==0
        gam_u=gam_new;
    else
        gam_l=gam_new;
    end
    gam_new=(gam_u+gam_l)/2;
    err=gam_u-gam_l;
end
gam_u
