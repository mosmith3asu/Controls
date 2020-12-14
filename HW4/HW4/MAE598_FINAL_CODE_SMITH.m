%enforcing quadratic stability using the polytopic framework%% Problem 3
% Definition of variables
clear all
clc
options = sdpsettings('solver','mosek','verbose',0);


%% P3.a
%Define given variables
g=9.8;
lh=1*g; %meters
ls=0.5;%meters
m0=60; %kg
k0=20; %Nm/rad

nm=0.2;
nk=0.1;
alpha=0.5;
% Define 9-matrix rep%%%%%%%%%%%%%%%%%%%%%
A=[0 1; (k0*ls+m0*lh-alpha*ls*k0) 0];
B1 = [0 0 0 0; 
    -nm nm/(lh^2*m0) nk/(lh^2*m0) nk/(lh^2*m0)
    ];
B2 = [0; -alpha*k0*ls/(lh^2*m0)];
B3= [0; -alpha*k0*ls/(lh^2*m0)];


C1=[(k0*ls+m0*lh-alpha*ls*k0)/(lh^2*m0) 0; lh*m0 0; 0 0; ls*k0 0];
C2 = [1 0;
    0 0];
C3 = [1 0];

D11 = [-nm nm/(lh^2*m0) nk/(lh^2*m0) nk/(lh^2*m0); 0 0 0 0; 0 0 0 0; 0 0 0 0];
D12 = [-alpha*k0*ls/(lh^2*m0); 0; -ls*k0; 0];
D13=[0;0;0;0];

D21 = [0 0 0 0;0 0 0 0];
D22 = [0;0];
D23= [0;1];

D31=[0 0 0 0];
D32=0;
D33=0;

%%
disp('###################################################')
disp('Sturctured singular value')
eta=.001;
n=2; % number of states
np=4; % number of uncertain outputs
nq=4; % number of uncertain inputs

% Declare Scalings


% =========================================
% Bisection gamma
tol=.01;
gam_feasable=100;
gam_l=0;
gam_new=gam_feasable;
err=gam_feasable;

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
    F=[F;X>=eta*eye(n)];
    MAT=[A'*X+X*A X*B1;B1'*X -Th]+1/gam_new/gam_new*[C1 D11]'*Th*[C1 D11];
    Ftemp=[F;MAT<=0];
    DIAG=optimize(Ftemp,[],options);
    if DIAG.problem==0
        gam_feasable=gam_new;
    else
        gam_l=gam_new;
    end
    gam_new=(gam_feasable+gam_l)/2;
    err=gam_feasable-gam_l;
    
    
end

disp('SSV open loop')
gam_feasable
mu = sqrt(gam_feasable)
disp('###################################################')
%% Controller Sythesis
% Define 9-matrix rep%%%%%%%%%%%%%%%%%%%%%
disp('###################################################')
disp('Controller Synthesis')
% Controlled Formulation

eta=.001;
n=2; % number of states
np=4; % number of uncertain outputs
nq=4; % number of uncertain inputs
nc=size(D13,2);

%==================================================
% Part d- Now do Hinfinity optimal robust control
disp('Hinfinity optimal robust control')
nro=size(C2,1);

% Declare Scalings
gamma=sdpvar(1);
th1=sdpvar(1);
th2=sdpvar(1);
Th=diag([th1;th1;th2;th2],0);
F=[];
F=[F;Th>=0];
X=sdpvar(n);
Z=sdpvar(nc,n,'full');
F=[F;X>=eta*eye(n)];
MAT=[A*X+B3*Z+X*A'+Z'*B3'+B2*B2'+B1*Th*B1' (C2*X+D23*Z)'    X*C1'+Z'*D13' ;
    C2*X+D23*Z            -gamma*eye(nro)  zeros(nro,np);
    C1*X+D13*Z             zeros(np,nro) -Th];
F=[F;MAT<=0];
DIAG=optimize(F,gamma,options);
sympref('FloatingPointOutput',true);

Zn=value(Z);
% Latex_out = latex(sym(Zn))
Xn=value(X);
% Latex_out = latex(sym(Xn))
K=Zn*inv(Xn)
% Latex_out = latex(sym(K))
Thn=value(Th)
% Latex_out = latex(sym(Thn))
% gamman=sqrt(value(gamma))

%% SSV

% =========================================
% Close the loop with stabilizing controller K
Acl=A+B3*K; C1cl=C1+D12*K;

% =========================================
% Declare Constraint on Theta
F=[];
F=[F;Th>=0];

% =========================================
% Bisection gamma
tol=.01;
gam_feasable=100;
gam_l=0;
gam_new=gam_feasable;
err=gam_feasable;

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
    F=[F;X>=eta*eye(n)];
    MAT=[Acl'*X+X*Acl X*B1;B1'*X -Th]+1/gam_new/gam_new*[C1cl D11]'*Th*[C1cl D11];
    Ftemp=[F;MAT<=0];
    DIAG=optimize(Ftemp,[],options);
    if DIAG.problem==0
        gam_feasable=gam_new;
    else
        gam_l=gam_new;
    end
    gam_new=(gam_feasable+gam_l)/2;
    err=gam_feasable-gam_l;
end
disp('###################################################')
disp('SSV closed loop controller')
gam_feasable
mu = sqrt(gam_feasable)
disp('###################################################')