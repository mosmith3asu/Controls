%%1A
clear all
clc
disp("####################################################################")
disp("1A")
disp("####################################################################")
A=[
    1   1   0   1   0   1;
    -1  -1  -1  0   0   1;
    1   0   1   -1  1   1;
    -1  1   -1  -1  0   0;
    -1  -1  1   1   1   -1;
    0   -1  0   0   -1  -1];
B1=[
    0   -1  -1;
    0   0   0;
    -1  1   1;
    -1  0   0;
    0   0   1;
    -1  1   1];
C1=[
    0   1   0   -1  -1  -1;
    0   0   0   -1  0   0;
    1   0   0   0   -1  0;
    ];
C2= zeros(size(C1));
B2a =[
    0   0   0;
    -1  0   1;
    -1  1   0;
    1   -1  0;
    -1  0   -1;
    0   1   1];
B2b = [
    0   0   0;
    -1  0   1;
    -1  1   0;
    1   1   0;
    1   0   1;
    0   -3  -1];
D11 = [1 2 3; 0 0 0; 0 0 0];
D12 = [0 0 0; 0 0 0; 0 0 0];
D21 = [0 0 0; 0 0 0; 0 0 0];
D22 = [0 0 0; 0 0 0; 0 0 0];


% measure numbers of inputs and outputs 

eta=.0001;    % degree of strict positivity   
ns=size(A,1);   % number of states
nc=size(B2a,2);  % number of actuators
nd=size(B1,2);  % number of external inputs
nr=size(C1,1);  % number of regulated outputs
C2t=eye(ns); D21t=zeros(ns,nd); D22t=zeros(ns,nc);
nm=size(C2,1);  % number of sensors


% LMI for Quadratic Polytopic H?-OptimalState-Feedback Control
% Lecture 13 Theorem 10.


% Declare the variables
gamma=sdpvar(1);               % represents the bound on the H-infinity norm of the CL system.
Y=sdpvar(ns);
Z=sdpvar(nc,ns,'full');

% declare constraints
m11 = [Y*A'+A*Y+Z'*(B2a+B2b)'+(B2a+B2b)*Z];
m21 = [B1'];
m31 = [C1*Y+D12*Z];
m22 = [-gamma*eye(3)];
m32 = [D11];
m33 = [-gamma*eye(3)];

MAT=[m11 m21' m31';
    m21 m22 m32';
    m31 m32 m33];

F=[MAT<=0];
F=[F;Y>=eta*eye(ns)];

OPTIONS = sdpsettings('solver','mosek','verbose',0);

% Solve the LMI, minimizing gamma
optimize(F,gamma,OPTIONS);
gamman=value(gamma)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% retrieve decision variables
Yn=value(Y);
Zn=value(Z);
F=Zn*inv(Yn)

%% 1B
disp("####################################################################")
disp("1B")
disp("####################################################################")

% measure numbers of inputs and outputs 
B2= B2a;
% measure numbers of inputs and outputs 

eta=.0001;    % degree of strict positivity   
ns=size(A,1);   % number of states
nc=size(B2,2);  % number of actuators
nd=size(B1,2);  % number of external inputs
nr=size(C1,1);  % number of regulated outputs
C2t=eye(ns); D21t=zeros(ns,nd); D22t=zeros(ns,nc);
nm=size(C2,1);  % number of sensors


% H-infinity State Feedback Controller Synthesis


% measure numbers of inputs and outputs 

eta=.0001;    % degree of strict positivity   
ns=size(A,1);   % number of states
nc=size(B2a,2);  % number of actuators
nd=size(B1,2);  % number of external inputs
nr=size(C1,1);  % number of regulated outputs
C2t=eye(ns); D21t=zeros(ns,nd); D22t=zeros(ns,nc);
nm=size(C2,1);  % number of sensors


% LMI for Quadratic Polytopic H?-OptimalState-Feedback Control
% Lecture 13 Theorem 10.


% Declare the variables
gamma=sdpvar(1);               % represents the bound on the H-infinity norm of the CL system.
Y=sdpvar(ns);
Z=sdpvar(nc,ns,'full');

% declare constraints
m11 = [Y*A'+A*Y+Z'*(B2a)'+(B2a)*Z];
m21 = [B1'];
m31 = [C1*Y+D12*Z];
m22 = [-gamma*eye(3)];
m32 = [D11];
m33 = [-gamma*eye(3)];

MAT=[m11 m21' m31';
    m21 m22 m32';
    m31 m32 m33];

F=[MAT<=0];
F=[F;Y>=eta*eye(ns)];

OPTIONS = sdpsettings('solver','mosek','verbose',0);

% Solve the LMI, minimizing gamma
optimize(F,gamma,OPTIONS);
gammanA=value(gamma)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% retrieve decision variables
Yn=value(Y);
Zn=value(Z);
F=Zn*inv(Yn)
