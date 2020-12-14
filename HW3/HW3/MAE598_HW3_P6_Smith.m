%% Homework 3 Problem 6
% Regulator Problem and Optimal State-Feedback
clc
clear all
ops = sdpsettings('solver','mosek','verbose',0);
empty=[];

eta=1e-5;                       % Constant to allow for stricly > calculation during P>0

%% System Definition
A= [
-1  1   0   1   0   1;
-1  -2  -1  0   0   1;
1   0   -2  -1  1   1;
-1  1   -1  -2  0   0;
-1  -1  1   1   -2  -1;
0   -1  0   0   -1  -3
];

B =[
0   -1  -1;
0   0   0;
-1  1   1;
-1  0   0;
0   0   1;
-1  1   1
];

zB = zeros(6,3); % Zero matrix in shape of B

C = [
0   1   0   -1  -1  -1;
0   0   0   -1  0   0;
1   0   0   0   -1  0
];
zC = zeros(3,6); % Zero matrix in shape of B

D = [
0   0   0;
0   0   0;
0   0   0
];
zD = zeros(3);


%% Part (a)
disp("#####Problem 6 Part A Output:######")
disp("Use plant and the regulator problem framework to cosntruct 9-matrix representaiton")
% Regulator Framework 9 Matrix Representation
% (work in document)
P = [
A B zB B;
C D zD D;
zC zD zD eye(3);
C D eye(3) D
];
disp("P=")
disp(P)


%% Part (b)
disp("#####Problem 6 Part B Output:######")
disp("(i) LMI to formulate H_inf optimal state-feedback problem")

% LMI for state feedback (Lecture 9 Theorem 5)

B1 = [B zB];        %6x6
B2 = B;             %6x3
C1 = [C; zC];       %6x6
D11 = [D zD; zD zD];%6x6
D12 = [D; eye(3)];       %6x3

gamma = sdpvar(1,1);
Y=sdpvar(6);
Z = sdpvar(3,6);

optimize_var = gamma;
LMI1 = [Y];
LMI2 = [
    Y*A'+A*Y+Z'*B2'+B2*Z B1 Y*C1'+Z'*D12';
    B1' -gamma*eye(6) D11';
    C1*Y+D12*Z D11 -gamma*eye(6);
    ];

Fun = [
    LMI1>=eta*eye(6);
    LMI2<=-eta*eye(18);
    ];

optimize(Fun,optimize_var,ops);          % Run the optimization

gammaf = value(gamma);          % Return feasable solution gamma
Yf = value(Y);                  % Return feasable solution Y
Zf = value(Z);                  % Return feasable solution Z
F = Zf*inv(Yf);                   % Construct F


disp("F=")
disp(F)
sympref('FloatingPointOutput',true);
Latex_out = latex(sym(F))
disp("(ii) Determine the closed-loop H2 gain")
disp("H_inf gain = ")
disp(gammaf)


%% Part (c)#################################################################################
disp("#####Problem 6 Part C Output:######")
%Define the closed loop system for the Hinf-optimal controller
Acl = A+B2*F;
Bcl = B1;
Ccl = C1+D12*F;
Dcl = D11;

disp("(i) Does it exist")
disp("The H2 gain DNE if the Dcl matrix is nonzero.")
disp("Since Dcl = 0, the H2 gain exists")
disp("Dcl=")
disp(Dcl)

disp("(ii) Use an LMI to determine H2 gain of closed loop system (if it exists)")

X=sdpvar(size(A,1));
gamma = sdpvar(1,1);

optimize_var = gamma;
LMI1=[X];
LMI2 = [Acl*X+X*Acl'+Bcl*Bcl'];
LMI3=[trace(Ccl*X*Ccl')];

Fun = [
    LMI1>=eta*eye(size(A,1));
    LMI2<=eta*eye(6);
    LMI3<=gamma];              % Combine into single constraint function

optimize(Fun,optimize_var,ops);          % Run the optimization

%gammaf = value(gamma)          % Return feasable solution gamma
Xf = value(X);                  % Return feasable solution X
gammaf = sqrt(value(gamma));         % Return feasable solution gamma

disp("X=")
disp(Xf)
sympref('FloatingPointOutput',true);
Latex_out = latex(sym(Xf))
disp("H2 gain = ")
disp(gammaf)
%sys = ss(Acl,Bcl,Ccl,Dcl);
%g = norm(sys,2)
%% Part (d)*
disp("#####Problem 6 Part D Output:######")
disp("% (i) Use an LMI to formulate and solve the H2-optimal state-feedback")
gamma = sdpvar(1,1);
X=sdpvar(size(A,1));
W = sdpvar(6);
Z = sdpvar(3,6);

optimize_var = gamma;
LMI1 = [X];
LMI2 = [[A B2]*[X;Z] + [X Z']*[A';B2']+ B1*B1'];
LMI3 = [X (C1*X+D12*Z)';C1*X+D12*Z W];
LMI4 = [trace(W)];
Fun = [
    LMI1>=eta*eye(6);
    LMI2<=-eta*eye(6);
    LMI3>=eta*eye(12);
    LMI4<=gamma];              % Combine into single constraint function

optimize(Fun,optimize_var,ops);          % Run the optimization

gammaf = sqrt(value(gamma));          % Return feasable solution gamma
Xf = value(X);                  % Return feasable solution X
Wf = value(W);                  % Return feasable solution W
Zf = value(Z);                  % Return feasable solution Z
K= Zf*inv(Xf);                   % Construct K

disp("K=")
disp(K)
sympref('FloatingPointOutput',true);
Latex_out = latex(sym(K))
disp("(ii) Determine the closed-loop H2 gain")
disp("H2 gain = ")
disp(gammaf)

%% Part (e)*
disp("#####Problem 6 Part E Output:######")
disp("(i) Use an LMI to determine the Hinf gain of the closed loop system")

% Define closed loop system
Acl = A+B2*K;
Bcl = B1;
Ccl = C1+D12*K;
Dcl = D11;

%Dialated KYP Lemma

% Define variables
gamma = sdpvar(1,1);
X=sdpvar(6);
optimize_var = gamma;

% Construct LMI constraints
LMI1 = [X];
LMI2 = [
    (Acl'*X+X*Acl)  (X*Bcl) (Ccl');
    (Bcl'*X)  (-gamma*eye(6))  (Dcl');
    (Ccl)  (Dcl)  (-gamma*eye(6))
    ];
Fun = [
    LMI1>=eta*eye(6);
    LMI2<=eta*eye(size(LMI2))];              % Combine into single constraint function
optimize(Fun,optimize_var,ops);          % Run the optimization

% Return feasable solutions
gammaf = sqrt(value(gamma));          % Return feasable solution gamma
Xf = value(X);                  % Return feasable solution X

% Display results
disp("X=")
disp(Xf)
sympref('FloatingPointOutput',true);
Latex_out = latex(sym(Xf))
disp("Hinf gain of closed loop system")
disp(gammaf)

