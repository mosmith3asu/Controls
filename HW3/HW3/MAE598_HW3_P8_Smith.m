%% Homework 3 Problem 6
% Mixed Norm Optimization
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
P = [
A B zB B;
C D zD D;
zC zD zD eye(3);
C D eye(3) D
];

% Define the 9-natrux representation components
B1 = [B zB];        %6x6
B2 = B;             %6x3
C1 = [C; zC];       %6x6
C2 = C;
D11 = [D zD; zD zD];%6x6
D12 = [D; eye(3)];       %6x3
D21 = [D eye(3)];
D22 = D;
n=6;

%% Part (a)
disp("#####Problem 8 Part A Output:######")
disp("(i) Reformulate the Hinf Output feedback problem so that it minimizes ||S(P,K)||^2_Hinf")

%% Part (b)
disp("#####Problem 8 Part B Output:######")
disp("(i) Use an LMI to formulate and solve the optimal output-feedback problem minimizing both the H2 and Hinf gains, ")
disp("    giving equal weight to each. min_K ||S(P;K)||^2_H2 +||S(P;K)||^2_Hinf")

X1=sdpvar(6);
Y1=sdpvar(6);
Z=sdpvar(6);
An=sdpvar(6,6);
Bn=sdpvar(6,3,'full');
Cn=sdpvar(3,6,'full');
Dn=sdpvar(3,3);
beta1=sdpvar(1);
beta2=sdpvar(1);

%H2 constraints
LMI1 = [(D11+D12*Dn*D21)];
LMI2 = [trace(Z)];
LMI3=[(A*Y1+Y1*A'+B2*Cn+Cn'*B2')  ((A'+An+(B2*Dn*C2)')') ((B1+B2*Dn*D21));
    (A'+An+(B2*Dn*C2)') (X1*A+A'*X1+Bn*C2+C2'*Bn')  (X1*B1+Bn*D21);
    ((B1+B2*Dn*D21)')   ((X1*B1+Bn*D21)') (-eye(6))];
    
LMI4=[(Y1)  (eye(6))  ((C1*Y1+D12*Cn)');
    (eye(6))  (X1)  ((C1+D12*Dn*C2)');
    (C1*Y1+D12*Cn)	(C1+D12*Dn*C2)	(Z)];

% Hinf constraint
LMI5=[(A*Y1+Y1*A'+B2*Cn+Cn'*B2')  ((A'+An+(B2*Dn*C2)')')	(B1+B2*Dn*D21)	((C1*Y1+D12*Cn)');
    (A'+An+(B2*Dn*C2)')	(X1*A+A'*X1+Bn*C2+C2'*Bn')	(X1*B1+Bn*D21)	((C1+D12*Dn*C2)');
    ((B1+B2*Dn*D21)')	(X1*B1+Bn*D21)'	(-beta2*eye(6))	((D11+D12*Dn*D21)');
    (C1*Y1+D12*Cn)	(C1+D12*Dn*C2)	(D11+D12*Dn*D21)	(-beta2*eye(6))];
     
Fun = [
    LMI1==0;
    LMI2<= beta1;
    LMI3 <=0;
    LMI4 >=0;
    LMI5 <=0];              % Combine into single constraint function

beta1_weight = 1.0;   % weight for min_K ||S(P;K)||^2_H2 
beta2_weight = 1.0;   % weight for min_K||S(P;K)||^2_Hinf
optimize_var= beta1_weight*beta1 + beta2_weight*beta2;

optimize(Fun,optimize_var,ops);          % Run the optimization


% Return feasable results
X1f = value(X1);                  % Return feasable solution X1
Y1f = value(Y1);                  % Return feasable solution Y1
Anf = value(An);                  % Return feasable solution An
Bnf = value(Bn);                  % Return feasable solution Bn
Cnf = value(Cn);                  % Return feasable solution Cn
Dnf = value(Dn);                  % Return feasable solution Dn
betaf = value(optimize_var);                  % Return feasable solution beta
beta1f = value(beta1);                  % Return feasable solution beta1
beta2f = value(beta2);                  % Return feasable solution beta2

gamma1=sqrt(beta1f);              
gamma2=beta2f;             

H2_LMI = gamma1;
Hinf_LMI = gamma2;

disp("Minimized H2 Gain = ")
disp(H2_LMI)
disp("Minimized Hinf Gain = ")
disp(Hinf_LMI)

%% Part (c)
disp("#####Problem 8 Part C Output:######")
disp("(i) Construct the corresponding controller")
Y2 = eye(n);
X2 = eye(n)-X1f*Y1f;

K2 = (inv([X2 X1f*B2; zeros(3,n) eye(3)]))*...
    ([Anf Bnf; Cnf Dnf]- [X1f*A*Y1f zeros(6,3); zeros(3,9)])*...
    (inv([Y2' zeros(6,3); C2*Y1f eye(3)]));


Ak2 = K2(1:6,1:6);
Bk2 = K2(1:6,7:9);
Ck2 = K2(7:9,1:6);
Dk2 = K2(7:9,7:9);

Dk = inv(eye(3)+Dk2*D22)*Dk2;
Bk = Bk2*(eye(3)-D22*Dk);
Ck = (eye(3)-Dk*D22)*Ck2;
Ak = Ak2-Bk*inv(eye(3)-D22*Dk)*D22*Ck;

K = [Ak, Bk;Ck, Dk];
disp("K=")
disp(K)
sympref('FloatingPointOutput',true);
Latex_out = latex(sym(K))

% Define closed loop system
Q = inv(eye(3)- D22*Dk);
Acl = [A zeros(6);zeros(6) Ak ] + [B2, zeros(6,3);zeros(6,3),Bk]...
      *inv([eye(3), -Dk; -D22, eye(3)])*[zeros(3,6), Ck; C2,zeros(3,6)];
Bcl = [B1+ B2*Dk*Q*D21; ...
      Bk*Q*D21];
Ccl = [C1, zeros(6)]+ [D12, zeros(6,3)]*inv([eye(3), -Dk;...
      -D22, eye(3)])*[zeros(3,6), Ck; C2,zeros(3,6)];
Dcl = D11+D12*Dk*Q*D21;

sys = ss(Acl, Bcl, Ccl, 0);

disp("(ii) Determine and compare the resulting H2 gain to the gains predicted by the LMI")
H2_ctrl = norm(sys,2);
disp("(Returned from Matlab System Analysis)")
disp("H_2 gain = ")
disp(H2_ctrl)
disp("Therefore the controller approx achieves the predicted gain since the following terms are approximatly equal:")
disp(strcat("   H2_(LMI) = H2_(controller) => ",string(H2_LMI), "=", string(H2_ctrl)))

disp("(iii) Determine and compare the resulting Hinf gain to the gains predicted by the LMI")
Hinf_ctrl = norm(sys,inf);
disp("(Returned from Matlab System Analysis)")
disp("H_inf gain = ")
disp(Hinf_ctrl)
disp("Therefore the controller approx achieves the predicted gain since the following terms are approximatly equal:")
disp(strcat("   Hinf_(LMI) = Hinf_(controller) => ",string(Hinf_LMI), "=", string(Hinf_ctrl)))

%% Part (d)
%disp("#####Problem 6 Part D Output:######")
%disp("(i) Compare these gains to those produced by pure H2-optimal output feedback.")
%disp("(ii) Compare these gains to those produced by pure Hinf-optimal output feedback. ")









