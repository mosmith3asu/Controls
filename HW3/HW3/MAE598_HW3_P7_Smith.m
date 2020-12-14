%% Homework 3 Problem 6
% Optimal Output Feedback
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
n=6;

% Define the 9-natrux representation components
B1 = [B zB];        %6x6
B2 = B;             %6x3
C1 = [C; zC];       %6x6
C2 = C;
D11 = [D zD; zD zD];%6x6
D12 = [D; eye(3)];       %6x3
D21 = [D eye(3)];
D22 = D;

%% Part (a)*
disp("#####Problem 7 Part A Output:######")
disp("(i) Use an LMI to formulate and solve the Hinf-optimal ouput-feedback problem. What Hinf gain do you  find?")


% Define variables
X1 = sdpvar(n);    % Define Z
Y1 = sdpvar(n);    % Define Y
An = sdpvar(n);    % Define An
Bn = sdpvar(n,3);    % Define Bn
Cn = sdpvar(3,n);    % Define Cn
Dn = sdpvar(3,3);    % Define Dn
gamma = sdpvar(1);

% Define LMIs and Constraints
optimize_var = gamma;
LMI1 = [X1 eye(n); eye(n) Y1];  
LMI2 = [(A*Y1+Y1*A'+B2*Cn+Cn'*B2') (A'+An+(B2*Dn*C2)')' ((B1+B2*Dn*D21)')' (C1*Y1+D12*Cn)';
       (A'+An+(B2*Dn*C2)') (X1*A+A'*X1+Bn*C2+C2'*Bn')  ((X1*B1+Bn*D21)')' (C1+D12*Dn*C2)';
       (B1+B2*Dn*D21)' (X1*B1+Bn*D21)'  (-gamma*eye(6))  (D11+D12*Dn*D21)';
       (C1*Y1+D12*Cn) (C1+D12*Dn*C2)  (D11+D12*Dn*D21)  (-gamma*eye(6))];
Fun = [
    LMI1>=eta*eye(2*n);
    LMI2<=-eta*eye(4*n)];              % Combine into single constraint function

% Optimize
optimize(Fun,optimize_var,ops);          % Run the optimization

% Feasable Results
X1f = value(X1);                  % Return feasable solution X1
Y1f = value(Y1);                  % Return feasable solution Y1
Anf = value(An);                  % Return feasable solution An
Bnf = value(Bn);                  % Return feasable solution Bn
Cnf = value(Cn);                  % Return feasable solution Cn
Dnf = value(Dn);                  % Return feasable solution Dn
gammaf = value(gamma);          % Return feasable solution gamma

Hinf_LMI = gammaf;

%(ii) What Hinf gain do you  find? 
disp("(Returned from LMI)")
disp("Hinf Gain=")
disp(Hinf_LMI) 

%% Part (b)
disp("#####Problem 7 Part B Output:######")
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
disp("(ii) show that it achieves the predicted closed-loop Hinf gain.") 
Q = inv(eye(3)- D22*Dk);
Acl = [A zeros(6);zeros(6) Ak ] + [B2, zeros(6,3);zeros(6,3),Bk]...
      *inv([eye(3), -Dk; -D22, eye(3)])*[zeros(3,6), Ck; C2,zeros(3,6)];
Bcl = [B1+ B2*Dk*Q*D21; ...
      Bk*Q*D21];
Ccl = [C1, zeros(6)]+ [D12, zeros(6,3)]*inv([eye(3), -Dk;...
      -D22, eye(3)])*[zeros(3,6), Ck; C2,zeros(3,6)];
Dcl = D11+D12*Dk*Q*D21;


sys = ss(Acl, Bcl, Ccl, Dcl);
Hinf_ctrl = norm(sys,inf);
disp("(Returned from Matlab System Analysis)")
disp("H_inf gain = ")
disp(Hinf_ctrl)
disp("Therefore the controller approx achieves the predicted gain since the following terms are approximatly equal:")
disp(strcat("   Hinf_(LMI) = Hinf_(controller) => ",string(Hinf_LMI), "=", string(Hinf_ctrl)))

%% Part (c)
disp("#####Problem 7 Part C Output:######")
disp("(i) Use an LMI to formulate and solve the H2-optimal ouput-feedback problem. What H2 gain diid you find?")

% Define variables
X1 = sdpvar(6);    % Define X1
Y1 = sdpvar(6);    % Define Y1
Z = sdpvar(6);    % Define Z
An = sdpvar(6,6,'full');    % Define An
Bn = sdpvar(6,3,'full');    % Define Bn
Cn = sdpvar(3,6,'full');    % Define Cn
Dn = sdpvar(3,3,'full');    % Define Dn
gamma = sdpvar(1,1);

% Define LMIs and Constraints
optimize_var = gamma;
LMI1=[(A*Y1 + Y1*A' + B2*Cn + (Cn')*B2') (A'+An+(B2*Dn*C2)')' ((B1+B2*Dn*D21)')';
      (A'+An+(B2*Dn*C2)')  (X1*A + (A')*X1 + Bn*C2 + (C2')*Bn')  ((X1*B1+Bn*D21)')';
      ((B1+B2*Dn*D21)') ((X1*B1+Bn*D21)')  (-eye(6))];  
  
LMI2 = [(Y1) (eye(6))' (C1*Y1+D12*Cn)';...
      (eye(6)) (X1)  (C1+D12*Dn*C2)';...
      (C1*Y1+D12*Cn) (C1+D12*Dn*C2)  (Z)]; 
LMI3 = trace(Z);
LMI4 = D11+D12*Dn*D21;

Fun =  [
    LMI1<=eta*eye(18);
    LMI2>=eta*eye(18);
    LMI3<=gamma;
    LMI4==0];

% Optimize
optimize(Fun,optimize_var,ops);          % Run the optimization

% Feasable Results
X1f = value(X1);                  % Return feasable solution X1
Y1f = value(Y1);                  % Return feasable solution Y1
Anf = value(An);                  % Return feasable solution An
Bnf = value(Bn);                  % Return feasable solution Bn
Cnf = value(Cn);                  % Return feasable solution Cn
Dnf = value(Dn);                  % Return feasable solution Dn
gammaf = sqrt(value(gamma));          % Return feasable solution gamma
H2_LMI = gammaf;

%disp("(ii) What H2 gain do you  find? ")

disp("(Returned from LMI) H_2 gain = ")
disp(H2_LMI)

%% Part (d)
disp("#####Problem 7 Part D Output:######")
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
disp("(ii) show that it achieves the predicted closed-loop H_2 gain.") 
Q = inv(eye(3)- D22*Dk);
Acl = [A zeros(6);zeros(6) Ak ] + [B2, zeros(6,3);zeros(6,3),Bk]...
      *inv([eye(3), -Dk; -D22, eye(3)])*[zeros(3,6), Ck; C2,zeros(3,6)];
Bcl = [B1+ B2*Dk*Q*D21; ...
      Bk*Q*D21];
Ccl = [C1, zeros(6)]+ [D12, zeros(6,3)]*inv([eye(3), -Dk;...
      -D22, eye(3)])*[zeros(3,6), Ck; C2,zeros(3,6)];
Dcl = D11+D12*Dk*Q*D21;


sys = ss(Acl, Bcl, Ccl, 0);
H2_ctrl = norm(sys,2);

disp("(Returned from Matlab System Analysis) H_2 gain = ")
disp(H2_ctrl)
disp("Therefore the controller approx achieves the predicted gain since the following terms are approximatly equal:")
disp(strcat("   H2_(LMI) = H2_(controller) => ",string(H2_LMI), "=", string(H2_ctrl)))











