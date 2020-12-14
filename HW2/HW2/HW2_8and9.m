clc;
clear all;
%% HW 2 Problem 7 Main Function Script

% Use MOSEK as Solver
ops = sdpsettings('solver','mosek','verbose',0);
empty=[];

% Variable Definitions
A=[ 0 1 0;-30 0 20;2000 0 -2000];
B=[0;0.1;0];
C=[0 0 1];
D = 0;

n=length(A);                    % Get the size of matrix A
eta=1e-5;                       % Constant to allow for stricly > calculation during P>0
stricly_greater0 =eta*eye(n); 	% Variable to estimate strictly zero constraint

%% Problem 8.a
% Use an LMI to design a feedback controller so that the closed-loop achieves
% a settling time of 10s or less, an overshoot of 10% or less and a rise-time 
% of 2s or greater. Print K and the eigenvalues of A+BK. 
disp("HW2 Problem 8.a Feedback Controller Design - Matlab Output")

ts = 10.0;        % Settling time
Mp = 0.10;        % Percent oversoot
tr = 2.0;         % Rise time

a = 4.6/ts;         % Simplified variable for LMI
r = 1.8/tr;         % Simplified variable for LMI
c = -log(Mp)/pi;    % Simplified variable for LMI

% LMI for D-Stabilization (Lecture 5 Lemma 25)
Z = sdpvar(1,3);    % Define Z
P = sdpvar(3,3);    % Define P
X = sdpvar(3,3);    % Define P

%LMI1 = [Z>=eta]; 
LMI1 = [P>=stricly_greater0];
LMI2 = [[(-r*P) (A*P + B*Z); 
    (A*P + B*Z)' (-r*P)] <= 0];                 %rise time constraint

LMI3 = [A*P + B*Z + (A*P+B*Z)' + 2*a*P <= 0];   %settling time constraint

LMI4 = [[(A*P + B*Z + (A*P + B*Z)') c*(A*P + B*Z + (A*P + B*Z)'); 
    c*((A*P + B*Z)' - (A*P + B*Z)) (A*P + B*Z + (A*P + B*Z)')] <=0]; % overshoot constraint

F = [LMI1;LMI2;LMI3;LMI4];      % Combine into single constraint function
optimize(F,empty,ops);          % Run the optimization
Z_sol = value(Z);               % Return feasable solution Z
P_sol = value(P);               % Return feasable solution P
K=Z_sol*inv(P_sol);             % Construct K

disp("|   The calculated value of K=")
disp(K)
disp("|   and the eigan values of A+BK=")
K_eigs = eigs(A+B*K);

disp(K_eigs)
if K_eigs<=0
    disp("|   Since eigs(A+BK)<0, A+BK is Hurwitz")
else
    disp("|   Since eigs(A+BK)>0, A+BK is NOT Hurwitz")
end
fprintf("\n\n\n\n\n")

%% Problem 8.b
% What is the state-space representation of the system with the controller K? 
% What are the eigenvalues of this closed-loop system? 
% Do these eigenvalues make sense? Why/Why not? 
disp("HW2 Problem 8.b State Space Representation of Controller - Matlab Output")
% x_dot(t) = 
disp("|   applying u(t) = Kx(t) + v(t)")
disp("|   where K = state feedback gain & v = auxillery input")
disp("|   The state space representation changes from:")
disp("|      x'(t) = Ax(t) +Bu(t)")
disp("|      y(t) = Cx(t) + Du(t)")
disp("|   to:")
disp("|      x'(t) = Ax(t) + B*(Kx(t)+v(t))=(A+BK)x(t)+Bv(t)")
disp("|      y(t) = Cx(t) + D*(Kx(t)+v(t)) =(C+DK)x(t)+Dv(t)")
disp("|   and since D=0 in this system C+DK = C")
disp("|      y(t) = Cx(t)")
disp("|")
disp("|   Therefore the state space representation of this system is:")
disp("|      x'(t)=(A+BK)x(t)")
disp("|      y(t)=Cx(t)")
disp("|   where (A+BK) =")
disp(A+B*K)
fprintf('\n')
disp("|   and C =")
disp(C)
disp("|   Eigs of closed loop system:")
disp(eigs(A+B*K))
fprintf("\n\n\n\n\n")

%% Problem 8.c
% Plot the system's closed loop step response response using the controller K. 
% Estimate the rise-time and overshoot achieved. 
% Does it achieve the speci cations? 
disp("HW2 Problem 8.c Controller Step Response - Matlab Output")


sys = ss((A+B*K),B,C,D);            % Define System
hold off;
step(sys);                          % plot step input
title('Problem 8.c - Step Response of Feedback Controller') % Configure plot
disp("|   Step response plotted")
hold on;
step_info = stepinfo(sys) 
response_tr = step_info.RiseTime;
response_ts = step_info.SettlingTime;
response_os = step_info.Overshoot;


disp("|   The feedback controller:")
if response_ts < ts
    disp("|   did achieve the sepcification for Settling Time < 10s ")
else
    disp("|   did NOT achieve the sepcification for Settling Time < 10s ")
end
if response_os < Mp
    %fprintf("did achieve the sepcification for Overshoot since"
    disp("|   did achieve the sepcification for Overshoot < 10% ")
else
    disp("|   did NOT achieve the sepcification for Overshoot < 10%")
end
if response_tr > tr
    disp("|   did achieve the sepcification for Rise Time > 2s")
else
    disp("|   did NOT achieve the sepcification for Rise Time > 2s")
end
fprintf("\n\n\n\n\n")
%% HW 2 Problem 9 Main Function Script

% Use MOSEK as Solver
ops = sdpsettings('solver','mosek','verbose',0);
empty=[];

% Variable Definitions
A=[ 0 1 0;-30 0 20;2000 0 -2000];
B=[0;0.1;0];
C=[0 0 1];
D = 0;

A_T = transpose(A);
B_T = transpose(B);
BB_T = B*B_T;

n=length(A);                        % Get the size of matrix A
eta=0.00001;                       % Constant to allow for stricly > calculation during P>0
stricly_greater0 =eta*eye(n);       % Variable to estimate strictly zero constraint
strict_less0 = -eta*eye(n);
%% Problem 9.a
% Use an LMI to design an observer so that the closed-loop system achieves a 
% settling time of 10s or less, an overshoot of 10% or less and a rise-time 
% of 2s or greater. Print the observer L and the eigenvalues of A+LC. 
disp("HW2 Problem 9.a Feedback Controller Design - Matlab Output")

ts = 10.0;        % Settling time
Mp = 10.0;        % Percent oversoot
tr = 2.0;         % Rise time

a = 4.6/ts;         % Simplified variable for LMI
r = 1.8/tr;         % Simplified variable for LMI
c = -log(Mp)/pi;    % Simplified variable for LMI

% Lecture 6 Lemma 14
X = sdpvar(n);
Z = sdpvar(3,1);    % Define Z
P = sdpvar(n);
LMI1 = [
    [-r*P (P*A + Z*C)';
    (P*A + Z*C)   -r*P] <=0  ];

LMI2 = [(P*A + Z*C)' + P*A + Z*C + 2*a*P <=0 ];
LMI3 = [
    [((P*A + Z*C)'+P*A + Z*C) c*((P*A + Z*C)' - (P*A + Z*C));
    c*(P*A + Z*C - (P*A + Z*C)') ((P*A + Z*C)' + P*A + Z*C) ] <= 0];


F = [LMI1;LMI2;LMI3];
optimize(F,empty,ops);              % Run the optimization
X_sol = value(X);
Z_sol = value(Z);
P_sol = value(P);

L = inv(P_sol)*Z_sol;

disp("|   Observer L = ")
disp(L)
disp("|   Eigan Values of A+LC")
disp(eigs(A+L*C))

if real(eigs(A+L*C))<=0
    disp("|   Since Re(eigs(A+LC))<0, A+LC is Hurwitz")
else
    disp("|   Since Re(eigs(A+LC))>0, A+LC is NOT Hurwitz")
end
fprintf("\n\n\n\n\n")
%% Problem 9.b
% Plot the system's closed loop step response response using the observer L. 
% Estimate the rise-time and overshoot achieved. 
% Does it achieve the specifications? 
disp("HW2 Problem 9.b Observer Based Closed Loop Step Response - Matlab Output")
At = [ A+B*K -B*K
      zeros(size(A)) A+L*C ];
Bt = [ B;
      zeros(size(B))];
Ct = [ C    zeros(size(C)) ];
sys = ss(At,Bt,Ct,D);

step(sys);                          % plot step input
%title('Problem 9.b - Step Response of Observer Controller') % Configure plot
disp("|   Step response plotted")

step_info = stepinfo(sys)       
response_tr = step_info.RiseTime;
response_ts = step_info.SettlingTime;
response_os = step_info.Overshoot;

% a settling time of 10s or less, an overshoot of 10% or less and a rise-time 
% of 2s or greater.
disp("The feedback controller:")
if response_ts < ts
    disp("did achieve the sepcification for Settling Time < 10s ")
else
    disp("did NOT achieve the sepcification for Settling Time < 10s ")
end
if response_os < Mp
    %fprintf("did achieve the sepcification for Overshoot since"
    disp("did achieve the sepcification for Overshoot < 10% ")
else
    disp("did NOT achieve the sepcification for Overshoot < 10%")
end
if response_tr > tr
    disp("did achieve the sepcification for Rise Time > 2s")
else
    disp("did NOT achieve the sepcification for Rise Time > 2s")
end
fprintf("\n\n\n\n\n")
%% Problem 9.c
% If you were to create a Observer-based controller, what would you predict 
% the eigenvalues of the closed loop system would be? Why? 
disp("HW2 Problem 9.c Observer-based controller Eigan Values - Matlab Output")
disp("|   I would expect that the real part of the eigan values of the Observer-based ")
disp("|   controller would be negative.  This is because designing an observer ")
disp("|   requires that the observer dynamics are Hurwitz so that the error dynamics converge.")
disp("|   Also, because of this condition, we specify constraints such that observer gain L makes A+LC hurwitz ")
disp("|   during controller synthesis. However, an observer-based controller implies that a ")
disp("|   feedback gain K must also be specified which creates another condition for system")
disp("|   stability: A+BK is Hurwitz. If both A+BK and A+LC are stable then this implies closed-loop")
disp("|   system is stable due to the principle of separation of estimation and control in linear systems.")



