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

%% Problem 7.a
% Use the Lyapunov LMI with Q = I to determine if the system is 
% open-loop stable. Cofirm this result by determining the 
% open loop poles of the system.

disp("HW2 Problem 7.a Stability Analysis - Matlab Output")

eta=1e-20;                   % Constant to allow for stricly > calculation during P>0
stricly_greater0 =eta*eye(n); 	% Variable to estimate strictly zero constraint

% LYAPUNOV LMI FOR SYSTEM STABILITY (Lecture 4 Lemma 8)
P = sdpvar(n);                  % Make P a variable (No extra parameters assumes symetric)  
LMI1 = [P>=stricly_greater0];   % First Lyapunov LMI (Canidate Lyapunov Function)
LMI2 = [A'*P + P*A <= 0];       % Second Lyapunov LMI
F =[LMI1; LMI2];                % Group constraints
optimize(F,empty,ops);          % Run the optimization
P_sol = value(P);               % Return P as varaible
%P_sol = round(P_sol,2);        % Round for errors

stable = P_sol>=stricly_greater0; % Test stability

% DISPLAY RESULTS
disp("|Lyapunov LMI Stability:") 
disp("|   A is Hurwitz if and only if there exists P>0 such that A_TP +PA<=0")
if stable
    disp("|   ANS: The system is stable since P>0")
else
    disp("|   ANS: The system is NOT stable since P<=0")
end
disp("|P Matrix")
disp(P_sol)

% POLE ANALYSIS FOR POLE STABILITY
% sys = ss(A,B,C,D);          % Define system
% poles = pole(sys);          % Return Poles
% real_poles = real(poles);   % Get real part of poles
poles = eigs(A);
real_poles = real(poles);

stable = real_poles < 0;        % Check stability condition

% DISPLAY RESULTS
disp('|Pole Analysis of Stability:')
disp('|   A is stable when the polse are located on the left plane')
disp('|   ( A is Hurwitz when Re(eigs(A))<0 )')
if stable
    disp('|   ANS: system is stable since the real part of poles <0')
else
    disp('|   ANS: system is unstable since the real part of poles >0')
end
disp("|   Poles:")
disp(poles)
fprintf("\n\n\n\n\n")

%% Problem 7.b
% Plot the step response of the open loop system. Does it have inherently
% nice properties or would it be useful to add a feedback controller to
% improve the response?
disp("HW2 Problem 7.b Step Response - Matlab Output")
sys = ss(A,B,C,D);          % Define system
step(sys);                  % Plot system step responses
title('Problem 7.b - Step Response of Open-Loop Controller')
disp("|Step response plotted")
disp("|Step Response Properties:")
disp("|   The system's settling time and overshoot are very large")
disp("|   and including a feedback controller would be helpful in improving")
disp("|   these properties in the step response")
fprintf("\n\n\n\n\n")


%% Problem 7.c
% Use an LMI to determine if the system is stabilizable. 
disp("HW2 Problem 7.c Stabilizability - Matlab Output")

eta=1e-3;                               % Constant to allow for stricly > calculation during P>0
stricly_greater0 =eta*eye(n);           % Variable to estimate strictly zero constraint


% STABALIZABILITY LMI (Lecture 5 Lemma 16)
X = sdpvar(3,3);                        % Create X variable (grammina)
gamma = sdpvar(3,3);                    % Create gamma variable
LMI1 = [X >= stricly_greater0];         % First LMI Constraint
LMI2 = [gamma >= stricly_greater0];     % Second LMI Constraint
LMI3 = [A*X + X*A' - gamma*B*B' <= 0];  % Third LMI Constraint
F=[LMI1;LMI2;LMI3];                     % Contain LMIs in single variable
 
optimize(F,empty,ops);                  % Run the optimization
X_sol = value(X);                       % Return X as varaible
%X_sol = round(X_sol,2);                % Round to eliminate calc errors
gamma_sol = value(gamma);
 
test1 = round(X_sol,0)>=stricly_greater0;               % Test initial LMI1 constraint
test2 = round(gamma_sol,0)>=stricly_greater0;           % Test initial LMI2 constraint

% DISPLAY RESULTS
disp("|Stabilizability of System:")
disp("|   (A,B) is stabilizable if and only if there exists a X>0,gamma>0")
disp("|   such that A X+X A^T - gamma BB^T < 0")

if all(test1(:)) && all(test2(:)) % Check if constraint is still satisfied
    disp("|   ANS: Since there exists X,gamma>0, (A,B) is stabilizable")
else
    disp("|   ANS: Since there DNE an X,gamma>0, (A,B) is unstabilizable")
end 

disp("|X")
disp(X_sol)
disp("|gamma")
disp(gamma_sol)
fprintf("\n\n\n\n\n")
%% Problem 7.d
% Determine if the system is controllable (it is not necessary to use an LMI for this part). 
disp("HW2 Problem 7.d Controllability - Matlab Output")

% CONTROLLABILITY MATRIX ANALYSIS (Definition 9 Lecture 6)
% The system (A,B) is controllable if C_{AB}= Im C(A,B) =R^n
CAB = [B A*B A^2*B];    % Define controllability matrix
%C_AB = image(CAB)      % How to find the image of CAB? ###############
CAB_rank = rank(CAB);   % Find out the rank of CAB
test1 = CAB_rank==n;    % Test controllable condition

%DISPLAY RESULTS
disp("Controllability of System")
disp("|   The system (A,B) is controllable if C_{AB}= Im C(A,B) =R^n")
if test1
    disp("|   ANS:The system is controllable since C_{AB}= Im C(A,B)= R^n")
else
    disp("|   ANS:The system is NOT controllable since C_{AB}= Im C(A,B)= R^n")
end
disp("|Controllability Matrix")
disp(CAB)
disp("|With rank:")
disp(CAB_rank)
fprintf("\n\n\n\n\n")

%% Problem 7.e
% Use an LMI to find the infinite time controllability grammian
disp("HW2 Problem 7.e Controllability Grammian - Matlab Output")

eta=1e-15;                   % Constant to allow for stricly > calculation during P>0
stricly_greater0 =eta*eye(n); 	% Variable to estimate strictly zero constraint

% CONTROLLABILITY GRAMIAN LMI (Lecture 5 Lemma 14)
W = sdpvar(3,3);                    % Make W a variable

LMI1 = [W >= stricly_greater0];     % First LMI Constraint
LMI2 = [A*W + W*A' + B*B'==0];      % Second LMI Constraint
F = [LMI1;LMI2];                    % Pack Constratints
optimize(F,empty,ops);              % Run the optimization
W_sol = value(W);                   % Return W as varaible

test1 = round(W_sol,2)>=stricly_greater0;    % Test initial LMI constraints

% DISPLAY RESULTS
disp("|   Controllability Grammian W")
disp("|   if (A,B) is controllable, W>0 is the unique solution to")
disp("|   AW+WA^T+BB^T = 0")
if test1
    disp("|   The controllable grammian W:")
else
    disp("|   unable to find grammian")
end
disp("|Controllability Grammian of (A,B)")
disp(W_sol)
%W_check = gram(sys,'c')
fprintf("\n\n\n\n\n")

%% Problem 7.f
% Use your controllability gramian to determine whether it will take more
% energy to reach the state xa = [ 0 1 1 ]T or xb = [ 1 0 0 ]T . 
% Carefully explain your conclusion. What is the minimum energy for each?
disp("HW2 Problem 7.f Controllability Grammian Energy - Matlab Output")
xa = [0; 1; 1];         % define state xa
xb = [1; 0; 0];         % define state xb
xa_T = transpose(xa);   % define state xa^T
xb_T = transpose(xb);   % define state xb^T
W_inv = inv(W_sol);     % define inverse of grammian

% FIND THE INPUT SIZE OF EACH STATES
size_u_xa = sqrt(xa_T*W_inv*xa);  % calculate size of input for state xa
size_u_xb = sqrt(xb_T*W_inv*xb);  % calculate size of input for state xb

% DISPLAY RESULTS
disp("|Energy to reac state xa,xb:")
disp("|   The size of input (u) to reach state xd is equal to")
disp("|   ||u||^2_L2 = xd^T *W^-1 *xd")  
disp("|   this quantity refers to the energy needed by an input, which can be expressed")  
disp("|   as the magnitude of the input to reach state xd in a fixed ammount of time.")  
disp("|   If a state is harder to reach, it will take a larger input to drive the system to that state")  
disp("|   in the same ammount fo time as a state that is easire to reach.")  

if size_u_xa>size_u_xb
    disp("|   ANS: state xa requires more energy to reach since the size of the input")
    disp("|   ANS: is larger than xb ")
else 
    disp("|   ANS: state xb requires more energy to reach since the size of the input")
    disp("|   ANS: is larger than xa ")
end
disp("|the size of each input is as follows")
disp("|xa:")
disp(size_u_xa)
disp("|xb:")
disp(size_u_xb)
fprintf("\n\n\n\n\n")
%% Problem 7.g
% Use an LMI to determine if the system is detectable. 
disp("HW2 Problem 7.g Detectability - Matlab Output")

eta=1e-3;                       % Constant to allow for stricly > calculation during P>0
stricly_greater0 =eta*eye(n); 	% Variable to estimate strictly zero constraint

% Theorem 11.An observer exists if and only if(C,A)is detectable
disp("|(A,B) is detectable if and only if there exists X>0 such that")
disp("|  AX+XA'-C'C<0")
X = sdpvar(3,3);                            % Make W a variable
LMI1 = [X >= stricly_greater0];             % First LMI Constraint
LMI2 = [A*X + X*A' - C'*C<=0];              % Second LMI Constraint
F = [LMI1;LMI2];
optimize(F,empty,ops);                      % Run the optimization
X_sol = value(X)
disp("|Since we were able to find X>0, (A,B) is detectable")
fprintf("\n\n\n\n\n")

%% Problem 7.h
% Use an LMI to find the infinite time observability grammian. 
% OBSERVABILITY GRAMMIAN (Lecture 6 Lemma 8)
%(C,A) is observable iff Y>0 is the unique solution to A^TY+YA+C^TC=0
disp("HW2 Problem 7.h Observability Grammian - Matlab Output")

eta=1e-15;                          % Constant to allow for stricly > calculation during P>0
stricly_greater0 =eta*eye(n);       % Variable to estimate strictly zero constraint

% Lecture 6 Lemma 8
Y = sdpvar(3,3);                    % Make W a variable

LMI1 = [Y >= stricly_greater0];     % First LMI Constraint
LMI2 = [A'*Y + Y*A + C'*C==0];      % Second LMI Constraint
F = [LMI1;LMI2];

optimize(F,empty,ops);              % Run the optimization
Y_sol = value(Y);                   % Return W as varaible


% DISPLAY RESULTS
disp("|   Observability Grammian Y")
disp("|   if (A,B) is controllable, W>0 is the unique solution to")
disp("|   A^TY+YA+C^TC = 0")

disp("|Controllability Grammian of (A,B)")
disp(Y_sol) 
%Y_check = gram(sys,'o')              % Check result with matlab sol
fprintf("\n\n\n\n\n")
%% Problem 7.i
% Use your observability grammian to determine which of the initial 
% conditions,xa0 = [ 0 1 1 ]T or xb0 = [ 1 0 0 ]T , will produce 
% the larger output. Carefully explain your conclusion. What is the 
% size of the output for each? 
disp("HW2 Problem 7.i Observability Grammian Energy - Matlab Output")
xa = [0; 1; 1];         % define state xa
xb = [1; 0; 0];         % define state xb
xa_T = transpose(xa);   % define state xa^T
xb_T = transpose(xb);   % define state xb^T

% FIND THE INPUT SIZE OF EACH STATES
size_y_xa = sqrt(xa_T*Y_sol*xa);  % calculate size of output for state xa
size_y_xb = sqrt(xb_T*Y_sol*xb);  % calculate size of output for state xb

% DISPLAY RESULTS
disp("|Size of output xa,xb:")
disp("|   The size of output xd is equal to")
disp("|   ||y||^2_L2 = xd^T*Y*xd")  
disp("   this quantity decribes what conditions (xd) will produce what size of outputs from the system")
disp("   However, states that produce a larger output are more weakly observable states.") 

if size_y_xa>size_y_xb
    disp("|   ANS: state xb will produce the largest output since the ellipsoid")
    disp("|   ANS: EXPLAINATION GIVEN ON LECTURE 6 SLIDE 4 ")
else 
    disp("|   ANS: state xb will produce the largest output since the ellipsoid")
    disp("|   ANS: EXPLAINATION GIVEN ON LECTURE 6 SLIDE 4 ")
end
disp("|the size of each input is as follows")
disp("|xa:")
disp(size_y_xa)
disp("|xb:")
disp(size_y_xb)
