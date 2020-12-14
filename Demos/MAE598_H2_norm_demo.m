% Demo from M. Peet for lecture 11 in MAE 598


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




