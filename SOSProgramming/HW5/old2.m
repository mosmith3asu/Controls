% pvar x y
% p=4*x^4+4*x^3*y-7*x^2*y^2-2*x*y^3+10*y^4;
% prog=sosprogram([x y]);
% prog=sosineq(prog,p);
% prog=sossolve(prog);
clear all 
clc

%% Global Stability
pvar x y gam
% =============================================
% First, initialize the sum of squares program
prog=sosprogram([x y]);
% =============================================
% Declare decision variable gam too
prog = sosdecvar(prog,[gam]);
% =============================================
% Next, define SOSP constraints
f=x^4+y^4-2*y*x^3-3*y^2*x^2+150*x^2+150*y^2;
% f=[4*x^3-6*y*x^2-6*y^2*x+300*x;
%     4*y^3-2*x^3-6*y*x^2+300*y];
% prog = sosineq(prog,x^2-12^2);
% prog = sosineq(prog,y^2-12^2);
prog = sosineq(prog,(f-gam),[-12,12;-12,12]);
% =============================================
% Set objective : maximize gam
prog = sossetobj(prog,-gam);
% =============================================
% And call solver
prog = sossolve(prog);
% =============================================
% Finally, get solution
SOLgamma = sosgetsol(prog,gam)
SOLx = sosgetsol(prog,x)
%%  SOSTOOLS Code:Global Stability
pvar x y gam
f=[4*x^3-6*y*x^2-6*y^2*x+300*x;
    4*y^3-2*x^3-6*y*x^2+300*y];
prog=sosprogram([x y]);
% =============================================
% The Lyapunov function V(x):
%[prog,V] = sospolyvar(prog,[x^2; y^2],'wscoeff');

Z=monomials([x,y],0:4);
[prog,V]=sossosvar(prog,Z);
% =============================================
% Next, define SOSP constraints
V=V+.0001*(x^2+y^2);
% Constraint 1 : V(x) - (x1^2 + x2^2 + x3^2) >= 0
prog = sosineq(prog,V-(x^2+y^2));
prog = sosineq(prog,x^2-12^2);
prog=soseq(prog,subs(V,[x; y],[0; 0]));
% Constraint 2: -dV/dx*(x3^2+1)*f >= 0
nablaV=[diff(V,x);diff(V,y)];
prog=sosineq(prog,-nablaV'*f);
% =============================================
% And call solver
prog=sossolve(prog);
Vn=sosgetsol(prog,V)
Vn=sosgetsol(prog,[x,y])
%% SOSTOOLS Code:Find a Local Lyapunov Function
% pvar x y
% mu=1; r=2.8;
% g=r?(x2+y2);
% f= [?y;?mu?(1?x2)?y+x];
% prog=sosprogram([x y]);
% Z2=monomials([x y],0:2);
% Z4=monomials([x y],0:4);
% [prog,V]=sossosvar(prog,Z2);
% V=V+.0001?(x4+y4);
% prog=soseq(prog,subs(V,[x, y]’,[0, 0]’));
% nablaV=[diff(V,x);diff(V,y)];
% [prog,s]=sossosvar(prog,Z2);
% prog=sosineq(prog,-nablaV’*f-s*g);
% prog=sossolve(prog);
% Vn=sosgetsol(prog,V)