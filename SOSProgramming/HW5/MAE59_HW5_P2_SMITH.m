clear all
clc
echo off

%%  SOSTOOLS Code:Global Stability

% =============================================
% First, initialize the sum of squares program
pvar x y gam


prog=sosprogram([x y]);
% =============================================
% Declare decision variable gam too
prog = sosdecvar(prog,[gam])
% =============================================
% Next, define SOSP constraints
% f=x^4+y^4-2*y*x^3-3*y^2*x^2+150*x^2+150*y^2;
f=[4*x^3-6*y*x^2-6*y^2*x+300*x;
    4*y^3-2*x^3-6*y*x^2+300*y];
%prog = sosineq(prog,(f-gam),[-12,12;-12,12]);
%prog = sossetobj(prog,-gam);
%prog = sossolve(prog);
%SOLgamma = sosgetsol(prog,gam)

Z=monomials([x,y],0:8);
[prog,V]=sossosvar(prog,Z);
V=V+.0001*(x^4+y^4);
prog=soseq(prog,subs(V,[x; y],[0; 0]));
% prog = sosineq(prog,(f-gam));
nablaV=[diff(V,x);diff(V,y)];
prog=sosineq(prog,-nablaV'*f);
% % =============================================
% % The Lyapunov function V(x):
% %[prog,V] = sospolyvar(prog,[x^2; y^2],'wscoeff');
% 
% 
% % =============================================
% % Next, define SOSP constraints
% 
% % Constraint 1 : V(x) - (x1^2 + x2^2 + x3^2) >= 0
% %prog = sosineq(prog,(f-gam),[-12,12;-12,12]);
% %prog = sosineq(prog,V-(x^2+y^2));
% %prog = sosineq(prog,x^2-12^2);
% %prog=soseq(prog,subs(V,[x; y],[0; 0]));
% % Constraint 2: -dV/dx*(x3^2+1)*f >= 0
% 
% % =============================================
% % And call solver
prog=sossolve(prog);
Vn=sosgetsol(prog,V)
% XY=sosgetsol(prog,x,y)