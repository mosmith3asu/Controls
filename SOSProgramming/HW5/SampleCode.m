pvar x y gam
prog=sosprogram([x y]);
% f= [-y-1.5*x^2-.5*x^3;3*x-y];
f=x^4+y^4-2*y*x^3-3*y^2*x^2+150*x^2+150*y^2;
%f=[4*x^3-6*y*x^2-6*y^2*x+300*x; 4*y^3-2*x^3-6*y*x^2+300*y];



[prog,V]=sossosvar(prog,Z);
V=V+.0001*(x^4+y^4);
nablaV=[diff(V,x);diff(V,y)];
prog=sosineq(prog,-nablaV'*f+gam);
%prog=soseq(prog,subs(V,[x; y],[0; 0]));


prog=sosineq(prog,-nablaV'*f);

prog=sossolve(prog);
Vn=sosgetsol(prog,V)


% 
% pvar x y gam;
% prog=sosprogram([x y]);
% prog = sosdecvar(prog,[gam]);
% 
% f=x^4+y^4-2*y*x^3-3*y^2*x^2+150*x^2+150*y^2;
% g=50000;
% prog=sosineq(prog,(f-gam),[-12,12;-12,12]); 
% prog = sossetobj(prog,[x,y]);
% 
% 
% prog=sossolve(prog);
% Gamma=sosgetsol(prog,x)