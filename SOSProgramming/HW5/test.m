clear all
clc
pvar x1 x2 x3 x4 x5;
vartable = [x1; x2; x3; x4; x5];
prog = sosprogram(vartable);

gamma = 50000;
f = 2.5 - .5*x1*x2 - .5*x2*x3 - .5*x3*x4 - .5*x4*x5 - .5*x5*x1;
bc1 = x1^2 - 1 ;
bc2 = x2^2 - 1 ;
bc3 = x3^2 - 1 ;
bc4 = x4^2 - 1 ;
bc5 = x5^2 - 1 ;
Z = monomials(vartable,0);
for i = 1:5
    [prog, p{1+i}]= sospolyvar(prog,Z);
end
expr = (gamma-f)+p{1}*bc1+p{2}*bc2+p{3}*bc3+p{4}*bc4+p{5}*bc5;
prog = sosineq(prog,expr);
prog = sossolve(prog);