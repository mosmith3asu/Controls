%% Problem 3
% Definition of variables
clear all
clc
options = sdpsettings('solver','mosek','verbose',0);


%% P3.a
disp("##########################")
disp("##########3.a#############")
disp("##########################")
% Define given variables
tz = 1;
JxJy_Jz = 0.7501;
g = 9.8;
V = 1401;
b7 = -.001;
wx = 10;

a1 = [1.593	1.485	1.269	1.13	0.896	0.559	0.398];
a1_prime= [0.285	0.192	0.147	0.118	0.069	0.055	0.043];
a2 = [260.559	266.415	196.737	137.385	129.201	66.338	51.003];
a3 = [185.488	182.532	176.932	160.894	138.591	78.404	53.84];
a4 = [1.506	1.295	1.169	1.13	1.061	0.599	0.421];
a5 = [0.298	0.243	0.217	0.191	0.165	0.105	0.078];


Y = sdpvar(3,3);
Z = sdpvar(1,3);
gamma = sdpvar(1);

eps = 1e-5;

F = [];
N=size(a1,2);
for i=1:1:N
A=[-a4(i) 1 -a5(i);
    ((-a1_prime(i)*a4(i))-a2(i)) (a1_prime(i)-a1(i)) ((-a1_prime(i)*a5(i))-a3(i));
    0 0 -(1/tz)];
B1 = (wx/57.3)*[
    -1 0;-a1_prime(i) JxJy_Jz; 
    0 0];
B2 = [0;0;(1/tz)]; 
C = (1/(57.3*g))*[(57.3*g) 0 0;V*a4(i) 0 V*a5(i)];
D12 = [0;0];
D11 = (1/(57.3*g))*[0 0; V*b7 0];


M11 = Y*A' + A*Y + Z'*B2' + B2*Z;
M21 = B1';
M31 = C*Y + D12*Z;
M32 = D11;
M22 = -gamma*eye(2);
M33 = -gamma*eye(2);
M = [M11 M21' M31';
    M21 M22 M32';
    M31 M32 M33];
F = [F; M <= -eps*eye(size(M))];
end


F=[F;Y>=eps*eye(size(Y))];
optimize(F,gamma,options);

gamma = value(gamma)

%Controller
K = value(Z)*(value(Y'))


%% 3.b
disp("##########################")
disp("##########3.b#############")
disp("##########################")
disp("in writeup")
%% 3.c
disp("##########################")
disp("##########3.c#############")
disp("##########################")
a1_prime = 0.2;
a1=[1.593,0.398];
a2=[260.559,51.003];
a3=[185.488,53.84];
a4=[1.506,0.421];
a5=[0.298,0.078];

% form 2^k aggregate combination of all verticies
 pre_a_aggr=combvec(a1,a2,a3,a4,a5)';
 a_aggr=zeros(size(pre_a_aggr));
 n=5; % number of inputs
 for i=1:n
     a_aggr(:,i)=pre_a_aggr(:,n-i+1);
 end
a_aggr = a_aggr';

sympref('FloatingPointOutput',true);
Latex_out = latex(sym(a_aggr))

% Extract new verticies
a1 = a_aggr(1,:);
a2 = a_aggr(2,:);
a3 = a_aggr(3,:);
a4 = a_aggr(4,:);
a5 = a_aggr(5,:);


N=size(a1,2);
for i=1:1:N
A=[-a4(i) 1 -a5(i);((-a1_prime*a4(i))-a2(i)) (a1_prime-a1(i)) ((-a1_prime*a5(i))-a3(i));0 0 -(1/tz)];
B1 = (wx/57.3)*[-1 0;-a1_prime JxJy_Jz; 0 0];
B2 = [0;0;(1/tz)]; 
C = (1/(57.3*g))*[(57.3*g) 0 0;V*a4(i) 0 V*a5(i)];
D12 = [0;0]; 
D11 = (1/(57.3*g))*[0 0; V*b7 0];


M11 = Y*A' + A*Y + Z'*B2' + B2*Z;
M21 = B1';
M31 = C*Y + D12*Z;
M32 = D11;
M22 = -gamma*eye(2);
M33 = -gamma*eye(2);
M = [M11 M21' M31';
    M21 M22 M32';
    M31 M32 M33];
F = [F; M <= -eps*eye(size(M))];
end


F=[F;Y>=eps*eye(size(Y))];
optimize(F,gamma,options);

gamma = value(gamma)

K = value(Z)*(value(Y'))



