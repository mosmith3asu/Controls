%% Problem 3
% Definition of variables
clear all
clc
options = sdpsettings('solver','mosek','verbose',0);
Ks=200;
g=9.8;
mh=60;
Bh=2.5;
Ie=0.36;
Be=0.607;
%Bc=0.01*Kc;
kf=0.034;
l=4;
%M=[Ih 0; 0 Ie]

%% P3.a
disp("##########################")
disp("##########3.a#############")
disp("##########################")
% Define given variables
lh=1;
nm=0.1;
nc=1;
nk=0.1;

tvp_Kc=[1 1000];
tvp_mh=[60 100];
b=0.01;
%tvp_Bc=-b*tvp_Kc;



a1=[1.593	1.485	1.269	1.13	0.896	0.559	0.398];

P = sdpvar(4,4);
Y = sdpvar(4,4);
Z = sdpvar(1,2);
gamma = sdpvar(1);

eps = 1e-5;

F = [];
N=size(tvp_Kc,2);
for i=1:1:N
 
    Kc=[0 0; 0 tvp_Kc(i)];
    Bc=-b*Kc;
    mh=tvp_mh(i);
%     
%     Th=0;
%     Te=0;
    Ih=mh*l^2;
    
    M=[Ih 0; 0 Ie]
    B=[Bh 0; 0 Be]
    K=[Ks -Ks; -Ks Ks];
    %G=[-mh*g*l*sin(q);0];
    G=[-mh*g*l 0;0 0];
    
    A0=[zeros(2) ones(2); -inv(M)*K+M*G -inv(M)*B]
    Avar=[zeros(2) zeros(2); -inv(M)*Kc -inv(M)*Bc]
    
    A=A0+Avar;
    
    F=[F;A'*P+A*P<=eps*eye(4)];
%     B2=[0 0; Ks/Ih 0]
%     B1=[0; Th/Ih]
%     C1=[0 0; Ks/Ie 0];
%     C2=[1 0];
%     D11=[0 1; -Ks/Ie -Be/Ie]
%     D12=[0; Te/Ie];
%     D21=[0 0];
%     D22=[0];
%     
%     M11 = Y*A' + A*Y + Z'*B2' + B2*Z;
%     M21 = B1';
%     M31 = C*Y + D12*Z;
%     M32 = D11;
%     M22 = -gamma*eye(2);
%     M33 = -gamma*eye(2);
%     M = [M11 M21' M31';
%         M21 M22 M32';
%         M31 M32 M33];
%     F = [F; M <= -eps*eye(size(M))];
    
    
end   
% for i=1:1:N
% A=[0 1/(lh*m0(i)^-1);-k0(i) -c0(i)];
% B1 = [-nm 0 0;0 nk nc];
% B2 = [0;1]; 
% C1 = [0 1/(lh*m0(i)); -k0(i) 0 ; 0 -c0(i)];
% C2 = [1 0];
% C= C1;
% D12 = [0;0;0];
% D11 = [-nm 0 0; 0 0 0; 0 0 0];
% 
% 
% M11 = Y*A' + A*Y + Z'*B2' + B2*Z;
% M21 = B1';
% M31 = C*Y + D12*Z;
% M32 = D11;
% M22 = -gamma*eye(2);
% M33 = -gamma*eye(2);
% M = [M11 M21' M31';
%     M21 M22 M32';
%     M31 M32 M33];
% F = [F; M <= -eps*eye(size(M))];
% end

F=[F;P>=eps*eye(size(Y))];
optimize(F,[],options);
P = value(P)
% F=[F;Y>=eps*eye(size(Y))];
% optimize(F,gamma,options);
% 
% gamma = value(gamma)
% 
% %Controller
% K = value(Z)*(value(Y'))


%% 3.b
% disp("##########################")
% disp("##########3.b#############")
% disp("##########################")
% disp("in writeup")
% %% 3.c
% disp("##########################")
% disp("##########3.c#############")
% disp("##########################")
% a1_prime = 0.2;
% a1=[1.593,0.398];
% a2=[260.559,51.003];
% a3=[185.488,53.84];
% a4=[1.506,0.421];
% a5=[0.298,0.078];
% 
% % form 2^k aggregate combination of all verticies
%  pre_a_aggr=combvec(a1,a2,a3,a4,a5)';
%  a_aggr=zeros(size(pre_a_aggr));
%  n=5; % number of inputs
%  for i=1:n
%      a_aggr(:,i)=pre_a_aggr(:,n-i+1);
%  end
% a_aggr = a_aggr';
% 
% sympref('FloatingPointOutput',true);
% Latex_out = latex(sym(a_aggr))
% 
% % Extract new verticies
% a1 = a_aggr(1,:);
% a2 = a_aggr(2,:);
% a3 = a_aggr(3,:);
% a4 = a_aggr(4,:);
% a5 = a_aggr(5,:);
% 
% 
% N=size(a1,2);
% for i=1:1:N
% A=[-a4(i) 1 -a5(i);((-a1_prime*a4(i))-a2(i)) (a1_prime-a1(i)) ((-a1_prime*a5(i))-a3(i));0 0 -(1/tz)];
% B1 = (wx/57.3)*[-1 0;-a1_prime JxJy_Jz; 0 0];
% B2 = [0;0;(1/tz)]; 
% C = (1/(57.3*g))*[(57.3*g) 0 0;V*a4(i) 0 V*a5(i)];
% D12 = [0;0]; 
% D11 = (1/(57.3*g))*[0 0; V*b7 0];
% 
% 
% M11 = Y*A' + A*Y + Z'*B2' + B2*Z;
% M21 = B1';
% M31 = C*Y + D12*Z;
% M32 = D11;
% M22 = -gamma*eye(2);
% M33 = -gamma*eye(2);
% M = [M11 M21' M31';
%     M21 M22 M32';
%     M31 M32 M33];
% F = [F; M <= -eps*eye(size(M))];
% end
% 
% 
% F=[F;Y>=eps*eye(size(Y))];
% optimize(F,gamma,options);
% 
% gamma = value(gamma)
% 
% K = value(Z)*(value(Y'))
% 


