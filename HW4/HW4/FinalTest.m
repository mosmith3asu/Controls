clear all

%m0=3;c0=1;k0=2;etam=.1;etac=.2;etak=.4;
m0=1;c0=.1;k0=1;etam=.1;etac=1;etak=.2;
%m0=3;c0=1;k0=2;etam=.4;etac=.2;etak=.4; % values from paper
% Open-Loop Formulation
A=[0 1;-k0/m0 -c0/m0];
B1=[0 0 0;
    -etam -etac/m0 -etak/m0]
B2=[0;
    1/m0];
C1=[-k0/m0 -c0/m0
    0 c0;
    k0 0];
C2=[1 0];
D11=[-etam -etac/m0 -etak/m0;
    0 0 0;
    0 0 0]
D12=[1/m0;
    0;
    0];
D21=[0 0 0];
D22=[0];

eta=.001
n=2; % number of states
np=3; % number of uncertain outputs
nq=3; % number of uncertain inputs

% Declare Scalings
th1=sdpvar(1);
th2=sdpvar(1);
th3=sdpvar(1);
Th=[th1 0 0;0 th2 0; 0 0 th3]
%Th=diag([th1;th2;th3],0)
%Th=blkdiag(blkdiag(th1,th2),th3)
F=[];
F=[F;Th>=0];

% Lyap variable
tol=.1
gam_u=100;
gam_l=0
X=sdpvar(n);
F=[F;X>=eta*eye(n)]
%B=B1;C=C
gam_new=gam_u;
err=gam_u
while err>tol
    MAT=[A'*X+X*A X*B1;B1'*X -Th]+1/gam_new/gam_new*[C1 D11]'*Th*[C1 D11];
    Ftemp=[F;MAT<=0]
    DIAG=optimize(Ftemp)
    if DIAG.problem==0
        gam_u=gam_new
    else
        gam_l=gam_new
    end
    gam_new=(gam_u+gam_l)/2
    err=gam_u-gam_l
end

% Stabilizing Controller Synthesis

% Controlled Formulation
A=[0 1;-k0/m0 -c0/m0];
B1=[0 0 0;
    -etam -etac/m0 -etak/m0]
B2=[0;
    1/m0];
B3=[0;
    1/m0];
C1=[-k0 -c0/m0
    0 c0;
    k0 0];
C2=[1 0;
    0 0];
C3=[1 0];
D11=[-etam -etac/m0 -etak/m0;
    0 0 0;
    0 0 0]
D12=[1/m0;
    0;
    0];
D13=[0;
    0;
    0];
D21=[0 0 0;
    0 0 0];
D22=[0;
    0];
D23=[0;
    1];
D31=[0 0 0];
D32=0;
D33=0;

eta=.001
n=2; % number of states
np=3; % number of uncertain outputs
nq=3; % number of uncertain inputs
nc=size(D13,2);

% Declare Scalings
th1=sdpvar(1);
th2=sdpvar(1);
th3=sdpvar(1);
Th=[th1 0 0;0 th2 0; 0 0 th3]
F=[];
%F=[F;Th>=0];
X=sdpvar(n);
Z=sdpvar(nc,n,'full')
F=[F;X>=eta*eye(n)]
MAT=[A*X+B3*Z+X*A'+Z'*B3' X*C1'+Z'*D13';C1*X+D13*Z -Th]+[B1; D11]*Th*[B1' D11'];
    F=[F;MAT<=0]
    DIAG=optimize(F)
Zn=value(Z)
Xn=value(X)
K=Zn*inv(Xn)

% Part c - Close the loop with stabilizing controller and recalculate mu

Acl=A+B3*K; C1cl=C1+D12*K;
% Declare Scalings
F=[];
%F=[F;Th>=0];

% Lyap variable
tol=.01
gam_u=100;
gam_l=0
gam_new=gam_u;
err=gam_u
while err>tol
    clear X F MAT th1 th2 th3 Th
th1=sdpvar(1);
th2=sdpvar(1);
th3=sdpvar(1);
    Th=[th1 0 0;0 th2 0; 0 0 th3];
    X=sdpvar(n);
    F=[];
    F=[F;X>=eta*eye(n)];
    MAT=[Acl'*X+X*Acl X*B1;B1'*X -Th]+1/gam_new/gam_new*[C1cl D11]'*Th*[C1cl D11];
    Ftemp=[F;MAT<=0]
    DIAG=optimize(Ftemp)
    if DIAG.problem==0
        gam_u=gam_new
    else
        gam_l=gam_new
    end
    gam_new=(gam_u+gam_l)/2
    err=gam_u-gam_l
end


% Part d- Now do Hinfinity optimal robust control

nro=size(C2,1);

% Declare Scalings
gamma=sdpvar(1);
th1=sdpvar(1);
th2=sdpvar(1);
th3=sdpvar(1);
Th=[th1 0 0;0 th2 0; 0 0 th3]
F=[];
%F=[F;Th>=0];
X=sdpvar(n);
Z=sdpvar(nc,n,'full')
F=[F;X>=eta*eye(n)]
MAT=[A*X+B3*Z+X*A'+Z'*B3'+B1*Th*B1' (C2*X+D23*Z)'    X*C1'+Z'*D13'+B1*Th*D11' ;
    C2*X+D23*Z            -gamma*eye(nro)  zeros(nro,np); 
    C1*X+D13*Z+D11*Th*B1'             zeros(np,nro) D11*Th*D11'-Th];
    F=[F;MAT<=0]
    DIAG=optimize(F,gamma)
Zn=value(Z)
Xn=value(X)
K=Zn*inv(Xn)
Thn=value(Th)

value(gamma)
% Part d- Now do Hinfinity optimal robust control

nro=size(C2,1);

% Declare Scalings
gamma=sdpvar(1);
th1=sdpvar(1);
th2=sdpvar(1);
th3=sdpvar(1);
Th=[th1 0 0;0 th2 0; 0 0 th3]
F=[];
F=[F;Th>=0];
X=sdpvar(n);
Z=sdpvar(nc,n,'full')
F=[F;X>=eta*eye(n)]
MAT=[A*X+B3*Z+X*A'+Z'*B3'+B1*Th*B1' (C2*X+D23*Z)'    X*C1'+Z'*D13' ;
    C2*X+D23*Z            -gamma*eye(nro)  zeros(nro,np); 
    C1*X+D13*Z             zeros(np,nro) -Th];
    F=[F;MAT<=0]
    DIAG=optimize(F,gamma)
Zn=value(Z)
Xn=value(X)
K=Zn*inv(Xn)
Thn=value(Th)

value(gamma)


%Thn=eye(3);



% % Now do D-K iteration
% 
% % Construct System to Control with nominal scalings
% Thg=sqrt(blkdiag(Thn,eye(nro))); Thgr=sqrt(blkdiag(Thn,eye(1))) 
% Adk=A;B1dk=[B1 B2]*inv(Thgr); B2dk=B3;
% C1dk=Thg*[C1;C2]; C2dk=C3; D11dk=Thg*[D11 D12;D21 D22]*inv(Thgr); D12dk=Thg*[D13;D23]; D21dk=[D31 D32]*inv(Thgr); D22dk=D33;
% 
% % Now compute the Hinfinity optimal dynamic-feedback controller
% 
% eps=.000001;    % degree of strict positivity   
% ns=size(Adk,1);   % number of states
% nc=size(B2dk,2);  % number of actuators
% nm=size(C2dk,1);  % number of sensors
% nd=size(B1dk,2);  % number of external inputs
% no=size(C1dk,1);  % number of regulated outputs
% 
% 
% % H-infinity Dynamic Output Feedback Controller Synthesis
% 
% 
% % Declare the variables
% gamma=sdpvar(1);               % represents the bound on the H-infinity norm of the CL system.
% X1=sdpvar(ns);
% Y1=sdpvar(ns);
% An=sdpvar(ns,ns,'full')
% Cn=sdpvar(nc,ns,'full');
% Dn=sdpvar(nc,nm,'full');
% Bn=sdpvar(ns,nm,'full');
% 
% % declare constraints
% F=[X1>=eps*eye(ns)]
% F=[F;Y1>=eps*eye(ns)]
% F=[F;[X1 eye(ns); eye(ns) Y1]>=0];
% MAT=[Adk*Y1+Y1*Adk'+B2dk*Cn+Cn'*B2dk'  (Adk'+An+(B2dk*Dn*C2dk)')'        B1dk+B2dk*Dn*D21dk           (C1dk*Y1+D12dk*Cn)'; 
%      Adk'+An+(B2dk*Dn*C2dk)'         X1*Adk+Adk'*X1+Bn*C2dk+C2dk'*Bn'    X1*B1dk+Bn*D21dk           (C1dk+D12dk*Dn*C2dk)'  ;
%      (B1dk+B2dk*Dn*D21dk)'           (X1*B1dk+Bn*D21dk)'             -gamma*eye(nd)          (D11dk+D12dk*Dn*D21dk)'  ;
%      C1dk*Y1+D12dk*Cn              C1dk+D12dk*Dn*C2dk                D11dk+D12dk*Dn*D21dk         -gamma*eye(no)];
% 
% F=[F;MAT<=0];
% OPTIONS = sdpsettings('solver','sedumi')
% 
% % Solve the LMI, minimizing gamma
% optimize(F,gamma,OPTIONS)
% gamman=value(gamma)
% 
% 
% % retrieve decision variables
% X1n=value(X1); Y1n=value(Y1); Ann=value(An);Bnn=value(Bn);Cnn=value(Cn);Dnn=value(Dn);
% temp1=[Ann Bnn; Cnn Dnn]-[X1n*Adk*Y1n zeros(ns,nm); zeros(nc,ns) zeros(nc,nm)];
% 
% % Choose X2, Y2, so that X2*Y2=I-X1*Y1;
% Y2n=eye(ns);X2n=eye(ns)-X1n*Y1n;
% 
% % Reverse variable substitution
% temp2=inv([X2n X1n*B2dk;zeros(nc,ns) eye(nc)])*temp1*inv([Y2n' zeros(ns,nm); C2dk*Y1n eye(nm)]);
% Ak2=temp2(1:ns,1:ns);Bk2=temp2(1:ns,(ns+1):(ns+nm));Ck2=temp2((ns+1):(ns+nc), 1:ns);Dk2=temp2((ns+1):(ns+nc), (ns+1):(ns+nm));
% Dk=inv(eye(nc)-Dk2*D22dk)*Dk2;
% Bk=Bk2*(eye(nm)-D22dk*Dk);
% Ck=(eye(nc)-Dk*D22dk)*Ck2;
% Ak=Ak2-Bk*inv(eye(nm)-D22dk*Dk)*D22dk*Ck;
% 
% % Close the loop with Lower LFT
% plant=ss(Adk,[B1dk B2dk],[C1dk;C2dk],[D11dk D12dk; D21dk D22dk]);
% controller=ss(Ak,Bk,Ck,Dk);
% sys_cl=lft(plant,controller);
% 
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % now reformulate the upper plant and find the scalings which minimize mu
% 
% [Ath,Bth,Cth,Dth] = ssdata(sys_cl);
% 
% % Declare Scalings
% th1=sdpvar(1);
% th2=sdpvar(1);
% th3=sdpvar(1);
% Th=[th1 0 0;0 th2 0; 0 0 th3];
% %Thr=blkdiag(Th,eye(2))
% Thl=sqrt(blkdiag(Thn,eye(nro))); Thr=sqrt(blkdiag(Thn,eye(1))) ;
% %Th=diag([th1;th2;th3],0)
% %Th=blkdiag(blkdiag(th1,th2),th3)
% F=[];
% F=[F;Th>=0];
% 
% % Lyap variable
% tol=.01;
% gam_u=100;
% gam_l=0;
% ntemp=size(Ath,1);
% X=sdpvar(ntemp);
% F=[F;X>=eta*eye(ntemp)];
% %B=B1;C=C
% gam_new=gam_u;
% err=gam_u;
% while err>tol
%     MAT=[Ath'*X+X*Ath X*Bth;Bth'*X -Thr]+1/gam_new/gam_new*[Cth Dth]'*Thl*[Cth Dth];
%     Ftemp=[F;MAT<=0]
%     DIAG=optimize(Ftemp)
%     if DIAG.problem==0
%         gam_u=gam_new
%     else
%         gam_l=gam_new
%     end
%     gam_new=(gam_u+gam_l)/2
%     err=gam_u-gam_l
% end
% 
% Thn=sqrt(value(Th));
% %Thg=sqrt(blkdiag(Thn,eye(nro))); Thgr=sqrt(blkdiag(Thn,eye(1))) 
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Retry D-K iteration using only uncertain outputs
% Thg=sqrt(Thn)
% % Construct System to Control with nominal scalings
% Adk=A;B1dk=[B1]*inv(Thg); B2dk=B3;
% C1dk=Thg*[C1]; C2dk=C3; D11dk=Thg*[D11]*inv(Thg); D12dk=Thg*[D13]; D21dk=[D31]*inv(Thg); D22dk=D33;
% 
% % Now compute the Hinfinity optimal dynamic-feedback controller
% 
% eps=.000001;    % degree of strict positivity   
% ns=size(Adk,1);   % number of states
% nc=size(B2dk,2);  % number of actuators
% nm=size(C2dk,1);  % number of sensors
% nd=size(B1dk,2);  % number of external inputs
% no=size(C1dk,1);  % number of regulated outputs
% 
% 
% % H-infinity Dynamic Output Feedback Controller Synthesis
% 
% 
% % Declare the variables
% gamma=sdpvar(1);               % represents the bound on the H-infinity norm of the CL system.
% X1=sdpvar(ns);
% Y1=sdpvar(ns);
% An=sdpvar(ns,ns,'full')
% Cn=sdpvar(nc,ns,'full');
% Dn=sdpvar(nc,nm,'full');
% Bn=sdpvar(ns,nm,'full');
% 
% % declare constraints
% F=[X1>=eps*eye(ns)]
% F=[F;Y1>=eps*eye(ns)]
% F=[F;[X1 eye(ns); eye(ns) Y1]>=0];
% MAT=[Adk*Y1+Y1*Adk'+B2dk*Cn+Cn'*B2dk'  (Adk'+An+(B2dk*Dn*C2dk)')'        B1dk+B2dk*Dn*D21dk           (C1dk*Y1+D12dk*Cn)'; 
%      Adk'+An+(B2dk*Dn*C2dk)'         X1*Adk+Adk'*X1+Bn*C2dk+C2dk'*Bn'    X1*B1dk+Bn*D21dk           (C1dk+D12dk*Dn*C2dk)'  ;
%      (B1dk+B2dk*Dn*D21dk)'           (X1*B1dk+Bn*D21dk)'             -gamma*eye(nd)          (D11dk+D12dk*Dn*D21dk)'  ;
%      C1dk*Y1+D12dk*Cn              C1dk+D12dk*Dn*C2dk                D11dk+D12dk*Dn*D21dk         -gamma*eye(no)];
% 
% F=[F;MAT<=0];
% OPTIONS = sdpsettings('solver','sedumi')
% 
% % Solve the LMI, minimizing gamma
% optimize(F,gamma,OPTIONS)
% gamman=value(gamma)
% 
% 
% % retrieve decision variables
% X1n=value(X1); Y1n=value(Y1); Ann=value(An);Bnn=value(Bn);Cnn=value(Cn);Dnn=value(Dn);
% temp1=[Ann Bnn; Cnn Dnn]-[X1n*Adk*Y1n zeros(ns,nm); zeros(nc,ns) zeros(nc,nm)];
% 
% % Choose X2, Y2, so that X2*Y2=I-X1*Y1;
% Y2n=eye(ns);X2n=eye(ns)-X1n*Y1n;
% 
% % Reverse variable substitution
% temp2=inv([X2n X1n*B2dk;zeros(nc,ns) eye(nc)])*temp1*inv([Y2n' zeros(ns,nm); C2dk*Y1n eye(nm)]);
% Ak2=temp2(1:ns,1:ns);Bk2=temp2(1:ns,(ns+1):(ns+nm));Ck2=temp2((ns+1):(ns+nc), 1:ns);Dk2=temp2((ns+1):(ns+nc), (ns+1):(ns+nm));
% Dk=inv(eye(nc)-Dk2*D22dk)*Dk2;
% Bk=Bk2*(eye(nm)-D22dk*Dk);
% Ck=(eye(nc)-Dk*D22dk)*Ck2;
% Ak=Ak2-Bk*inv(eye(nm)-D22dk*Dk)*D22dk*Ck;
% 
% % Close the loop with Lower LFT
% plant=ss(Adk,[B1dk B2dk],[C1dk;C2dk],[D11dk D12dk; D21dk D22dk]);
% controller=ss(Ak,Bk,Ck,Dk);
% sys_cl=lft(plant,controller);
% 
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % now reformulate the upper plant and find the scalings which minimize mu
% 
% [Ath,Bth,Cth,Dth] = ssdata(sys_cl);
% 
% % Declare Scalings
% th1=sdpvar(1);
% th2=sdpvar(1);
% th3=sdpvar(1);
% Th=[th1 0 0;0 th2 0; 0 0 th3];
% %Th=diag([th1;th2;th3],0)
% %Th=blkdiag(blkdiag(th1,th2),th3)
% F=[];
% F=[F;Th>=0];
% 
% % Lyap variable
% tol=.01;
% gam_u=100;
% gam_l=0;
% ntemp=size(Ath,1);
% X=sdpvar(ntemp);
% F=[F;X>=eta*eye(ntemp)];
% %B=B1;C=C
% gam_new=gam_u;
% err=gam_u;
% while err>tol
%     MAT=[Ath'*X+X*Ath X*Bth;Bth'*X -Th]+1/gam_new/gam_new*[Cth Dth]'*Th*[Cth Dth];
%     Ftemp=[F;MAT<=0];
%     DIAG=optimize(Ftemp);
%     if DIAG.problem==0
%         gam_u=gam_new;
%     else
%         gam_l=gam_new;
%     end
%     gam_new=(gam_u+gam_l)/2
%     err=gam_u-gam_l;
% end
% 
% Thn=value(Th);
% %Thg=sqrt(blkdiag(Thn,eye(nro))); Thgr=sqrt(blkdiag(Thn,eye(1))) 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 



