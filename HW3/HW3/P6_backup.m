%% Part (b)
disp("#####Problem 6 Part B Output:######")
disp("(i) LMI to formulate H_inf optimal state-feedback problem")
% LMI for state feedback (Lecture 9 Theorem 5)
% 
% % Define variables
% Z = sdpvar(3,6);    % Define Z
% Y = sdpvar(6,6);    % Define Y
% gamma = sdpvar(1,1);
% 
% % Define LMIs and Constraints
% optimize_var = gamma;
% LMI1 = [Y>=eta*eye(6)];  
% LMI2 = [[
% (Y*A' + A*Y + Z'*B2'+ B2*Z)    (B1)    (Y*C1'+Z'*D12');
% (B1')   (-gamma*eye(6))  (D11');
% (C1*Y+D12*Z)    (D11)   (-gamma*eye(6))]<=0 ];
% Fun = [LMI1;LMI2];              % Combine into single constraint function
% 
% % Optimize
% optimize(Fun,optimize_var,ops);          % Run the optimization
% 
% % Feasable Results
% Zf = value(Z);                  % Return feasable solution Z
% Yf = value(Y);                  % Return feasable solution Y
% gammaf = value(gamma);          % Return feasable solution gamma
% 
% % Construct solutions
% F=Zf*inv(Yf);             % Construct F
% 
% disp("F=")
% disp(F)
% disp("(ii) What is the closed-loop H_inf gain?")
% H_inf_norm = gammaf 

%% Part (e)*
disp("#####Problem 6 Part E Output:######")
disp("(i) Use an LMI to determine the Hinf gain of the closed loop system")
disp("INCOMPLETE")

% X = sdpvar(6,6);
% LMI1 = [X>=eta*eye(6)]
% LMI2_1 = [
% (A'*X+X*A)  (X*B);
% (B'*X)  (-gamma*eye(3))
% ];
% LMI2_2 = [C'; D'];
% LMI2_3 = [C D];
% LMI2 = [LMI2_1 + 1/gamma*(LMI2_2)*LMI2_3 <= 0];
% 
% Fun = [LMI1; LMI2];
% 
% % Optimize
% optimize(Fun,empty,ops);          % Run the optimization
% 
% % Feasable Results
% Xf = value(X);                  % Return feasable solution X
% gammaf = value(gamma);          % Return feasable solution gamma
% 
% disp("H_inf gain of closed loop system = ")
% disp(gammaf)

% gamma = sdpvar(1,1);
% Y=sdpvar(6);
% Z = sdpvar(3,6);
% 
% optimize_var = gamma;
% LMI1 = [Y];
% LMI2 = [
%     Y*A'+A*Y+Z'*B2'+B2*Z B1 Y*C1'+Z'*D12';
%     B1' -gamma*eye(6) D11';
%     C1*Y+D12*Z D11 -gamma*eye(6);
%     ];
% 
% Fun = [
%     LMI1>=eta*eye(6);
%     LMI2<=-eta*eye(18);
%     ];
% 
% optimize(Fun,optimize_var,ops);          % Run the optimization
% 
% gammaf = value(gamma);          % Return feasable solution gamma
% Yf = value(Y);                  % Return feasable solution Y
% Zf = value(Z);                  % Return feasable solution Z
% F = Zf*inv(Yf);                   % Construct F
% 
% 
% disp("F=")
% disp(F)
% % (ii) Determine the closed-loop H2 gain
% disp("H_inf gain = ")
% disp(gammaf)
% 

