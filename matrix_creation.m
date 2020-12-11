% Remember 
% n state A
% m number of inputs 
% p number of outputs
% m_w number of state and output disturbances
% m_n number of disturbances in the input

%state eq
% A n x n
A = [0.3, 0.4; 0, 0.2];
% Bu n x m
Bu = [1; 7];
% Bw n x m_w
Bw = [0.001, 0; 0.002, 0];        
% Du m x m_n
Du = 0.3;

%output eq
% C p x n 
C = [1, 1]; 
% Dw p x m_w
Dw = [0, 0.4];
D = 0;

tau = 3;            %delay steps

%check conditions for original system
O_m = obsv(A, C);
if( (length(A) - rank(O_m)) == 0)
    disp('starting system observable');
end
R_m = ctrb(A, Bw);
if( (length(A) - rank(R_m)) == 0)
    disp('starting system reachable');
end

% 
% 
% condi = [ A-eye(size(A))*exp(i*theta), Bw;
%   C, Dw]
% test = rank(condi);

%enlarge 
[~, m_n] = size(Du);
[p, ~] = size(C);
[A_tau, Bw_tau, Bu_tau, C_tau, Dw_tau, D_tau] = enlarge(A, Bw, Bu, C, Dw, Du, D, tau);

%check conditions for enlarged system
%our_condition = Bu*Du
O_m_tau = obsv(A_tau, C_tau);
if( length(A_tau) - rank(O_m_tau) == 0 )
    disp('final system observable');
end
R_m_tau = ctrb(A_tau, Bw_tau);
if( length(A_tau) - rank(R_m_tau) == 0 )
    disp('final system reachable');
end

syms theta
condi = [ A_tau-eye(size(A_tau))*exp(i*theta), Bw_tau;
  C_tau, Dw_tau];


%compute kalman
variance_mat = [Bw_tau; Dw_tau]*[Bw_tau; Dw_tau]';
Q_tilde = Bw_tau(:, 2:3)*Bw_tau(:, 2:3)';
R_tilde = Dw_tau(:, 2:3)*Dw_tau(:, 2:3)';
%%
[Y,K,L]=idare(A_tau',C_tau',Q_tilde,R_tilde);
K = -K'
Y
L


    
%% simulink preparation
% the real system evolves without -Bu*Du ...
n = size(A, 1);
p = size(C, 1);
m = size(Bu, 2);
B_simu = [Bu, Bw];
D_simu = [zeros(p,m), Dw];
C_simu_state = eye( n, n);
m_en =size(B_simu, 2);
D_simu_state = zeros(n, m_en);
% 
% A_KF = A_tau + K*C_tau;
% B_KF = [Bu_tau, -K];
% [n_kf, ~] = size(A_KF);
% [~, m_kf] = size(B_KF);
% C_KF = eye( n_kf, n_kf );
% D_KF = zeros(n_kf, m_kf);
% 
% [n_kf, ~] = size(A);
% A_KF_red = A + K(1:n_kf)*C;
% B_KF_red = [Bu, -K(1:n_kf)];
% [~, m_kf] = size(B_KF_red);
% C_KF_red = eye( n_kf, n_kf);
% D_KF_red = zeros(n_kf, m_kf);


sim('disturbed_kalman');
%% functions
function [A_tau, Bw_tau, Bu_tau, C_tau, Dw_tau, D_tau] = enlarge(A, Bw, Bu, C, Dw, Du, D, tau)
%get dimensions 
[n, m, p, m_w ] = check_dimensions(A, Bw, Bu, C, Dw, D);
m_n = size(Du, 2);
if( n ~=0 ) %if everything is ok
    
    A_tau = [A, Bu, zeros(n, tau-1);
        zeros(tau-1, n+m) eye(tau-1);
        zeros(1, n+tau)];
    Bu_tau = [zeros(n+tau-1, m);
        eye(m)];
    Bw_tau = [zeros(n+(tau-1)*m, m_n), [Bw; zeros((tau-1)*m, m_w)];
              -Du, zeros(m, m_w)];
    
    C_tau = [C, zeros(p, tau*m);  
             zeros(tau*m,n), eye(tau*m)];
    
    Dw_tau = [zeros(p, m_n),        Dw;
              ones(tau, m_n)*Du,    zeros(tau, m_w)];
    D_tau = D;
end
end

function [n, m, p, m_w ] = check_dimensions(A, Bw, Bu, C, Dw, D)
% return the dimension of the system if wverything is ok,n = 0 when not
% n state dimension
% m input number
% p output number
% m_w disturbances number

[A_r,  A_c]  = size(A);
[Bw_r, Bw_c] = size(Bw);
[Bu_r, Bu_c] = size(Bu);
[C_r,  C_c]  = size(C);
[Dw_r, Dw_c] = size(Dw);
[D_r, D_c] = size(D);

%states are all the same size
c1 = (A_r == A_c && A_c == Bw_r && Bw_r == Bu_r && Bu_r == C_c); % states = n;
c2 = (Bw_c == Dw_c);    % disturbances = m_w;
c3 = (C_r == Dw_r && C_r == D_r);     % outputs = p;
c4 = (Bu_c == D_c);             % inputs = m

if( c1 && c2 && c3 && c4)
    n = A_r;
    m = Bu_c;
    p = C_r;
    m_w = Bw_c;
else
    n = 0;
end

end