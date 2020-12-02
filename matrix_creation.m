% Remember 
% n state A
% m number of inputs 
% p number of outputs
% m_w number of state and output disturbances
% m_n number of disturbances in the input

%state eq
% A n x n
A = [0.3, 0.4; 0, 0.5];
% Bu n x m
Bu = [1; 7];
% Bw n x m_w
Bw = [1, 2.3; 4.1, 1.7];        
% Du m x m_n
Du = 0.5 ;

%output eq
% C p x n 
C = [1, 1]; 
% Dw p x m_w
Dw = [0.4, 0.3];
D = 0;

tau = 5;            %delay steps

%check conditions for original system


%enlarge 
[~, m_n] = size(Du);
[p, ~] = size(C);
[A_tau, Bw_tau, Bu_tau, C_tau, Dw_tau, D_tau] = enlarge(A, [-Bu*Du,  Bw], Bu, C, [zeros(p, m_n), Dw], D, tau);

%check conditions for enlarged system
Bu*Du

%compute kalman
Q_tilde = Bw_tau*Bw_tau';
R_tilde = Dw_tau*Dw_tau';

[Y,K,L]=idare(A_tau',C_tau',Q_tilde,R_tilde);
K = -K'
Y
L
    
%% simulink preparation
% the real system evolves without -Bu*Du ...
[A_simu, Bw_simu, Bu_simu, C_simu, Dw_simu, Du_simu] = enlarge(A, Bw, Bu, C, Dw, D, tau);
B_simu = [Bu_simu, Bw_simu];
D_simu = [Du_simu, Dw_simu];
[n_en, ~] = size(A_simu);
[~, m_en] =size(B_simu);
C_simu_state = eye( n_en, n_en);
D_simu_state = zeros(n_en, m_en);

A_KF = A_tau + K*C_tau;
B_KF = [Bu_tau, -K];
[n_kf, ~] = size(A_KF);
[~, m_kf] =size(B_KF);
C_KF = eye( n_kf, n_kf );
D_KF = zeros(n_kf, m_kf);
sim('disturbed_kalman');
%% functions
function [A_tau, Bw_tau, Bu_tau, C_tau, Dw_tau, D_tau] = enlarge(A, Bw, Bu, C, Dw, D, tau)
%get dimensions 
[n, m, p, m_w ] = check_dimensions(A, Bw, Bu, C, Dw, D);

if( n ~=0 ) %if everything is ok
    
    A_tau = [A, Bu, zeros(n, tau-1);
        zeros(tau-1, n+m) eye(tau-1);
        zeros(1, n+tau)];
    Bu_tau = [zeros(n+tau-1, m);
        eye(m)];
    Bw_tau = [ Bw;
        zeros(tau, m_w)];
    
    C_tau = [C, zeros(p, tau)];
    Dw_tau = Dw;
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