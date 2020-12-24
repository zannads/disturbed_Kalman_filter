clear
% Remember 
% n state A
% m number of inputs 
% p number of outputs
% m_w number of state and output disturbances
% m_n number of disturbances in the input

%state eq
% A n x n
A = [2, 0.4; 0, 0.2];
% Bu n x m
Bu = [1; 7];
% Bw n x m_w
Bw = [0.0001, 0; 0.0001, 0];        
% Du m x m_n
Du = 0.005;

%output eq
% C p x n 
C = [1, 1]; 
% Dw p x m_w
Dw = [0, 0.007];
D = [0];

tau = 3;            %delay steps

%check conditions for original system
O_m = obsv(A, C);
if( (length(A) - rank(O_m)) == 0)
    disp('starting system observable');
end

Rn_m = ctrb(A, Bw);
if( (length(A) - rank(Rn_m)) == 0)
    disp('starting system reachable from disturbance');
end

R_m = ctrb(A, Bu);
if( (length(A) - rank(R_m)) == 0)
    disp('starting system reachable from input');
end

%enlarge 
n = size(A, 1);
p = size(C, 1);
m = size(Bu, 2);
[~, m_n] = size(Du);
[A_tau, Bw_tau, Bu_tau, C_tau, Dw_tau, D_tau] = enlarge(A, Bw, Bu, C, Dw, Du, D, tau);

%check conditions for enlarged system
%our_condition = Bu*Du
% O_m_tau = obsv(A_tau, C_tau);
% if( length(A_tau) - rank(O_m_tau) == 0 )
%     disp('final system observable');
% end
% R_m_tau = ctrb(A_tau, Bu_tau);
% if( length(A_tau) - rank(R_m_tau) == 0 )
%     disp('final system reachable');
% end
% 
% syms thetaclc
% 
% condi = [ A_tau-eye(size(A_tau))*exp(i*theta), Bw_tau;
%   C_tau, Dw_tau];


%compute kalman
variance_mat = [Bw_tau; Dw_tau]*[Bw_tau; Dw_tau]';
Q_tilde = Bw_tau*Bw_tau';
R_tilde = Dw_tau*Dw_tau';
%%

[Y, L, A_LC, INFO]=idare(A_tau',C_tau',Q_tilde,R_tilde);
L = -L';
%Y

Q_tilde_ = Bw*Bw' + Bu*(Du*Du')*Bu'/(tau+1);
R_tilde_ = Dw*Dw';

[Y_, L_, A_LC_, INFO_]=idare(A',C',Q_tilde_,R_tilde_);
L_ = -L_';
%Y_

% if you do it wrong but considering more variance
Q_tilde_err = Bw*Bw' + Bu*(Du*Du')*Bu';
% if you do it wrong without considering the new variance
Q_tilde_err = Bw*Bw';
R_tilde_err = Dw*Dw';
%Q_tilde_err( size(A_tau,1), size(A_tau,1) ) = 0;
%R_tilde_err( size(C_tau,1), size(C_tau,1) ) = 0;

[Y_err, L_err, A_LC_err, INFO_err]=idare(A',C',Q_tilde_err,R_tilde_err);
L_err = -L_err';

%% partial pole placement 
T = [A^tau, A^2*Bu, A*Bu, Bu];
F = place (A, Bu, [0.5, 0.2]);
K = F*T;

eig(A_tau-Bu_tau*K)
%% simulink preparation
% the real system evolves without -Bu*Du ...
x0 = [1;0];

B_simu = [Bu, Bw];
D_simu = [zeros(p,m), Dw];
C_simu_state = eye( n, n);
m_en =size(B_simu, 2);
D_simu_state = zeros(n, m_en);

 A_KF = A_tau + L*C_tau;
 B_KF = [Bu_tau, -L];
 [n_kf, ~] = size(A_KF);
 [~, m_kf] = size(B_KF);
 C_KF = eye( n_kf, n_kf );
 D_KF = zeros(n_kf, m_kf);
 
[n_kf, ~] = size(A);
A_KF_red = A + L_*C;
B_KF_red = [Bu, -L_];
[~, m_kf] = size(B_KF_red);
C_KF_red = eye( n_kf, n_kf);
D_KF_red = zeros(n_kf, m_kf);

%non funziona perchè questa formulazione è con zeri nel cerchio unitario
A_KF_err = A + L_err*C;
B_KF_err = [Bu, -L_err];
[n_kf, ~] = size(A_KF_err);
[~, m_kf] = size(B_KF_err);
C_KF_err = eye( n_kf, n_kf);
D_KF_err = zeros(n_kf, m_kf);


sim('disturbed_kalman');
%% functions
function [A_tau, Bw_tau, Bu_tau, C_tau, Dw_tau, D_tau] = enlarge(A, Bw, Bu, C, Dw, Du, D, tau)
%get dimensions 
[n, m, p, m_w ] = check_dimensions(A, Bw, Bu, C, Dw, D);
m_n = size(Du, 2);
if( n ~=0 ) %if everything is ok
    
    A_tau = [A,     Bu,                 zeros(n, m*(tau-1) );
            zeros(m*(tau-1), n+m),      eye(m*(tau-1));
            zeros(1, n+m*tau)];
    Bu_tau = [zeros(n+m*(tau-1), m);
              eye(m)];
          
    Bw_tau = zeros(n + tau*m, m_w + tau*m +m_n);
    Bw_tau( 1: size(Bw,1), 1:size(Bw,2) ) = Bw;
    Bw_tau( (end-size(Du,1)+1):end, (end-size(Du,2)+1):end) = Du;
    
    C_tau = [C,                 zeros(p, tau*m);  
             zeros(tau*m,n),    eye(tau*m)];
    
    %Dw_tau = zeros(p+tau*m, m_w + tau*m +m_n);
    %Dw_tau( 1: size(Dw,1), 1:size(Dw,2) ) = Dw;
    Dw_tau = Dw;
    for idx = 1: tau
        Dw_tau = blkdiag(Dw_tau, Du);
    end
    Dw_tau = [Dw_tau, zeros(p+tau*m, m_n)];
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