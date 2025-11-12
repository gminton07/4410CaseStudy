%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Homework 5
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Author: Holt Stewart and Gabe Minton
% Class: ESE 4410
% Date: 11/05/2025

clc; clear; close all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Part 2.4 Numerical Specification
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
syms kf kfa kca Cc Cf gamma uf up 

A = [-(kf+kfa)/Cf kf/Cf 0; ... 
        kf/Cc -(kf+kca)/Cc 0; ...
            0 0 0];
B = [gamma*uf/Cf gamma*up/Cf; ...
        0 0; ...
            -1 0];

Wc = [B A*B A^2*B];

disp('Controllability Matrix')
disp(Wc)
disp('Rank(Wc)= ')
disp(rank(Wc))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Part 2.6 Numerical Specification
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Define Numerical Values
Cf = 500;       % J/C°
Cc = 2500;      % J/C°
kf = 150;       % W/°C
kfa = 20;       % W/°C
kca = 25;       % W/°C
gamma = 1000;   % W/(g/s*fan)
up0 = 0.4;      % g/s
uf0 = 0.6;      % unitless
Tamb = 25;      % °C

% Define A, B, C, and D
A = [-(kf+kfa)/Cf kf/Cf 0; ... 
        kf/Cc -(kf+kca)/Cc 0; ...
            0 0 0];
B = [gamma*uf0/Cf gamma*up0/Cf; ...
        0 0; ...
            -1 0];
C = [0 1 0];
D = 0;

% Define Steady State Temperatures for the Previously Given Parameters
Tc0 = 30; % °C
Tf0 = 31; % °C
mp0 = 0; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Part 2.7 Design Challenges
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Manually Verify Correct Matrices
disp('A Matrix:')
disp(A)
disp('B Matrix:')
disp(B)
disp('C Matrix:')
disp(C)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Part 2.7.1 Step Response and Model Validation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Define Input (u_p: Unit Step; u_f: Constant)
u = @(t) [1; 0];

f = @(t, x) A*x + B*u(t);

% Initial Condition
x0 = [0; 0; 0];

% ODE45 Solver
tspan = [0 400];
[t, x] = ode45(f, tspan, x0);

% Output y=Cx
Tc = C * x.';

% Plot Tc in the State Space
figure;
plot(t, Tc, 'LineWidth', 2);
grid on;
xlabel('Time (s)');
ylabel('State Space T_c (°C)');
title('State Space Unit Step Response');

% Nonlinear Model Solution 

% Define Inputs
up = @(t) up0 + 1*(t>=0);
uf = @(T) uf0;
q_comb = @(up,uf) gamma*up.*uf;

% Create Model Function
f_model = @(t, S) [ ...
    ( -kf*(S(1)-S(2)) - kfa*(S(1)-Tamb) + q_comb(up(t), uf(t)) )/Cf;  % Tf_dot
    (  kf*(S(1)-S(2)) - kca*(S(2)-Tamb) )/Cc;                         % Tc_dot
    - up(t) ];                                                         % mp_dot

% ODE45 Solutions
S0 = [Tf0; Tc0; mp0];
[t_nl, S] = ode45(f_model, tspan, S0);
% Tf_nl = Xnl(:,1);
Tc_nl = S(:,2);

% Plot Overlayed States
figure;
plot(t, Tc+Tc0, 'LineWidth', 2); hold on;
plot(t_nl, Tc_nl, '--', 'LineWidth', 2);
grid on; xlabel('Time (s)'); ylabel('T_c (°C)');
title('Unit Step of Linearized and Nonlinear System');

% Time Constant Calculation
Tc_final = Tc(end);
per67 = Tc_final*(1-exp(-1));
disp('Value at 67%:')
disp(per67)
index = find(Tc >= per67, 1, 'first');
disp('Time Constant (Seconds):')
disp(index)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Part 2.7.2 Integral Action
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Define New Timespan
tspan = [0 30];

% A_Tilda Matrix (4, 4)
At = [A [0; 0; 0];
      -C 0];
% B_Tilda Matrix (4, 1)
Bt = [B;
      0 0];
% C_Tilda Matrix (1, 4))
Ct = [C 0];


% Controllability Matrix
Wc = [Bt At*Bt At^2*Bt At^3*Bt];
Rank_Wc = rank(Wc);
fprintf('Rank(W_c) = %d \n', Rank_Wc)

% Desired Poles
p = [-1 -2 -3 -4];

% Pole Solver Using 'place()'
K = place(At, Bt, p);

% Verify Eigenvalues of (At-Bt*K)
eig_values = eig(At - Bt*K);
disp('Eigenvalues: ')
disp(eig_values)

% Desired Temperature
Tc_sp = 110;
r = Tc_sp - Tc0;

% Define Kr (Feed Forward)
alpha = 0;
a = 0; 
b = 1;
Kr = alpha*[a; b];

% Define Augmented State Space
f_aug = @(t, xt) (At - Bt*K)*xt + (Bt*Kr + [0; 0; 0; 1])*r;

% Initial Conditions
xt0 = [0; 0; 0; 0];

% ODE45 Solver
[t_aug, xt] = ode45(f_aug, tspan, xt0);

% Output Chamber Temperature
y = (C*xt(:,1:3).').' + Tc0;

% Plot
figure;
plot(t_aug, y, 'LineWidth', 2);
grid on;
xlabel('Time (s)');
ylabel('T_c (°C)');
title('Augmented Model with Integral Action');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Part 2.7.3 Disturbance Rejection
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Disturbance (kca)
kca_t = @(t) kca*(1 + 0.5*(t>=10 & t<20));

% Time-varying A
A_time = @(t) [ -(kf+kfa)/Cf kf/Cf 0;
                 kf/Cc -(kf + kca_t(t))/Cc 0;
                 0 0 0];

% Time-varying A TILDA
At_time = @(t) [ A_time(t) [0; 0; 0];
                 -C 0 ];

% Define Augmented State Space
f_aug = @(t, xt) (At_time(t) - Bt*K)*xt + (Bt*Kr + [0; 0; 0; 1])*r;

% ODE45 Solver
[t_aug_dist, xt] = ode45(f_aug, tspan, xt0);

% Output (absolute Tc)
y_dist = xt(:,1:3)*C.' + Tc0;

% Plot
figure;
plot(t_aug, y, 'LineWidth', 2); hold on;
plot(t_aug_dist, y_dist, 'LineWidth', 2);
grid on;
xlabel('Time (s)');
ylabel('T_c (°C)');
legend('No Disturbance', 'Disturbance')
title('Augmented Model with Integral Action and Disturbance');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Part 2.7.4 Observer Design
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Part 2.7.5 Optimizing for Changing Setpoints
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Part 2.7.6 Robustness
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Part 2.7.7 Advanced Control System
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Part 2.7.4 Observer Design (corrected & short)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
