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

C = [0 1 0];

Wc = [B A*B A^2*B];

disp('Controllability Matrix')
disp(Wc)
disp('Rank(Wc)= ')
disp(rank(Wc))

Wo = [C; C*A; C*A^2];

disp('Observability Matrix')
disp(Wo)
disp('Rank(Wo)= ')
disp(rank(Wo))


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

% Define b and q(u)
b = [kfa*Tamb/Cf; kca*Tamb/Cc; 0];
q = [gamma*up0*uf0/Cf; 0; -up0];

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

% Define A, B, C, and D Without Mass Since it Does Not Effect Dynamics
A = [-(kf+kfa)/Cf kf/Cf; ... 
        kf/Cc -(kf+kca)/Cc];
B = [gamma*uf0/Cf gamma*up0/Cf; ...
        0 0];
C = [0 1];
D = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Part 2.7.1 Step Response and Model Validation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Define Input (u_p: Unit Step; u_f: Constant)
u = @(t) [0+1*(t>0); 0];

f = @(t, x) A*x + B*u(t);

% Initial Condition
x0 = [0; 0];

% ODE45 Solver
tspan = [0 400];
[t, x] = ode45(f, tspan, x0);

% Output y=Cx
Tc = C * x.';

% Nonlinear Model Solution 
% Define Inputs
up = @(t) up0 + 1*(t>=0);
uf = @(T) uf0;
q = @(up,uf) [gamma*up.*uf./Cf; 0];
b = [kfa*Tamb/Cf; kca*Tamb/Cc];


% Define the Nonlinear System
f_model = @(t, S) A*S + b + q(up(t),uf(t));
          
% ODE45 Solutions
S0 = [Tf0; Tc0];
[t_nl, S] = ode45(f_model, tspan, S0);
Tc_nl = S(:,2);

% Plot Nonlinear and Linearized System
figure;
plot(t, Tc+Tc0, 'b', 'LineWidth', 2); hold on;
plot(t_nl, Tc_nl, 'r--', 'LineWidth', 2);
grid on; xlabel('Time (s)'); ylabel('T_c (°C)');
legend('Ax+Bu', 'Ax+b+q(u)');
title('Unit Step of Linearized and Nonlinear System');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Part 2.7.2 Integral Action Linear System
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Define New Timespan
tspan = [0 200];

% A_Tilda Matrix (4, 4)
At = [A [0; 0];
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
p = [-0.1 -0.2 -0.3];

% Pole Solver Using 'place()'
K = place(At, Bt, p);
Kx = K(:,1:2); 
KI = K(:,3); 

% Verify Eigenvalues of (At-Bt*K)
eig_values = eig(At - Bt*K);
disp('Eigenvalues: ')
disp(eig_values)

% Desired Temperature
Tc_sp = 110;
r = Tc_sp - Tc0;

% Define Kr (Feed Forward)
alpha = 0;
a = 0.01; 
b = 0.1;
Kr = alpha*[a; b];

% Define Linear System
f_aug = @(t, xt) (At - Bt*K)*xt + (Bt*Kr + [0; 0; 1])*r;


% Initial Conditions
xt0 = [0; 0; 0];

% ODE45 Solver
[t_aug, xt] = ode45(f_aug, tspan, xt0);

% Output Chamber Temperature
Tc_lin = xt(:,2) + Tc0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Part 2.7.2 Integral Action Nonlinear System
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Nonlinear Inputs and Losses to Ambient Temp
q_fun = @(up, uf) [gamma*up.*uf./Cf; 0];
b_vec = [kfa*Tamb/Cf; kca*Tamb/Cc];

% Construct Differences
xdev   = @(X) [X(1) - Tf0; ...
               X(2) - Tc0];

% Augmented Matrix
xtilde = @(X) [xdev(X); X(3)];

% Define Inputs
uf = @(X) -K(1,:)*xtilde(X) + Kr(1)*Tc_sp;
up = @(X) -K(2,:)*xtilde(X) + Kr(2)*Tc_sp;

% Nonlinear Closed Loop System
f_cl = @(t, X) [ ...
      A * [X(1); X(2)] + ...                        % Linear Component
      b_vec + ...                                   % Ambient Losses
      q_fun(up0 + up(X), uf0 + uf(X)); ...          % Nonlinear Inputs
      Tc_sp - X(2)                                  % Integral Action
];

% ODE45 Solver For Nonlinear System
X0 = [Tf0; Tc0; 0];
[t_nl, X] = ode45(f_cl, tspan, X0);

% Tc From Nonlinear System
Tc_nl = X(:,2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Part 2.7.2 Integral Action Plots
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure;
plot(t_aug, Tc_lin, 'b', 'LineWidth', 2); hold on;
plot(t_nl, Tc_nl, 'r--', 'LineWidth', 2);
grid on;
xlabel('Time (s)');
ylabel('T_c (°C)');
legend('Linear Model with Integral Action', ...
       'Nonlinear Model with Integral Action');
title('Integral Action: Linear vs Nonlinear Model');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Part 2.7.3 Disturbance Rejection Linear System
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Disturbance (kca)
kca_t = @(t) kca*(1 + 0.5*(t>=10 & t<20));

% Time-varying A
A_time = @(t) [ -(kf+kfa)/Cf kf/Cf;
                 kf/Cc -(kf + kca_t(t))/Cc];

% Time-varying A TILDA
At_time = @(t) [ A_time(t) [0; 0];
                 -C 0 ];

% Define Augmented State Space
f_aug = @(t, xt) (At_time(t) - Bt*K)*xt + (Bt*Kr + [0; 0; 1])*r;

% ODE45 Solver
[t_aug_dist, xt] = ode45(f_aug, tspan, xt0);

% Output (absolute Tc)
Tc_dist_lin = xt(:,1:3)*Ct.' + Tc0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Part 2.7.3 Disturbance Rejection Nonlinear System
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Nonlinear Inputs
q_fun = @(up, uf) [gamma*up.*uf./Cf; 0];

% Time Varient Ambient Loss
b_vec = @(t) [kfa*Tamb/Cf; ...
              kca_t(t)*Tamb/Cc];

% Deviation From Steady State
xdev = @(X) [X(1) - Tf0; ...
             X(2) - Tc0];

% Augmented Deviation State
xtilde = @(X) [xdev(X); X(3)];

% Define Inputs
uf = @(X) -K(1,:)*xtilde(X) + Kr(1)*Tc_sp;
up = @(X) -K(2,:)*xtilde(X) + Kr(2)*Tc_sp;

% Nonlinear closed-loop system with time-varying kca(t)
% X = [Tf; Tc; xI]
f_cl = @(t, X) [ ...
      A_time(t) * [X(1); X(2)] + ...              % Linear Component
      b_vec(t) + ...                              % Ambient Terms
      q_fun( up0 + up(X), uf0 + uf(X) ); ...      % Nonlinear Input
      Tc_sp - X(2)                                % Integral Action
];

% ODE45 Solver for Nonlinear Disturbed System
X0 = [Tf0; Tc0; 0];
[t_dist_nl, X] = ode45(f_cl, tspan, X0);

% Chamber temperature from nonlinear system
Tc_dist_nl = X(:,2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Part 2.7.3 Disturbance Rejection Plots
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot of Disturbance and Nondisturbance in Linearized System
figure;
plot(t_aug, Tc_lin, 'b', 'LineWidth', 2); hold on;
plot(t_aug_dist, Tc_dist_lin, 'r--', 'LineWidth', 2);
grid on;
xlabel('Time (s)');
ylabel('T_c (°C)');
legend('No Disturbance', 'Disturbance')
title('Augmented Model with Integral Action and Disturbance');

% Disturbance for Linear and Nonlinear Systems
figure;
plot(t_aug_dist, Tc_dist_lin, 'b', 'LineWidth', 2); hold on;
plot(t_dist_nl, Tc_dist_nl, 'r--', 'LineWidth', 2);
grid on;
xlabel('Time (s)');
ylabel('T_c (°C)');
legend('Linear Model with Integral Action', ...
       'Nonlinear Model with Integral Action');
title('Integral Action: Linear vs Nonlinear Model With Disturbance');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Part 2.7.4 Observer Design Linear System
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Design Observer for Negative Eigenvalues
p_obs = [-2 -1];
L = place(A', C', p_obs).';

% Combine the Observer and the Integral Action
f_comb_lin = @(t, z) [ ...
    % xdot = A x + B u
    A*z(1:2) + B*(-Kx*z(4:5) - KI*z(3)); ...
    % xIdot = Tcsp - y = Tcsp - Cx
    r - C*z(1:2); ...
    % xhat_dot = A xhat + B u + L(y - C xhat)
    A*z(4:5) + B*(-Kx*z(4:5) - KI*z(3)) + L*(C*z(1:2) - C*z(4:5))
];

% Define Initial Conditions
z0 = [0; 0; 0; 0; 0];

% ODE45 Solver
[t_lin_obs, z_lin] = ode45(f_comb_lin, tspan, z0);

% View States
x_actual = z_lin(:,1:2);
xI     = z_lin(:,3);
x_hat  = z_lin(:,4:5);

% Chamber temperatures (absolute)
Tc_actual_lin = x_actual(:,2) + Tc0;
Tc_hat_lin  = x_hat(:,2)  + Tc0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Part 2.7.4 Observer Design Nonlinear System
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Nonlinear Input
q_fun = @(up, uf) [gamma * up .* uf ./ Cf; 0];

% Nonlinear Matrix
A_nl = [-(kf + kfa)/Cf, kf/Cf; ...
         kf/Cc, -(kf + kca)/Cc];

% Ambient Losses
b_nl = [kfa * Tamb / Cf; ...
        kca * Tamb / Cc];

% Initial Conditions
z0_nl = [Tf0; Tc0; 0; 0; 0];

% ODE45 Solver Using the Function 
[t_nl_obs, z_nl] = ode45( ...
        @(t,z) f_nl_obs(t, z, A_nl, b_nl, A, B, C, Kx, KI, L, Kr, q_fun, uf0, up0, Tc0, r), ...
        tspan, z0_nl);

% Results
Tf_nl   = z_nl(:,1);
Tc_nl   = z_nl(:,2);
xI_nl   = z_nl(:,3);
xhat_nl = z_nl(:,4:5);

Tf_hat = xhat_nl(:,1) + Tf0;
Tc_hat = xhat_nl(:,2) + Tc0;
e_Tf = Tf_nl - Tf_hat;
e_Tc = Tc_nl - Tc_hat;


function dz = f_nl_obs(t, z, A_nl, b_nl, A, B, C, Kx, KI, L, Kr, q_fun, uf0, up0, Tc0, r)
% Closed Loop Observer Function
% z = [ Tf; Tc; xI; Tf_hat_dev; Tc_hat_dev ]

% Define States
Tf = z(1);
Tc = z(2);
xI = z(3);
xhat = z(4:5);

x_nl = [Tf; Tc];
y = Tc - Tc0;

% Define the Input
u = -Kx * xhat - KI * xI + Kr * r;

% Individual Inputs
uf = uf0 + u(1);
up = up0 + u(2);

% Define Nonlinear Dynamics
xdot = A_nl * x_nl + b_nl + q_fun(up, uf);

% Define Integral Action
xIdot = r - y;

% Define the Model System
xhat_dot = A * xhat + B * u + L * (y - C * xhat);

% Augmented State Space
dz = [xdot;
      xIdot;
      xhat_dot];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Part 2.7.4 Observer Design Mass Extrapolation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Extract Number of Points
N = length(t_nl_obs);
% Define U
u_hist = zeros(N, 2);

for k = 1:N
    % xI at Time k
    xI_k   = xI_nl(k);
    % xhat at Time k
    xhat_k = xhat_nl(k, :).';
    % u at Time k
    u_k = -Kx * xhat_k - KI * xI_k + Kr * r;
    % Store Inputs at Time k
    u_hist(k, :) = u_k.';
end

uf_k = u_hist(:,1);
up_k = u_hist(:,2);

% Pellet Input and Fan Input
up_actual = up0 + up_k;
uf_actual = uf0 + uf_k;
% Time Steps
dt = t_nl_obs(2:end) - t_nl_obs(1:end-1);
% Total Mass
mp = mp0 - cumsum(up_actual(1:end-1) .* dt);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Part 2.7.4 Observer Design Plots
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Plot Combined Control System on the Linear System
figure;
plot(t_lin_obs, Tc_actual_lin, 'b', 'LineWidth', 2); hold on;
plot(t_lin_obs, Tc_hat_lin, 'r--', 'LineWidth', 2);
grid on;
xlabel('Time (s)');
ylabel('T_c (°C)');
legend('T_c (Actual)', 'T_c (Estimate)');
title('Luenberger Observer and Integral Action on Linear System');

% Plot Control System with the Nonlinear System
figure;
plot(t_nl_obs, Tc_nl, 'b', 'LineWidth', 2);  hold on;
plot(t_nl_obs, Tc_hat, 'r--', 'LineWidth', 2);
grid on;
xlabel('Time (sec)');
ylabel('Temperature (°C)');
title('Actual vs Observed T_c');
legend('Actual T_c', 'Observed T_c');

% Plot the Actual vs Observed Tf with Distrubances
figure;
plot(t_nl_obs, Tf_nl, 'b', 'LineWidth', 2);  hold on;
plot(t_nl_obs, Tf_hat, 'r--', 'LineWidth', 2);
grid on;
xlabel('Time (s)');
ylabel('Temperature (°C)');
title('Actual vs Observed T_f');
legend('Actual T_f', 'Observed T_f');

% Plot the Observer Error for Tc and Tf
figure;
plot(t_nl_obs, e_Tf, 'b', 'LineWidth', 1.5); hold on;
plot(t_nl_obs, e_Tc, 'r', 'LineWidth', 1.5);
xlabel('Time (s)');
ylabel('Estimation Error (°C)');
legend('T_f error', 'T_c error');
grid on;
title('Observer Estimation Error');

% Pellet Mass
figure; 
plot(t_nl_obs(1:end-1), mp, 'b', 'LineWidth', 2);
grid on;
xlabel('Time (s)');
ylabel('m_p (g)');
title('Mass of Pellets')

% Plot Inputs
figure;
plot(t_nl_obs, uf_actual, 'b', 'LineWidth', 2); hold on;
plot(t_nl_obs, up_actual, 'r--', 'LineWidth', 2);
grid on;
xlabel('Time (s)');
ylabel('Control Inputs');
legend('u_f','u_p');
title('Control Inputs u_f(t) and u_p(t)');



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Part 2.7.4 Observer Design With Disturbances
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Disturbance
kca_t = @(t) kca*(1 + 6.5*(t>=10 & t<20));

% Define Time Varying A Matrix
A_time = @(t) [ -(kf+kfa)/Cf kf/Cf;
                 kf/Cc -(kf + kca_t(t))/Cc];

% Define Time Varying Losses to Ambient Temp
b_vec = @(t) [kfa*Tamb/Cf; ...
              kca_t(t)*Tamb/Cc];

% Initial Conditions
z0_nl = [Tf0; Tc0; 0; 0; 0];

% Solve nonlinear + observer closed-loop system (no saturation for now)
[t_nl_obs, z_nl] = ode45(@(t,z) ...
        f_nl_obs_dist(t, z, A_time, b_vec, A, B, C, Kx, KI, L, Kr, q_fun, uf0, up0, Tc0, r), ...
        tspan, z0_nl);

% Results
Tf_nl   = z_nl(:,1);
Tc_nl   = z_nl(:,2);
xI_nl   = z_nl(:,3);
xhat_nl = z_nl(:,4:5);

Tf_hat = xhat_nl(:,1) + Tf0;
Tc_hat = xhat_nl(:,2) + Tc0;

e_Tf = Tf_nl - Tf_hat;
e_Tc = Tc_nl - Tc_hat;

function dz = f_nl_obs_dist(t, z, A_time, b_vec, A, B, C, Kx, KI, L, Kr, q_fun, uf0, up0, Tc0, r)
% Closed Loop Observer Funciton With Disturbane
% z = [ Tf; Tc; xI; Tf_hat_dev; Tc_hat_dev ]

% Define States
Tf   = z(1);
Tc   = z(2);
xI   = z(3);
xhat = z(4:5);
x_nl = [Tf; Tc];
y = Tc - Tc0;

% Define the Inputs
u = -Kx * xhat - KI * xI + Kr * r;

% Individual Inputs
uf = uf0 + u(1);
up = up0 + u(2);

% Nonlinear System
xdot = A_time(t) * x_nl + b_vec(t) +  q_fun(up, uf);

% Integral Action
xIdot = r - y;

% Observer System
xhat_dot = A * xhat + B * u + L * (y - C * xhat);

% Augmented State Space
dz = [xdot;
      xIdot;
      xhat_dot];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Part 2.7.4 Observer Design Mass Extrapolation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Extract Number of Points
N = length(t_nl_obs);
% Define U
u_hist = zeros(N, 2);

for k = 1:N
    % xI at Time k
    xI_k   = xI_nl(k);
    % xhat at Time k
    xhat_k = xhat_nl(k, :).';
    % u at Time k
    u_k = -Kx * xhat_k - KI * xI_k + Kr * r;
    % Store Inputs at Time k
    u_hist(k, :) = u_k.';
end

uf_k = u_hist(:,1);
up_k = u_hist(:,2);

% Pellet Input and Fan Input
up_actual = up0 + up_k;
uf_actual = uf0 + uf_k;
% Time Steps
dt = t_nl_obs(2:end) - t_nl_obs(1:end-1);
% Total Mass
mp = mp0 - cumsum(up_actual(1:end-1) .* dt);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Part 2.7.4 Observer Design Plots Disturbance
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Plot the Actual vs Observed Tc with Distrubances
figure;
plot(t_nl_obs, Tc_nl, 'b', 'LineWidth', 2);  hold on;
plot(t_nl_obs, Tc_hat, 'r--', 'LineWidth', 2);
grid on;
xlabel('Time (s)');
ylabel('Temperature (°C)');
title('Actual vs Observed T_c');
legend('Actual T_c', 'Observed T_c');

% Plot the Actual vs Observed Tf with Distrubances
figure;
plot(t_nl_obs, Tf_nl, 'b', 'LineWidth', 2);  hold on;
plot(t_nl_obs, Tf_hat, 'r--', 'LineWidth', 2);
grid on;
xlabel('Time (s)');
ylabel('Temperature (°C)');
title('Actual vs Observed T_f');
legend('Actual T_f', 'Observed T_f');

% Plot the Error 
figure;
plot(t_nl_obs, e_Tf, 'b','LineWidth', 1.5); hold on;
plot(t_nl_obs, e_Tc, 'r', 'LineWidth', 1.5);
xlabel('Time (s)');
ylabel('Estimation Error (°C)');
legend('Tf error', 'Tc error');
grid on;
title('Observer Estimation Error');

% Pellet Mass
figure; 
plot(t_nl_obs(1:end-1), mp, 'b', 'LineWidth', 2);
grid on;
xlabel('Time (s)');
ylabel('m_p (g)');
title('Mass of Pellets')

% Plot Inputs
figure;
plot(t_nl_obs, uf_actual, 'b', 'LineWidth', 2); hold on;
plot(t_nl_obs, up_actual, 'r--', 'LineWidth', 2);
grid on;
xlabel('Time (s)');
ylabel('Control Inputs');
legend('u_f','u_p');
title('Control Inputs u_f(t) and u_p(t)');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Part 2.7.5 Changing Setpoints 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Cooking Profile
Tc_sp_fun = @(t) ...
    (t >=   0 & t <  600).*90  + ...   % 0–600 s:  90 °C
    (t >= 600 & t < 1200).*110 + ...   % 600–1200 s: 110 °C
    (t >=1200 & t <=1800).*130;        % 1200–1800 s:130 °C

% Setpoint
r_fun = @(t) Tc_sp_fun(t) - Tc0;

% Define New TSPAN
tspan_prof = [0 1800];

% Initial Conditions
z0_prof = [Tf0; Tc0; 0; 0; 0];

% ODE45 Solver
[t_prof, z_prof] = ode45(@(t,z) ...
        f_profile(t, z, A_nl, b_nl, A, B, C, Kx, KI, L, Kr, q_fun, uf0, up0, Tc0, r_fun), ...
        tspan_prof, z0_prof);

% Results
Tf_prof   = z_prof(:,1);
Tc_prof   = z_prof(:,2);
xI_prof   = z_prof(:,3);
xhat_prof = z_prof(:,4:5);

Tf_hat = xhat_prof(:,1) + Tf0;
Tc_hat = xhat_prof(:,2) + Tc0;

e_Tf = Tf_prof - Tf_hat;
e_Tc = Tc_prof - Tc_hat;


function dz = f_profile(t, z, A_nl, b_nl, A, B, C, Kx, KI, L, Kr, q_fun, uf0, up0, Tc0, r_fun)
% Closed Loop System with Observer and Varying Set Points
% z = [ Tf; Tc; xI; Tf_hat_dev; Tc_hat_dev ]

% Define States
Tf   = z(1);
Tc   = z(2);
xI   = z(3);
xhat = z(4:5);
x_nl  = [Tf; Tc];
y = Tc - Tc0;

% Time Varying Setpoint
r_t = r_fun(t);

% Define the Input
u = -Kx * xhat - KI * xI + Kr * r_t;

% Individual Inputs
uf = uf0 + u(1);
up = up0 + u(2);

% Nonlinear System
xdot = A_nl * x_nl + b_nl + q_fun(up, uf);

% Integral Action
xIdot = r_t - y;

% Observer System
xhat_dot = A * xhat + B * u + L * (y - C * xhat);

% Augmented System
dz = [xdot;
      xIdot;
      xhat_dot];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Part 2.7.5 Changing Setpoints Mass Extrapolation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Extract Number of Points
N = length(t_prof);
% Define U
u_hist = zeros(N, 2);

for k = 1:N
    % xI at Time k
    xI_k   = xI_prof(k);
    % xhat at Time k
    xhat_k = xhat_prof(k, :).';
    % u at Time k
    u_k = -Kx * xhat_k - KI * xI_k + Kr * r;
    % Store Inputs at Time k
    u_hist(k, :) = u_k.';
end

uf_k = u_hist(:,1);
up_k = u_hist(:,2);

% Pellet Input and Fan Input
up_actual = up0 + up_k;
uf_actual = uf0 + uf_k;
% Time Steps
dt = t_prof(2:end) - t_prof(1:end-1);
% Total Mass
mp = mp0 - cumsum(up_actual(1:end-1) .* dt);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Part 2.7.5 Changing Setpoints Plots
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Plot Tc with The Profile
figure;
plot(t_prof, Tc_prof, 'b', 'LineWidth', 2);  hold on;
plot(t_prof, Tc_hat,'r--', 'LineWidth', 1.5);
grid on;
xlabel('Time (s)');
ylabel('T_c (°C)');
legend('Actual T_c','Estimated T_c');
title('Cooking Profile Tracking (No Disturbance)');

% Plot the Actual vs Observed Tf with Distrubances
figure;
plot(t_prof, Tf_prof, 'b', 'LineWidth', 2);  hold on;
plot(t_prof, Tf_hat, 'r--', 'LineWidth', 2);
grid on;
xlabel('Time (s)');
ylabel('Temperature (°C)');
title('Actual vs Observed T_f');
legend('Actual T_f', 'Observed T_f');

% Plot the Error 
figure;
plot(t_prof, e_Tf, 'b', 'LineWidth', 1.5); hold on;
plot(t_prof, e_Tc, 'r', 'LineWidth', 1.5);
xlabel('Time (s)');
ylabel('Estimation Error (°C)');
legend('Tf error', 'Tc error');
grid on;
title('Observer Estimation Error');

% Pellet Mass
figure; 
plot(t_prof(1:end-1), mp, 'b', 'LineWidth', 2);
grid on;
xlabel('Time (s)');
ylabel('m_p (g)');
title('Mass of Pellets')

% Plot Inputs
figure;
plot(t_prof, uf_actual, 'b', 'LineWidth', 2); hold on;
plot(t_prof, up_actual, 'r--', 'LineWidth', 2);
grid on;
xlabel('Time (s)');
ylabel('Control Inputs');
legend('u_f','u_p');
title('Control Inputs u_f(t) and u_p(t)');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Part 2.7.6 Robustness
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Cf = 500;       % J/C°
Cc = 300;       % J/C° (decrease causes the system to become unstable at higher temps)
kf = 200;       % W/°C (increased
kfa = 50;       % W/°C (increased)
kca = 50;       % W/°C (increased)
gamma = 1000;   % W/(g/s*fan)
up0 = 0.4;      % g/s
uf0 = 0.6;      % unitless
Tamb = 25;      % °C (decreased) (closer the ambient temp is to the set point the less error in Tf)


% Nonlinear Input
q_fun = @(up, uf) [gamma * up .* uf ./ Cf; 0];

% Nonlinear Matrix
A_nl = [-(kf + kfa)/Cf, kf/Cf; ...
         kf/Cc, -(kf + kca)/Cc];

% Ambient Losses
b_nl = [kfa * Tamb / Cf; ...
        kca * Tamb / Cc];

% Define A, B, C, and D Without Mass Since it Does Not Effect Dynamics
A = [-(kf+kfa)/Cf kf/Cf; ... 
        kf/Cc -(kf+kca)/Cc];
B = [gamma*uf0/Cf gamma*up0/Cf; ...
        0 0];
C = [0 1];
D = 0;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Part 2.7.6 Robustness Nonlinear System
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Cooking Profile
Tc_sp_fun = @(t) ...
    (t >=   0 & t <  600).*90  + ...   % 0–600 s:  90 °C
    (t >= 600 & t < 1200).*110 + ...   % 600–1200 s: 110 °C
    (t >=1200 & t <=1800).*130;        % 1200–1800 s:130 °C

% Setpoint
r_fun = @(t) Tc_sp_fun(t) - Tc0;

% Define New TSPAN
tspan_prof = [0 1800];

% Initial Conditions
z0_prof = [Tf0; Tc0; 0; 0; 0];

% ODE45 Solver
[t_prof, z_prof] = ode45(@(t,z) ...
        f_profile2(t, z, A_nl, b_nl, A, B, C, Kx, KI, L, Kr, q_fun, uf0, up0, Tc0, r_fun), ...
        tspan_prof, z0_prof);

% Results
Tf_prof   = z_prof(:,1);
Tc_prof   = z_prof(:,2);
xI_prof   = z_prof(:,3);
xhat_prof = z_prof(:,4:5);

Tf_hat = xhat_prof(:,1) + Tf0;
Tc_hat = xhat_prof(:,2) + Tc0;

e_Tf = Tf_prof - Tf_hat;
e_Tc = Tc_prof - Tc_hat;


function dz = f_profile2(t, z, A_nl, b_nl, A, B, C, Kx, KI, L, Kr, q_fun, uf0, up0, Tc0, r_fun)
% Closed Loop System with Observer and Varying Set Points
% z = [ Tf; Tc; xI; Tf_hat_dev; Tc_hat_dev ]

% Define States
Tf   = z(1);
Tc   = z(2);
xI   = z(3);
xhat = z(4:5);
x_nl  = [Tf; Tc];
y = Tc - Tc0;

% Time Varying Setpoint
r_t = r_fun(t);

% Define the Input
u = -Kx * xhat - KI * xI + Kr * r_t;

% Individual Inputs
uf = uf0 + u(1);
up = up0 + u(2);

% Nonlinear System
xdot = A_nl * x_nl + b_nl + q_fun(up, uf);

% Integral Action
xIdot = r_t - y;

% Observer System
xhat_dot = A * xhat + B * u + L * (y - C * xhat);

% Augmented System
dz = [xdot;
      xIdot;
      xhat_dot];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Part 2.7.6 Robustness Mass Extrapolation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Extract Number of Points
N = length(t_prof);
% Define U
u_hist = zeros(N, 2);

for k = 1:N
    % xI at Time k
    xI_k   = xI_prof(k);
    % xhat at Time k
    xhat_k = xhat_prof(k, :).';
    % u at Time k
    u_k = -Kx * xhat_k - KI * xI_k + Kr * r;
    % Store Inputs at Time k
    u_hist(k, :) = u_k.';
end

uf_k = u_hist(:,1);
up_k = u_hist(:,2);

% Pellet Input and Fan Input
up_actual = up0 + up_k;
uf_actual = uf0 + uf_k;
% Time Steps
dt = t_prof(2:end) - t_prof(1:end-1);
% Total Mass
mp = mp0 - cumsum(up_actual(1:end-1) .* dt);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Part 2.7.6 Robustness Plots
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Plot Tc with The Profile
figure;
plot(t_prof, Tc_prof, 'b', 'LineWidth', 2);  hold on;
plot(t_prof, Tc_hat,'r--', 'LineWidth', 1.5);
grid on;
xlabel('Time (s)');
ylabel('T_c (°C)');
legend('Actual T_c','Estimated T_c');
title('Cooking Profile Tracking (No Disturbance)');

% Plot the Actual vs Observed Tf with Distrubances
figure;
plot(t_prof, Tf_prof, 'b', 'LineWidth', 2);  hold on;
plot(t_prof, Tf_hat, 'r--', 'LineWidth', 2);
grid on;
xlabel('Time (s)');
ylabel('Temperature (°C)');
title('Actual vs Observed T_f');
legend('Actual T_f', 'Observed T_f');

% Plot the Error 
figure;
plot(t_prof, e_Tf, 'b', 'LineWidth', 1.5); hold on;
plot(t_prof, e_Tc, 'r', 'LineWidth', 1.5);
xlabel('Time (s)');
ylabel('Estimation Error (°C)');
legend('Tf error', 'Tc error');
grid on;
title('Observer Estimation Error');

% Pellet Mass
figure; 
plot(t_prof(1:end-1), mp, 'b', 'LineWidth', 2);
grid on;
xlabel('Time (s)');
ylabel('m_p (g)');
title('Mass of Pellets')

% Plot Inputs
figure;
plot(t_prof, uf_actual, 'b', 'LineWidth', 2); hold on;
plot(t_prof, up_actual, 'r--', 'LineWidth', 2);
grid on;
xlabel('Time (s)');
ylabel('Control Inputs');
legend('u_f','u_p');
title('Control Inputs u_f(t) and u_p(t)');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Part 2.7.7 Advanced Control System
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
