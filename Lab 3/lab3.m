%% 4.2 Symbolic Linearization and Controllability Analysis

% define numerical variables 
M = 1.1911;
m = 0.23;
l = 0.3302;
g = 9.8;
alpha1 = 1.7265; 
alpha2 = 10.5;

% syms z, z_dot,theta,theta_dot,u real
syms z real
syms z_dot real
syms theta real
syms theta_dot real
syms u real

% define non-linear system
z_ddot = (-m*l*sin(theta)*theta_dot^2 + m*g*sin(theta)*cos(theta) + alpha1*u - alpha2*z_dot) / (M + m*sin(theta)^2);
theta_ddot = (-m*l*sin(theta)*cos(theta)*theta_dot^2 + (M + m)*g*sin(theta) + (alpha1*u - alpha2*z_dot)*cos(theta)) / (l*(M + m*sin(theta)^2));

% define state and its derivative 
x = [z;z_dot;theta;theta_dot]
f = [z_dot;z_ddot;theta_dot;theta_ddot]

matlabFunction(eval(f), 'Vars', {x, u}, 'File', 'cartpend', 'Outputs',{'xdot'});
test_cartpend = cartpend([0 0 pi/2 0]', 1) % should be [0; 1.2149; 0; 21.679], there may be a typo in the lab document

% equilibrium state
x_bar = zeros(4, 1);
u_bar = 0;

A = subs(jacobian(f,x), [x;u], [x_bar;u_bar])
B = subs(jacobian(f,u), [x;u], [x_bar;u_bar])

% convert A and B entries to float for later use
% A = double(A);
% B = double(B);

% check controllability
Co = ctrb(A, B)

r = rank(Co) % should have full rank = 4 to be controllable. check the result = 4, so the system is controllable

%% 4.3 Eigenvalue Assignment
p_desired = [-5.0, -5.0, -5.0, -5.0];
K_acker = acker(double(A), double(B), p_desired) % = -14.5282  -17.7043   47.0341    8.3938

p_check = eval(eig(A - B*K_acker)) % this should be the same as the K from acker()

% use place()
p_perturbed = [-5, -5.0001, -4.9999, -4.9995];
K_place = place(double(A), double(B), p_perturbed) % = -14.2226  -17.5209   46.5737    8.3093
p_place_check = double(eig(A - B*K_place)) % this should get an almost identical result considering the perturbation on the desired eigenvalues

%% 4.4 Linear Quadratic Optimal Control

% Case 1
q1_1 = 2;
q2_1 = 1;
R_1 = 0.2;

Q_1 = [q1_1 0 0 0;
        0 0 0 0;
        0 0 q2_1 0;
        0 0 0 0]

K_lqr1 = lqr(double(A), double(B), Q_1, R_1) % = -3.1623  -13.7010   40.7104    7.2175

eig_lqr1 = double(eig(A - B*K_lqr1))

% Case 2
q1_2 = 0.1;
q2_2 = 1;
R_2 = 1;

Q_2 = [q1_2 0 0 0;
        0 0 0 0;
        0 0 q2_2 0;
        0 0 0 0]

K_lqr2 = lqr(double(A), double(B), Q_2, R_2) % -0.3162  -12.3180   37.5935    6.6541

eig_lqr2 = double(eig(A - B*K_lqr2))

% we compare eig_lqr1 with eig_lqr2
% eig_lqr1 = [-9.56, -6.04, -4.52, -0.52]
% eig_lqr2 = [-9.66, -5.74, -4.72, -0.05]
% the two systems have different eigenvalues, especially the final eigenvalue
% -0.52 vs. -0.052, they are approximately a factor of 10. The other 3
% eigenvalues are similar.

%% 4.5 Simulation of Linear Controllers and Linearized Plant

% linearization matrices
A = double(A)
B = double(B)
x0 = 2*[0; 0; pi/24; 0]; % initial condition

K = -K_acker; % K_acker is obtained using the acker() function and K_place is obtained using place() function
% We use K_acher here
out_place = sim('lab3_linear.slx', 10)

K = -K_lqr1;
out_LQR1 = sim('lab3_linear.slx', 10)

K = -K_lqr2;
out_LQR2 = sim('lab3_linear.slx', 10)

% plotting
figure(1)
subplot(311)
title('Comparison of the three controllers with linearized plant')
subtitle('Preparation')
ylabel('z')
hold on
plot(out_place.z)
plot(out_LQR1.z)
plot(out_LQR2.z)

subplot(312)
ylabel('theta')
hold on
plot(out_place.theta)
plot(out_LQR1.theta)
plot(out_LQR2.theta)

subplot(313)
ylabel('u')
hold on
plot(out_place.u)
plot(out_LQR1.u)
plot(out_LQR2.u)
legend('pole assignment','LQR1','LQR2','Location','NorthEast')

%% 5
Ulim = 5.875;
K = -K_acker;
x0 = [0; 0; pi/24; 0];

%% 5.1 no saturation pole assignment
Ulim = inf;
K = -K_acker;
n = 7; % start with 1, increment by 1
x0 = [0; 0; n*pi/24; 0];

out_place = sim('lab3_nonlinear.slx', 10)

% plotting
fig = figure(1);
fontsize(fig, scale=3)

subplot(311)
title('Comparison linear vs. nonlinear: pole assignment')
subtitle('z(t)')
ylabel('z [m]')
hold on
plot(out_place.z)

subplot(312)
subtitle('theta(t)')
ylabel('theta [rad]')
hold on
plot(out_place.theta)

subplot(313)
subtitle('u(t)')
ylabel('u [V]')
hold on
plot(out_place.u)

xlabel('t [s]')
legend('linear', 'nonlinear','Location','NorthEast')

%% 5.2 no sat LQR 1

Ulim = inf;
K = -K_lqr1;
n = 8; % start with 1, increment by 1
x0 = [0; 0; n*pi/24; 0];

out_place = sim('lab3_nonlinear.slx', 10)

% plotting
fig = figure(1);
fontsize(fig, scale=3)

subplot(311)
title('Comparison linear vs. nonlinear: K_L_Q_R_1')
subtitle('z(t)')
ylabel('z [m]')
hold on
plot(out_place.z)

subplot(312)
subtitle('theta(t)')
ylabel('theta [rad]')
hold on
plot(out_place.theta)

subplot(313)
subtitle('u(t)')
ylabel('u [V]')
hold on
plot(out_place.u)

xlabel('t [s]')
legend('linear', 'nonlinear','Location','NorthEast')

%% 5.3 no sat LQR 2

Ulim = inf;
K = -K_lqr2;
n = 9; % start with 1, increment by 1
x0 = [0; 0; n*pi/24; 0];

out_place = sim('lab3_nonlinear.slx', 10)

% plotting
fig = figure(1);
fontsize(fig, scale=3)

subplot(311)
title('Comparison linear vs. nonlinear: K_L_Q_R_2')
subtitle('z(t)')
ylabel('z [m]')
hold on
plot(out_place.z)

subplot(312)
subtitle('theta(t)')
ylabel('theta [rad]')
hold on
plot(out_place.theta)

subplot(313)
subtitle('u(t)')
ylabel('u [V]')
hold on
plot(out_place.u)

xlabel('t [s]')
legend('linear', 'nonlinear','Location','NorthEast')

%% 5.4 sat pole assignment

Ulim = 5.875;
K = -K_acker;
n = 2; % start with 1, increment by 1
x0 = [0; 0; n*pi/24; 0];

out_place = sim('lab3_nonlinear.slx', 10)

% plotting
fig = figure(1);
fontsize(fig, 50, "points")

subplot(311)
title('Comparison linear vs. nonlinear: K_p_l_a_c_e')
subtitle('z(t)')
ylabel('z [m]')
hold on
plot(out_place.z)

subplot(312)
subtitle('theta(t)')
ylabel('theta [rad]')
hold on
plot(out_place.theta)

subplot(313)
subtitle('u(t)')
ylabel('u [V]')
hold on
plot(out_place.u)

xlabel('t [s]')
legend('linear', 'nonlinear','Location','NorthEast')

%% 5.5 sat LQR 1

Ulim = 5.875;
K = -K_lqr1;
n = 2; % start with 1, increment by 1
x0 = [0; 0; n*pi/24; 0];

out_place = sim('lab3_nonlinear.slx', 10)

% plotting
fig = figure(1);
fontsize(fig, scale=3)

subplot(311)
title('Comparison linear vs. nonlinear: K_L_Q_R_1')
subtitle('z(t)')
ylabel('z [m]')
hold on
plot(out_place.z)

subplot(312)
subtitle('theta(t)')
ylabel('theta [rad]')
hold on
plot(out_place.theta)

subplot(313)
subtitle('u(t)')
ylabel('u [V]')
hold on
plot(out_place.u)

xlabel('t [s]')
legend('linear', 'nonlinear','Location','NorthEast')

%% 5.6 sat LQR 2

Ulim = 5.875;
K = -K_lqr2;
n = 2; % start with 1, increment by 1
x0 = [0; 0; n*pi/24; 0];

out_place = sim('lab3_nonlinear.slx', 10)

% plotting
fig = figure(1);
fontsize(fig, scale=3)

subplot(311)
title('Comparison linear vs. nonlinear: K_L_Q_R_2')
subtitle('z(t)')
ylabel('z [m]')
hold on
plot(out_place.z)

subplot(312)
subtitle('theta(t)')
ylabel('theta [rad]')
hold on
plot(out_place.theta)

subplot(313)
subtitle('u(t)')
ylabel('u [V]')
hold on
plot(out_place.u)

xlabel('t [s]')
legend('linear', 'nonlinear','Location','NorthEast')

