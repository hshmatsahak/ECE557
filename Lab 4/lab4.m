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
A = double(A);
B = double(B);

% Define C matrix
C = [1, 0, 0, 0; 0, 0, 1, 0];

% check Observa
Qo = obsv(A, C)
r = rank(Qo) % should have full rank = 4 to be controllable. check the result = 4, so the system is controllable


% Eigenvalue Assignment
p_desired = [-5, -5.0001, -4.9999, -4.9995];
K_place = place(double(A), double(B), p_desired) % = -14.4305  -17.6456   46.8869    8.3668
%flip the sign of K_place
K = -K_place
(eig(A+B*K));

% Tune p_desired to get  
p_desired = [-7.81, -7.82, -7.79, -7.78];
K_place = place(double(A), double(B), p_desired) 
%flip the sign of K_place
K = -K_place
(eig(A+B*K))

% Tuning K with LQR
% LQR parameters that satisfy prelab req
q1 = 2000;
q2 = 0.5;
R = 0.2; % K_lqr = [-100.00 -56.13 129.01 23.44]
% 

%%% working set in real experiment
% q1 = 8;
% q2 = 5;
% R = 0.2;

%%% testing set
% q1 = 10;
% q2 = 60;
% R = 0.1;

Q = [q1 0 0 0;
      0 0 0 0;
      0 0 q2 0;
      0 0 0 0]

K_lqr = lqr(double(A), double(B), Q, R) 
K = -K_lqr
eig_lqr = double(eig(A - B*K_lqr))

%% L eigenvalues = -10

%Determine L using duality with K
L_p_desired = [-10, -10.0001, -9.999, -10.0002];
A_transpose = A';
C_transpose = C';
K_L_place = place(double(A_transpose), double(C_transpose), L_p_desired);
L = K_L_place'

%Check the value of the assigned eigenvalues.
(eig(A-L*C));

%Assign the ctrl matrices
Actrl = A + B*K - L*C; 
Bctrl = [L -B*K];
Cctrl = K;
Dctrl = [0 0 -K];

%Define the initial condition vector and voltage saturation level. 
x0 = [0; 0; pi/24; 0]
Ulim = 5.875

out = sim('lab4_prep.slx', 30)

% plotting
figure(1)
subplot(311)
title('Comparison of the two controllers with linearized plant')
subtitle('Preparation')
ylabel('z')
hold on
plot(out.z)
legend('state feedback', 'output feedback', 'reference')
subplot(312)
ylabel('theta')
hold on
plot(out.theta)
legend('state feedback', 'output feedback')
subplot(313)
ylabel('u')
hold on
plot(out.u)
legend('state feedback','output feedback','Location','NorthEast')

%% L eigenvalues = -40

%Determine L using duality with K
L_p_desired = [-40, -40.0001, -39.999, -40.0002];
A_transpose = A';
C_transpose = C';
K_L_place = place(double(A_transpose), double(C_transpose), L_p_desired);
L = K_L_place'

%Check the value of the assigned eigenvalues.
(eig(A-L*C));

%Assign the ctrl matrices
Actrl = A + B*K - L*C; 
Bctrl = [L -B*K];
Cctrl = K;
Dctrl = [0 0 -K];

%Define the initial condition vector and voltage saturation level. 
x0 = [0; 0; pi/24; 0]
Ulim = 5.875

out = sim('lab4_prep.slx', 30)

% plotting
figure(1)
subplot(311)
title('Comparison of the two controllers with linearized plant')
subtitle('Preparation')
ylabel('z')
hold on
plot(out.z)
legend('state feedback', 'output feedback', 'reference')
subplot(312)
ylabel('theta')
hold on
plot(out.theta)
legend('state feedback', 'output feedback')
subplot(313)
ylabel('u')
hold on
plot(out.u)
legend('state feedback','output feedback','Location','NorthEast')

%% 5.2 - Output Feedback Control with integral action 
% define Abar and Bbar
Abar = [A zeros(4,1); -C1 0];
Bbar = [B; 0];
C1 = [1 0 0 0];

Q_c = ctrb(Abar, Bbar); % controllability matrix
rank(Q_c) % check that the system is controllable with full rank = 5

Q = [0.2 0 0 0 0; 
     0 0 0 0 0; 
     0 0 50 0 0; 
     0 0 0 0 0; 
     0 0 0 0 3];
R = 0.05;

Kbar = -lqr(Abar, Bbar, Q, R); % need to flip the sign
Kx = Kbar(1:end-1); % first 4 components
Kxi = Kbar(end); % last component

eigenvals = eig(Abar + Bbar * Kbar)
% all these eigenvalues are in the LHP
