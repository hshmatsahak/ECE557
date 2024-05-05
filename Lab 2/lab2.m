%% Define A, B, C and the Physical Parameters
alpha_1 = 1.7265; % units: N/V
alpha_2 = 10.5; % units: Ns/m
M = 0.8211; % units: Kg

A = [0 1; 0 -alpha_2/M];
B = [0; alpha_1/M];
C = [1 0];

%% 4.2 State Feedback Stablization
eigval_1 = [-5 -5]; % T_s1 = 1.167s
eigval_2 = [-5 -5]; % T_s2 = 2.921s; T_s1 / T_s2 = 0.4
eig_test = eigval_1; % the eigenvalues to test

[K1 K2] = state_feedback_design(A, B, eig_test);
K = [K1 K2]; % obtain K

eigenVal = eig(A+B*K);
disp(eigenVal == eig_test.'); % test if the eigenvalues are correct

%% 4.4 State Estimation and Output Feedback Control
eigval_3 = [-10 -10];
[L1,L2] = observer_design(A, eigval_3);
L = [L1 ; L2]; % obtain L

eigenVal_sys = eig(A - L*C);
disp(int8(eigenVal_sys) == eigval_3.'); % test if the eigenvalues are correct

%% Output Feedback Controller
p_feedback = [-50 -50]; % desired eigenvalues for (A + BK)
p_observer = [-50 -50]; % desired eigenvalues for (A - LC)

[Actrl, Bctrl, Cctrl, Dctrl] = output_feedback_controller(A, B, C, p_feedback, p_observer);

% display the 4 controller matrices
fprintf("Actrl= \n");
disp(Actrl);
fprintf("Bctrl= \n");
disp(Bctrl);
fprintf("Cctrl= \n");
disp(Cctrl);
fprintf("Dctrl= \n");
disp(Dctrl);

%% Functions
function [K1 K2] = state_feedback_design(A, B, p)
    % Inputs: matrix A, B; row vector p with the desired eigenvalues 
    % Outputs: gains K1, K2 for the state feedback controller

    K1 = -p(1)*p(2)/B(2);
    K2 = (-A(2,2) + p(1) + p(2))/B(2);
end

function [L1,L2] = observer_design(A, p)
    % Inputs: matrix A, row vector p with desired eigenvalues
    % Outputs: observer gains L1, L2

    L1 = A(2, 2) - (p(1) + p(2));
    L2 = L1*A(2, 2) + p(1)*p(2);
end

function [Actrl, Bctrl, Cctrl, Dctrl] = output_feedback_controller(A, B, C, p_feedback, p_observer)
    % Inputs: matrix A, B, C, and p_feedback: desired eigenvalues for
    % (A + BK), p_observer: desired eigenvalues for (A - LC)
    % Outputs: 4 controller matrices: Actrl, Bctrl, Cctrl, Dctrl

    [K1 K2] = state_feedback_design(A, B, p_feedback); % use p_feedback for state feedback
    K = [K1 K2]; % obtain state gains K1, K2

    [L1, L2] = observer_design(A, p_observer); % use p_observer for observer 
    L = [L1 ; L2]; % obtain observer gains L1, L2

    % calculate the 4 controller matrices
    Actrl = A + B*K - L*C; % 2*2
    Bctrl = [L -B*K]; % 2*3
    Cctrl = K; % 1*2
    Dctrl = [0 -K]; % 1*3
end
