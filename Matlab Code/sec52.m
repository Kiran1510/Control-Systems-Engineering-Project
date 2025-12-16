%% section 5.2: state estimation from output along the same open-loop trajectory as 5.1
clear; clc; close all;

% make all plot text use LaTeX
set(groot,'defaultTextInterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');
set(groot,'defaultAxesTickLabelInterpreter','latex');

% reduced-order CLTI system
A = [0.1883  1.4362;
     0.3178  2.7380];
B = [0;
     0.0146];
C = [0 1];

% operating point
Vx_star     = 20;
lambda_star = 0.2;
x_star      = [Vx_star; lambda_star];
u_star      = 725.398;

% time horizon + discretization
T  = 2.0;
dt = 1e-3;
t  = 0:dt:T;
N  = numel(t);

% initial condition (same as 5.1)
x0_abs = [21; 0.35];
x0     = x0_abs - x_star;
xf = [0; 0];

% minimum-energy open-loop input
W = zeros(2,2);
for k = 1:N
    s = t(k);
    E = expm(A*s);
    W = W + (E*(B*B')*E')*dt;
end
Winv = pinv(W);

delta = (xf - expm(A*T)*x0);

u_tilde = zeros(1,N);
for k = 1:N
    tau = t(k);
    u_tilde(k) = (B' * expm(A'*(T - tau)) * (Winv * delta));
end

% observer design
p_obs = [-10+4.8j, -10-4.8j];
L = place(A', C', p_obs).';

% simulate true plant + observer
x_true = zeros(2,N);
x_hat  = zeros(2,N);

x_true(:,1) = x0;

% intentionally wrong initial estimate
xhat0_abs  = [25; 0.05];
x_hat(:,1) = xhat0_abs - x_star;

for k = 1:N-1
    % true output (deviation)
    y_tilde = C*x_true(:,k);

    % true plant dynamics
    xdot = A*x_true(:,k) + B*u_tilde(k);
    x_true(:,k+1) = x_true(:,k) + xdot*dt;

    % observer dynamics
    yhat_tilde = C*x_hat(:,k);
    xhat_dot = A*x_hat(:,k) + B*u_tilde(k) + L*(y_tilde - yhat_tilde);
    x_hat(:,k+1) = x_hat(:,k) + xhat_dot*dt;
end

% convert to absolute for plotting
x_true_abs = x_true + x_star;
x_hat_abs  = x_hat  + x_star;

Vx     = x_true_abs(1,:);
lam    = x_true_abs(2,:);
Vxhat  = x_hat_abs(1,:);
lamhat = x_hat_abs(2,:);

eVx  = Vx - Vxhat;
elam = lam - lamhat;

% plotting results
figure('Position', [100 100 1000 700]);

subplot(2,2,1);
plot(t, Vx, 'LineWidth', 2); hold on;
plot(t, Vxhat, '--', 'LineWidth', 2);
yline(Vx_star, '--', 'LineWidth', 1.2);
grid on;
xlabel('Time (s)'); ylabel('$V_x$ (m/s)');
title('$V_x$: true vs estimated');
legend('$V_x$', '$\hat{V}_x$', '$V_x^\ast$', 'Location', 'best');

subplot(2,2,2);
plot(t, lam, 'LineWidth', 2); hold on;
plot(t, lamhat, '--', 'LineWidth', 2);
yline(lambda_star, '--', 'LineWidth', 1.2);
grid on;
xlabel('Time (s)'); ylabel('$\lambda$');
title('$\lambda$: true vs estimated');
legend('$\lambda$', '$\hat{\lambda}$', '$\lambda^\ast$', 'Location', 'best');

subplot(2,2,3);
plot(t, eVx, 'LineWidth', 2); hold on;
yline(0, '--', 'LineWidth', 1.2);
grid on;
xlabel('Time (s)'); ylabel('$e_{V_x}=V_x-\hat{V}_x$');
title('Estimation error in $V_x$');

subplot(2,2,4);
plot(t, elam, 'LineWidth', 2); hold on;
yline(0, '--', 'LineWidth', 1.2);
grid on;
xlabel('Time (s)'); ylabel('$e_{\lambda}=\lambda-\hat{\lambda}$');
title('Estimation error in $\lambda$');

sgtitle('Section 5.2: Observer estimates along the SAME open-loop trajectory as 5.1', ...
    'FontWeight', 'bold');