%% section 4.3: closed-loop estimator/controller (observer-based state feedback)
clear; clc; close all;

fprintf('=== section 4.3: observer-based state feedback ===\n\n');

% vehicle and tire parameters
m  = 342;
Jw = 1.13;
Rr = 0.33;
g  = 9.81;
FN = m*g;

% friction model parameters
c1 = 1.2801;
c2 = 23.99;
c3 = 0.52;
c4 = 0.03;

% equilibrium values
lambda_star = 0.2;
Vx_star     = 20;

mu_star = (c1*(1-exp(-c2*lambda_star)) - c3*lambda_star) * exp(-c4*Vx_star);
u_star = mu_star * FN * Rr * (1 + (1-lambda_star)*Jw/(m*Rr^2));

fprintf('u* = %.3f N*m, mu* = %.4f\n\n', u_star, mu_star);

% reduced-order observable subsystem: x_o = [Vx; lambda], y = lambda
Ao = [0.1883  1.4362;
      0.3178  2.7380];

Bo = [0;
      0.0146];

Co = [0 1];

% controller gain
K = [294.887  474.404];

Acl = Ao - Bo*K;
pc  = eig(Acl);

fprintf('K = [%.3f  %.3f]\n', K(1), K(2));
fprintf('controller poles:\n');
fprintf('  %.6f %+ .6fj\n', real(pc(1)), imag(pc(1)));
fprintf('  %.6f %+ .6fj\n\n', real(pc(2)), imag(pc(2)));

% observer gain
L = [401.871;
     22.926];

Aobs = Ao - L*Co;
pob  = eig(Aobs);

fprintf('L = [%.3f; %.3f]\n', L(1), L(2));
fprintf('observer poles:\n');
fprintf('  %.6f %+ .6fj\n', real(pob(1)), imag(pob(1)));
fprintf('  %.6f %+ .6fj\n\n', real(pob(2)), imag(pob(2)));

% control saturation limits
u_min = 0;
u_max = 1200;

% simulation setup
t_sim = 5.0;
dt    = 0.001;
t     = 0:dt:t_sim;
N     = numel(t);

% initial conditions
x0 = [19.5;    % initial velocity
      0.30];   % initial slip ratio

xhat0 = [Vx_star;
         x0(2)];

% allocate arrays
x    = zeros(2, N);
xhat = zeros(2, N);
u    = zeros(1, N);
y    = zeros(1, N);

x(:,1)    = x0;
xhat(:,1) = xhat0;

% simulation loop
x_star = [Vx_star; lambda_star];

for k = 1:N-1
    % measurement
    y(k) = Co*x(:,k);

    % control law using estimated state
    du_cmd = -K*(xhat(:,k) - x_star);
    u(k) = u_star + du_cmd;

    % apply saturation
    u(k) = min(max(u(k), u_min), u_max);

    % true plant dynamics
    xdot = Ao*(x(:,k) - x_star) + Bo*(u(k) - u_star);

    % observer dynamics
    yhat   = Co*xhat(:,k);
    xhatdot = Ao*(xhat(:,k) - x_star) + Bo*(u(k) - u_star) + L*(y(k) - yhat);

    % integrate states
    x(:,k+1)    = x(:,k)    + xdot*dt;
    xhat(:,k+1) = xhat(:,k) + xhatdot*dt;
end

y(end) = Co*x(:,end);
u(end) = u(end-1);

% performance metrics
lambda_traj = x(2,:);
lambda_ss   = lambda_star;

% settling time (2% band)
tol = 0.02*abs(lambda_ss);
idx = find(abs(lambda_traj - lambda_ss) <= tol, 1, 'first');
if isempty(idx)
    ts_actual = NaN;
else
    ts_actual = t(idx);
end

% overshoot and steady-state error
overshoot = (max(lambda_traj) - lambda_ss)/lambda_ss*100;
ss_error = abs(lambda_traj(end) - lambda_ss);

% torque range
umin = min(u);
umax = max(u);

fprintf('performance (observer-based):\n');
fprintf('  overshoot = %.2f %%\n', overshoot);
fprintf('  ss error  = %.6f\n', ss_error);
if isnan(ts_actual)
    fprintf('  ts(2%%)    = NaN (not within window)\n\n');
else
    fprintf('  ts(2%%)    = %.3f s\n\n', ts_actual);
end

fprintf('control effort:\n');
fprintf('  u* = %.2f N*m\n', u_star);
fprintf('  u range = [%.2f, %.2f] N*m\n', umin, umax);
fprintf('  saturation = [%.1f, %.1f] N*m\n\n', u_min, u_max);

% add comparison to section 4.1
fprintf('comparison to section 4.1 (perfect state feedback):\n');
fprintf('  4.1 overshoot: 176.56%%, 4.3 overshoot: %.2f%%\n', overshoot);
fprintf('  4.1 ts: 2.681s, 4.3 ts: %.3fs\n', ts_actual);
fprintf('  difference: %.2f%% less overshoot, %.3fs slower\n\n', 176.56-overshoot, ts_actual-2.681);

% plotting results - clean 2x2 layout
figure('Position', [100 100 1200 700]);

% velocity tracking
subplot(2,2,1);
plot(t, x(1,:), 'LineWidth', 2); hold on;
plot(t, xhat(1,:), '--', 'LineWidth', 2);
yline(Vx_star, 'r--', 'LineWidth', 1.5);
grid on;
xlabel('Time (s)');
ylabel('Velocity Vx (m/s)');
title('Velocity: true vs estimated');
legend('Vx (true)', 'Vxhat (estimated)', 'Vx* (equilibrium)', 'Location', 'best');

% slip ratio tracking
subplot(2,2,2);
plot(t, x(2,:), 'LineWidth', 2); hold on;
yline(lambda_star, 'r--', 'LineWidth', 1.5);
yline(lambda_star+tol, 'g:', 'LineWidth', 1.2);
yline(lambda_star-tol, 'g:', 'LineWidth', 1.2);
grid on;
xlabel('Time (s)');
ylabel('Slip Ratio lambda');
title('Slip tracking (2 percent band)');
legend('lambda', 'lambda*', 'pm 2 percent band', 'Location', 'best');

% control input
subplot(2,2,3);
plot(t, u, 'LineWidth', 2); hold on;
yline(u_star, 'r--', 'LineWidth', 1.5);
yline(u_min, 'k:', 'LineWidth', 1.2);
yline(u_max, 'k:', 'LineWidth', 1.2);
grid on;
xlabel('Time (s)');
ylabel('Torque u (Nm)');
title('Control input (torque)');
legend('u(t)', 'u*', 'limits', 'Location', 'best');

% tracking error
subplot(2,2,4);
e = x(2,:) - lambda_star;
plot(t, e, 'LineWidth', 2); hold on;
yline(0, 'r--', 'LineWidth', 1.5);
yline(tol, 'g:', 'LineWidth', 1.2);
yline(-tol, 'g:', 'LineWidth', 1.2);
grid on;
xlabel('Time (s)');
ylabel('error e = lambda - lambda*');
title('Tracking error');
legend('error', '0', 'pm 2 percent band', 'Location', 'best');

sgtitle('Section 4.3: Closed-loop Estimator/Controller (Observer-based State Feedback)', ...
    'FontSize', 14, 'FontWeight', 'bold');