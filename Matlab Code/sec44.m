%% section 4.4: apply reduced-order observer/controller to original nonlinear system
clear; clc; close all;

fprintf('=== section 4.4: apply observer/controller to nonlinear (CNL) plant ===\n\n');

% vehicle and tire parameters
m  = 342;
Jw = 1.13;
R  = 0.33;
g  = 9.81;
FN = m*g;

% friction curve parameters
c1 = 1.2801;
c2 = 23.99;
c3 = 0.52;
c4 = 0.03;

% operating point
lambda_star = 0.2;
Vx_star     = 20;

mu_star = mu_friction(lambda_star, Vx_star, c1,c2,c3,c4);
u_star  = mu_star * FN * R * ( 1 + (1-lambda_star)*Jw/(m*R^2) );

fprintf('operating point:\n');
fprintf('  Vx* = %.3f m/s, lambda* = %.3f\n', Vx_star, lambda_star);
fprintf('  mu* = %.4f\n', mu_star);
fprintf('  u*  = %.3f N*m\n\n', u_star);

% reduced-order linear subsystem
Ao = [0.1883  1.4362;
      0.3178  2.7380];

Bo = [0;
      0.0146];

Co = [0 1];

% controller design
ts   = 2.0;
zeta = 0.9;
wn   = 4/(zeta*ts);

pc1 = -zeta*wn + 1i*wn*sqrt(1-zeta^2);
pc2 = -zeta*wn - 1i*wn*sqrt(1-zeta^2);

K = place(Ao, Bo, [pc1 pc2]);
Acl = Ao - Bo*K;

fprintf('controller design:\n');
fprintf('  ts=%.2f, zeta=%.2f -> wn=%.3f\n', ts, zeta, wn);
fprintf('  K = [%.3f  %.3f]\n', K(1), K(2));
pcl = eig(Acl);
fprintf('  controller poles: %.6f%+.6fj, %.6f%+.6fj\n\n', real(pcl(1)), imag(pcl(1)), real(pcl(2)), imag(pcl(2)));

% observer design
alpha = 5;
pobs  = alpha * pcl;
L     = place(Ao.', Co.', pobs).';

fprintf('observer design:\n');
fprintf('  alpha = %.2f\n', alpha);
fprintf('  L = [%.3f; %.3f]\n', L(1), L(2));
pob = eig(Ao - L*Co);
fprintf('  observer poles: %.6f%+.6fj, %.6f%+.6fj\n\n', real(pob(1)), imag(pob(1)), real(pob(2)), imag(pob(2)));

% actuator saturation
u_min = 0;
u_max = 1200;

% simulation setup
t_sim = 5.0;
dt    = 0.001;
t     = 0:dt:t_sim;
N     = numel(t);

% initial conditions
Sx0 = 0;
Vx0 = 19.5;
lam0 = 0.30;

x = zeros(3, N);
x(:,1) = [Sx0; Vx0; lam0];

% observer state (deviation coordinates)
xhat = zeros(2, N);
xhat(:,1) = [0; 0];  % start with zero deviation estimate

% DEBUG: Check initial values
fprintf('=== DEBUG: Initial Conditions ===\n');
fprintf('True initial (absolute): Vx0=%.2f, lam0=%.2f\n', Vx0, lam0);
fprintf('Observer initial (deviation): xhat=[%.2f; %.2f]\n', xhat(1,1), xhat(2,1));
fprintf('Observer initial (absolute): Vxhat=%.2f, lamhat=%.2f\n', ...
    Vx_star+xhat(1,1), lambda_star+xhat(2,1));
fprintf('Expected Vxhat_tr(1) should be: %.2f\n\n', Vx_star+xhat(1,1));

u  = zeros(1, N);

% closed-loop simulation
for k = 1:N-1
    % measurement from nonlinear plant
    y      = x(3,k);
    ytilde = y - lambda_star;

    % current estimated deviation state
    dxhat = xhat(:,k);

    % state feedback using estimated deviation state
    du_k = -K * dxhat;
    u_k  = u_star + du_k;

    % saturation
    u_k  = min(max(u_k, u_min), u_max);
    du_k = u_k - u_star;

    u(k)  = u_k;

    % nonlinear plant dynamics
    Vx  = x(2,k);
    lam = x(3,k);

    Vx_safe = max(Vx, 0.1);

    mu = mu_friction(lam, Vx_safe, c1,c2,c3,c4);

    Sx_dot = Vx;
    Vx_dot = -mu * g;

    w = (1 - lam) * Vx / R;
    w_dot = (mu*FN*R - u_k) / Jw;
    lam_dot = -(R/Vx_safe)*w_dot + (R*w/(Vx_safe^2))*Vx_dot;

    xdot = [Sx_dot; Vx_dot; lam_dot];

    % integrate nonlinear plant
    x(:,k+1) = x(:,k) + xdot*dt;

    % linear observer (deviation coordinates)
    innov = ytilde - Co*dxhat;
    dxhat_dot = Ao*dxhat + Bo*du_k + L*innov;
    xhat(:,k+1) = xhat(:,k) + dxhat_dot*dt;
end

u(end) = u(end-1);

% DEBUG: Check after simulation
fprintf('=== DEBUG: After Simulation ===\n');
fprintf('xhat(1,1) = %.6f (should be 0)\n', xhat(1,1));
fprintf('xhat(2,1) = %.6f (should be 0)\n', xhat(2,1));
fprintf('xhat(1,end) = %.6f\n', xhat(1,end));
fprintf('xhat(2,end) = %.6f\n', xhat(2,end));
fprintf('\n');

% build estimated full trajectories for plotting
Vx_tr   = x(2,:);
lam_tr  = x(3,:);

Vxhat_tr  = Vx_star + xhat(1,:);
lamhat_tr = lambda_star + xhat(2,:);

% DEBUG: Check plotting arrays
fprintf('=== DEBUG: Plotting Arrays ===\n');
fprintf('Vxhat_tr(1) = %.6f (should be %.2f)\n', Vxhat_tr(1), Vx_star);
fprintf('Vxhat_tr(end) = %.6f (should be near %.2f)\n', Vxhat_tr(end), Vx_star);
fprintf('lamhat_tr(1) = %.6f (should be %.2f)\n', lamhat_tr(1), lambda_star);
fprintf('lamhat_tr(end) = %.6f (should be near %.2f)\n\n', lamhat_tr(end), lambda_star);

% performance metrics
tol = 0.02 * abs(lambda_star);

idx_settle = find(abs(lam_tr - lambda_star) <= tol, 1, 'first');
if isempty(idx_settle)
    ts_actual = NaN;
else
    ts_actual = t(idx_settle);
end

overshoot = (max(lam_tr) - lambda_star) / lambda_star * 100;
ss_err = abs(lam_tr(end) - lambda_star);

fprintf('performance (nonlinear + observer-based):\n');
fprintf('  overshoot = %.2f %%\n', overshoot);
fprintf('  ss error  = %.6f\n', ss_err);
if isnan(ts_actual)
    fprintf('  ts(2%%)    = NaN (not within window)\n');
else
    fprintf('  ts(2%%)    = %.3f s\n', ts_actual);
end

fprintf('\ncontrol effort:\n');
fprintf('  u* = %.2f N*m\n', u_star);
fprintf('  u range = [%.2f, %.2f] N*m\n', min(u), max(u));
fprintf('  sat = [%.1f, %.1f] N*m\n\n', u_min, u_max);

% add comparison to section 4.3
fprintf('comparison to section 4.3 (linear plant):\n');
fprintf('  4.3 overshoot: 171.36%%, 4.4 overshoot: %.2f%%\n', overshoot);
fprintf('  4.3 ts: 2.907s, 4.4 ts: %.3fs\n', ts_actual);
if ts_actual < 2.907
    fprintf('  improvement: %.2f%% less overshoot, %.3fs faster\n\n', 171.36-overshoot, 2.907-ts_actual);
else
    fprintf('  difference: %.2f%% less overshoot, %.3fs slower\n\n', 171.36-overshoot, ts_actual-2.907);
end

% plotting results - clean 2x2 layout
figure('Position',[100 100 1200 700]);

% velocity tracking
subplot(2,2,1);
plot(t, Vx_tr, 'LineWidth', 2); hold on;
plot(t, Vxhat_tr, '--', 'LineWidth', 2);
yline(Vx_star, 'r--', 'LineWidth', 1.5);
grid on;
xlabel('Time (s)'); ylabel('Velocity Vx (m/s)');
title('Velocity: true vs estimated');
legend('Vx (true)', 'Vxhat (estimated)', 'Vx* (equilibrium)', 'Location','best');

% slip ratio tracking
subplot(2,2,2);
plot(t, lam_tr, 'LineWidth', 2); hold on;
plot(t, lamhat_tr, '--', 'LineWidth', 2);
yline(lambda_star, 'r--', 'LineWidth', 1.5);
yline(lambda_star + tol, 'g:', 'LineWidth', 1.2);
yline(lambda_star - tol, 'g:', 'LineWidth', 1.2);
grid on;
xlabel('Time (s)'); ylabel('Slip ratio lambda');
title('Slip ratio vs Time (Nonlinear closed-loop)');
legend('lambda (true)', 'lambdahat (estimated)', 'lambda* (equilibrium)', 'pm 2 percent band', 'Location','best');

% control torque
subplot(2,2,3);
plot(t, u, 'LineWidth', 2); hold on;
yline(u_star, 'r--', 'LineWidth', 1.5);
yline(u_min, 'k:', 'LineWidth', 1.2);
yline(u_max, 'k:', 'LineWidth', 1.2);
grid on;
xlabel('Time (s)'); ylabel('Brake torque u (Nm)');
title('Control input (with saturation)');
legend('u(t)', 'u*', 'limits', 'Location','best');

% tracking error
subplot(2,2,4);
e = lam_tr - lambda_star;
plot(t, e, 'LineWidth', 2); hold on;
yline(0, 'r--', 'LineWidth', 1.5);
yline(tol, 'g:', 'LineWidth', 1.2);
yline(-tol, 'g:', 'LineWidth', 1.2);
grid on;
xlabel('Time (s)'); ylabel('Error e = lambda - lambda*');
title('Tracking error');
legend('error', '0', 'pm 2 percent band', 'Location','best');

sgtitle('Section 4.4: Observer-based Controller Applied to Nonlinear ABS Plant', ...
    'FontSize', 14, 'FontWeight','bold');

%% local function
function mu = mu_friction(lambda, Vx, c1,c2,c3,c4)
    mu = (c1*(1-exp(-c2*lambda)) - c3*lambda) * exp(-c4*Vx);
end