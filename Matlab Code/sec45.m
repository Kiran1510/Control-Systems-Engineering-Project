%% section 4.5: disturbance/parameter variation on CNL plant
clear; clc; close all;

fprintf('=== section 4.5: nonlinear plant + friction drop disturbance ===\n\n');

% vehicle and tire parameters
m  = 342;
Jw = 1.13;
R  = 0.33;
g  = 9.81;
FN = m*g;

% nominal (dry asphalt) friction parameters
c_dry = struct('c1',1.2801,'c2',23.99,'c3',0.52,'c4',0.03);

% disturbed (wet road) friction parameters
c_wet = struct('c1',0.70,'c2',15.00,'c3',0.35,'c4',0.03);

t_switch = 0.5;     % seconds: dry -> wet

% operating point
lambda_star = 0.2;
Vx_star     = 20;

mu_star = mu_friction(lambda_star, Vx_star, c_dry);
u_star = mu_star * FN * R * (1 + (1-lambda_star)*Jw/(m*R^2));

fprintf('operating point (designed about dry road):\n');
fprintf('  Vx* = %.3f m/s, lambda* = %.3f\n', Vx_star, lambda_star);
fprintf('  mu*(dry) = %.4f\n', mu_star);
fprintf('  u* = %.3f N*m\n\n', u_star);

% reduced-order controller/observer
Ao = [0.1883  1.4362;
      0.3178  2.7380];
Bo = [0; 0.0146];
Co = [0 1];

% controller design
ts   = 2.0;
zeta = 0.9;
wn   = 4/(zeta*ts);

pc1 = -zeta*wn + 1i*wn*sqrt(1-zeta^2);
pc2 = -zeta*wn - 1i*wn*sqrt(1-zeta^2);
K   = place(Ao, Bo, [pc1 pc2]);

Acl = Ao - Bo*K;
p_cl = eig(Acl);

% observer design
alpha = 5.0;
p_obs = alpha * p_cl;
L = place(Ao', Co', p_obs).';

fprintf('controller K = [%.3f  %.3f]\n', K(1), K(2));
fprintf('controller poles:\n  %.6f%+.6fj\n  %.6f%+.6fj\n', real(p_cl(1)), imag(p_cl(1)), real(p_cl(2)), imag(p_cl(2)));
fprintf('observer L = [%.3f; %.3f]\n', L(1), L(2));
p_ob = eig(Ao - L*Co);
fprintf('observer poles:\n  %.6f%+.6fj\n  %.6f%+.6fj\n\n', real(p_ob(1)), imag(p_ob(1)), real(p_ob(2)), imag(p_ob(2)));

% simulation setup
t_sim = 5.0;
dt    = 0.001;
t     = 0:dt:t_sim;
N     = numel(t);

% saturation
u_min = 0.0;
u_max = 1200.0;

% initial conditions
Sx0     = 0;
Vx0     = 19.5;
lambda0 = 0.30;

x_true   = zeros(3, N);
x_true(:,1) = [Sx0; Vx0; lambda0];

% observer state in deviation coordinates
xhat_dev = zeros(2, N);
xhat_dev(:,1) = [0; 0];

u = zeros(1, N);

% main simulation loop
for k = 1:N-1
    tk = t(k);

    Vx     = x_true(2,k);
    lambda = x_true(3,k);

    % choose friction params (disturbance)
    if tk < t_switch
        c_now = c_dry;
    else
        c_now = c_wet;
    end

    % nonlinear friction
    mu = mu_friction(lambda, Vx, c_now);

    % control law using estimated deviation state
    du = -K * xhat_dev(:,k);
    u_cmd = u_star + du;

    % saturate physical torque
    u_cmd = min(max(u_cmd, u_min), u_max);
    u(k)  = u_cmd;

    % nonlinear plant dynamics
    Vx_safe = max(Vx, 0.5);
    
    Sx_dot = Vx;
    Vx_dot = -mu * g;

    w = (1 - lambda) * Vx / R;
    w_dot = (mu*FN*R - u_cmd) / Jw;
    lambda_dot = -(R/Vx_safe)*w_dot + (R*w/(Vx_safe^2))*Vx_dot;

    % integrate true states
    x_true(:,k+1) = x_true(:,k) + dt * [Sx_dot; Vx_dot; lambda_dot];

    % observer (deviation coordinates)
    y_dev = lambda - lambda_star;
    u_dev = u_cmd - u_star;

    xhat_dev_dot = Ao*xhat_dev(:,k) + Bo*u_dev + L*(y_dev - Co*xhat_dev(:,k));
    xhat_dev(:,k+1) = xhat_dev(:,k) + dt*xhat_dev_dot;
end

u(end) = u(end-1);

% convert to absolute for plotting
lam_tr = x_true(3,:);
lamhat_tr = lambda_star + xhat_dev(2,:);

% performance metrics
tol = 0.02*abs(lambda_star);

settled_idx = find(abs(lam_tr - lambda_star) <= tol, 1, 'first');
if isempty(settled_idx)
    ts_actual = NaN;
else
    ts_actual = t(settled_idx);
end

overshoot = (max(lam_tr) - lambda_star) / lambda_star * 100;
ss_err = abs(lam_tr(end) - lambda_star);

fprintf('performance (4.5 disturbed nonlinear):\n');
fprintf('  overshoot = %.2f %%\n', overshoot);
fprintf('  ss error  = %.6f\n', ss_err);
if isnan(ts_actual)
    fprintf('  ts(2%%)    = NaN\n');
else
    fprintf('  ts(2%%)    = %.3f s\n', ts_actual);
end
fprintf('control effort:\n');
fprintf('  u range = [%.2f, %.2f] N*m\n\n', min(u), max(u));

% add comparison to section 4.4
fprintf('comparison to section 4.4 (no disturbance):\n');
fprintf('  4.4 overshoot: 50.17%%, 4.5 overshoot: %.2f%%\n', overshoot);
fprintf('  4.4 ts: 0.039s, 4.5 ts: ');
if isnan(ts_actual)
    fprintf('NaN\n');
else
    fprintf('%.3fs\n', ts_actual);
end
fprintf('  disturbance causes additional transient at t=%.2fs\n\n', t_switch);

% plotting results - clean 1x3 horizontal layout
figure('Position',[100 100 1400 400]);

% slip ratio tracking
subplot(1,3,1);
plot(t, lam_tr, 'LineWidth', 2); hold on;
plot(t, lamhat_tr, '--', 'LineWidth', 2);
yline(lambda_star, 'r--', 'LineWidth', 1.5);
yline(lambda_star+tol, 'g:', 'LineWidth', 1);
yline(lambda_star-tol, 'g:', 'LineWidth', 1);
xline(t_switch, 'color', [1 1 0], 'LineStyle', '--', 'LineWidth', 2);
grid on;
xlabel('Time (s)'); ylabel('lambda');
title('Slip ratio: true vs estimated');
legend('lambda (true)', 'lambdahat (estimated)', 'lambda* (equilibrium)', 'pm 2 percent band', 'friction drop', 'Location','best');

% control torque
subplot(1,3,2);
plot(t, u, 'LineWidth', 2); hold on;
yline(u_star, 'r--', 'LineWidth', 1.5);
yline(u_max, 'c:', 'LineWidth', 1.0);
yline(u_min, 'c:', 'LineWidth', 1.0);
xline(t_switch, 'color', [1 1 0], 'LineStyle', '--', 'LineWidth', 2);
grid on;
xlabel('Time (s)'); ylabel('u (Nm)');
title('Control input (with saturation)');
legend('u(t)', 'u*', 'limits', 'friction drop', 'Location','best');

% tracking error
subplot(1,3,3);
err = lam_tr - lambda_star;
plot(t, err, 'LineWidth', 2); hold on;
yline(0, 'r--', 'LineWidth', 1.5);
yline(tol, 'g:', 'LineWidth', 1);
yline(-tol, 'g:', 'LineWidth', 1);
xline(t_switch, 'color', [1 1 0], 'LineStyle', '--', 'LineWidth', 2);
grid on;
xlabel('Time (s)'); ylabel('error e = lambda - lambda*');
title('Tracking error');
legend('error','0','pm 2 percent band','friction drop','Location','best');

sgtitle(sprintf('Section 4.5: Disturbed CNL plant (dry->wet at t=%.2fs) with observer-based controller', t_switch), ...
    'FontSize', 14, 'FontWeight', 'bold');

%% local function
function mu = mu_friction(lambda, Vx, c)
    % keep lambda in reasonable range
    lam = max(min(lambda, 1.0), -0.2);
    
    mu = (c.c1*(1-exp(-c.c2*lam)) - c.c3*lam) * exp(-c.c4*Vx);
    
    % clamp mu to physical bounds
    mu = max(min(mu, 1.5), -1.5);
end