%% section 4.1: state feedback controller (reduced-order observable subsystem)
clear; clc; close all;

% vehicle and tire parameters
m = 342; Jw = 1.13; Rr = 0.33; g = 9.81; FN = m*g;
c1 = 1.2801; c2 = 23.99; c3 = 0.52; c4 = 0.03;

% equilibrium values
lambda_star = 0.2;
Vx_star = 20;
mu_star = (c1*(1-exp(-c2*lambda_star)) - c3*lambda_star)*exp(-c4*Vx_star);
u_star = mu_star*FN*Rr*(1 + (1-lambda_star)*Jw/(m*Rr^2));

% reduced-order observable subsystem (2 states: Vx, lambda)
Ao = [0.1883 1.4362; 0.3178 2.7380];
Bo = [0; 0.0146];

% design specifications
ts = 2.0;
zeta = 0.9;
wn = 4/(zeta*ts);
p1 = -zeta*wn + 1i*wn*sqrt(1-zeta^2);
p2 = -zeta*wn - 1i*wn*sqrt(1-zeta^2);
desired_poles = [p1 p2];

fprintf('=== section 4.1: reduced-order pole placement ===\n\n');
fprintf('design specs: ts=%.2fs, zeta=%.2f, wn=%.3f rad/s\n', ts, zeta, wn);
fprintf('desired poles: %.3f +/- %.3fj\n\n', real(p1), abs(imag(p1)));

% compute feedback gain
K = place(Ao, Bo, desired_poles);
Acl = Ao - Bo*K;
poles_cl = eig(Acl);

fprintf('feedback gain K = [%.3f, %.3f]\n', K(1), K(2));
fprintf('closed-loop poles: %.6f%+.6fj, %.6f%+.6fj\n\n', ...
    real(poles_cl(1)), imag(poles_cl(1)), real(poles_cl(2)), imag(poles_cl(2)));

% simulation setup
t_sim = 5.0; dt = 0.001; t = 0:dt:t_sim;
x0 = [0; 19.5; 0.30];
xo_star = [Vx_star; lambda_star];
x = zeros(3, length(t));
u = zeros(1, length(t));
x(:,1) = x0;

for k = 1:length(t)-1
    Vx = x(2,k);
    lam = x(3,k);
    xo = [Vx; lam];
    dxo = xo - xo_star;
    du = -K*dxo;
    u(k) = u_star + du;
    d_dxo = Acl*dxo;
    xdot = [Vx; d_dxo(1); d_dxo(2)];
    x(:,k+1) = x(:,k) + xdot*dt;
end
u(end) = u(end-1);

% performance metrics
lambda_traj = x(3,:);
tol = 0.02*abs(lambda_star);
settled_idx = find(abs(lambda_traj - lambda_star) <= tol, 1);
if isempty(settled_idx)
    ts_actual = NaN;
else
    ts_actual = t(settled_idx);
end
overshoot = (max(lambda_traj) - lambda_star)/lambda_star*100;
ss_error = abs(lambda_traj(end) - lambda_star);

fprintf('performance:\n');
if isnan(ts_actual)
    fprintf('  settling time: >5s\n');
else
    fprintf('  settling time: %.3fs\n', ts_actual);
end
fprintf('  overshoot: %.2f%%\n', overshoot);
fprintf('  ss error: %.6f\n\n', ss_error);

fprintf('control torque:\n');
fprintf('  u* = %.2f Nm\n', u_star);
fprintf('  max |u-u*| = %.2f Nm\n', max(abs(u-u_star)));
fprintf('  range: [%.2f, %.2f] Nm\n\n', min(u), max(u));

% plotting results
figure('Position', [100 100 1200 800]);

subplot(3,2,1);
plot(t, x(1,:), 'linewidth', 2); grid on;
xlabel('time (s)'); ylabel('Sx (m)'); title('position');

subplot(3,2,2);
plot(t, x(2,:), 'linewidth', 2); hold on;
yline(Vx_star, 'r--', 'linewidth', 1.5); grid on;
xlabel('time (s)'); ylabel('Vx (m/s)'); title('velocity');
legend('Vx', 'equilibrium');

subplot(3,2,3);
plot(t, x(3,:), 'linewidth', 2); hold on;
yline(lambda_star, 'r--', 'linewidth', 1.5);
yline(lambda_star+tol, 'g:'); yline(lambda_star-tol, 'g:'); grid on;
xlabel('time (s)'); ylabel('lambda'); title('slip ratio');
legend('lambda', 'reference', 'pm 2 percent');

subplot(3,2,4);
plot(t, u, 'linewidth', 2); hold on;
yline(u_star, 'r--', 'linewidth', 1.5); grid on;
xlabel('time (s)'); ylabel('u (Nm)'); title('control torque');
legend('u', 'equilibrium');

subplot(3,2,5);
plot(x(3,:), x(2,:), 'linewidth', 2); hold on;
plot(lambda_star, Vx_star, 'ro', 'markersize', 10, 'markerfacecolor', 'r');
plot(x(3,1), x(2,1), 'go', 'markersize', 8, 'markerfacecolor', 'g'); grid on;
xlabel('lambda'); ylabel('Vx (m/s)'); title('phase portrait');
legend('trajectory', 'equilibrium', 'initial');

subplot(3,2,6);
plot(t, x(3,:)-lambda_star, 'linewidth', 2); hold on;
yline(0, 'r--'); yline(tol, 'g:'); yline(-tol, 'g:'); grid on;
xlabel('time (s)'); ylabel('error'); title('tracking error');
legend('error', '0', 'pm 2 percent');

sgtitle('Section 4.1: Reduced-order Pole Placement Controller');