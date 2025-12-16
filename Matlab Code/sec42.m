%% reduced-order luenberger observer: convergence along a non-equilibrium trajectory
clear; clc; close all;

% vehicle and tire parameters
m  = 342;  Jw = 1.13;  Rr = 0.33;  g = 9.81;  FN = m*g;
c1 = 1.2801; c2 = 23.99; c3 = 0.52; c4 = 0.03;

% equilibrium values
lambda_star = 0.2;
Vx_star     = 20;
mu_star = (c1*(1-exp(-c2*lambda_star)) - c3*lambda_star) * exp(-c4*Vx_star);
u_star  = mu_star * FN * Rr * (1 + (1-lambda_star)*Jw/(m*Rr^2));

fprintf('=== reduced-order observer design ===\n');
fprintf('u* = %.3f N*m, mu* = %.4f\n\n', u_star, mu_star);

% reduced-order observable subsystem: x_o = [Vx; lambda], y = lambda
Ao = [0.1883  1.4362;
      0.3178  2.7380];
Bo = [0;
      0.0146];
Co = [0 1];

% controller design via pole placement
ts   = 2.0;
zeta = 0.9;
wn   = 4/(zeta*ts);
pc1 = -zeta*wn + 1i*wn*sqrt(1-zeta^2);
pc2 = -zeta*wn - 1i*wn*sqrt(1-zeta^2);
K   = place(Ao, Bo, [pc1 pc2]);
Acl = Ao - Bo*K;
p_cl = eig(Acl);

fprintf('controller poles:\n');
fprintf('  %.6f %+.6fj\n', real(p_cl(1)), imag(p_cl(1)));
fprintf('  %.6f %+.6fj\n', real(p_cl(2)), imag(p_cl(2)));
fprintf('K = [%.3f  %.3f]\n\n', K(1), K(2));

% observer design (faster than controller)
alpha = 5;
p_obs = alpha * p_cl;
L = place(Ao', Co', p_obs).';
Aobs = Ao - L*Co;
p_ob = eig(Aobs);

fprintf('observer poles (alpha=%.1f):\n', alpha);
fprintf('  %.6f %+.6fj\n', real(p_ob(1)), imag(p_ob(1)));
fprintf('  %.6f %+.6fj\n', real(p_ob(2)), imag(p_ob(2)));
fprintf('L = [%.3f; %.3f]\n\n', L(1), L(2));

% simulation setup
t_sim = 5.0;
dt    = 0.001;
t     = 0:dt:t_sim;
xo_star = [Vx_star; lambda_star];

% initial conditions
xo0 = [19.5; 0.30];      % true state (non-equilibrium)
xhat0 = [23.0; 0.10];    % initial estimate (intentionally wrong)

xo   = zeros(2, length(t));
xhat = zeros(2, length(t));
u    = zeros(1, length(t));
y    = zeros(1, length(t));
xo(:,1)   = xo0;
xhat(:,1) = xhat0;

for k = 1:length(t)-1
    % true deviation from equilibrium
    dx = xo(:,k) - xo_star;
    
    % control law
    du = -K * dx;
    u(k) = u_star + du;
    
    % measurement
    y(k) = Co * xo(:,k);
    
    % true system dynamics
    ddx = Acl * dx;
    xo(:,k+1) = xo(:,k) + ddx*dt;
    
    % observer update
    dxhat = xhat(:,k) - xo_star;
    dy    = y(k) - lambda_star;
    dyhat = (Co*xhat(:,k)) - lambda_star;
    ddxhat = Ao*dxhat + Bo*du + L*(dy - dyhat);
    xhat(:,k+1) = xhat(:,k) + ddxhat*dt;
end

u(end) = u(end-1);
y(end) = Co*xo(:,end);

% estimation error analysis
e = xo - xhat;
fprintf('final errors at t=%.1f s:\n', t_sim);
fprintf('  abs(Vx - Vxhat)      = %.8f\n', abs(e(1,end)));
fprintf('  abs(lambda - lambdahat) = %.8f\n\n', abs(e(2,end)));

% find convergence time (1% error band)
tol_convergence = 0.01 * Vx_star;
conv_idx = find(abs(e(1,:)) <= tol_convergence, 1, 'first');
if ~isempty(conv_idx)
    t_conv = t(conv_idx);
    fprintf('observer convergence time (1%% band): %.3f s\n\n', t_conv);
end

% plotting convergence
figure('Position', [120 120 1200 700]);

subplot(2,2,1);
plot(t, xo(1,:), 'LineWidth', 2); hold on;
plot(t, xhat(1,:), '--', 'LineWidth', 2);
yline(Vx_star, 'k:', 'LineWidth', 1.5);
if exist('t_conv','var')
    xline(t_conv, 'g--', 'LineWidth', 1, 'Alpha', 0.5);
end
grid on;
xlabel('Time (s)'); ylabel('Vx (m/s)');
title('Velocity: true vs estimated');
legend('Vx (true)','Vxhat (estimated)','Vx* (equilibrium)','Location','best');

subplot(2,2,2);
plot(t, xo(2,:), 'LineWidth', 2); hold on;
plot(t, xhat(2,:), '--', 'LineWidth', 2);
yline(lambda_star, 'k:', 'LineWidth', 1.5);
if exist('t_conv','var')
    xline(t_conv, 'g--', 'LineWidth', 1, 'Alpha', 0.5);
end
grid on;
xlabel('Time (s)'); ylabel('lambda');
title('Slip ratio: true vs estimated');
legend('lambda (true)','lambdahat (estimated)','lambda* (equilibrium)','Location','best');

subplot(2,2,3);
plot(t, e(1,:), 'LineWidth', 2); hold on;
yline(0, 'k--', 'LineWidth', 1.5);
if exist('t_conv','var')
    xline(t_conv, 'g--', 'LineWidth', 1, 'Alpha', 0.5);
end
grid on;
xlabel('Time (s)'); ylabel('Error (m/s)');
title('Estimation error in Vx');

subplot(2,2,4);
plot(t, e(2,:), 'LineWidth', 2); hold on;
yline(0, 'k--', 'LineWidth', 1.5);
if exist('t_conv','var')
    xline(t_conv, 'g--', 'LineWidth', 1, 'Alpha', 0.5);
end
grid on;
xlabel('Time (s)'); ylabel('Error');
title('Estimation error in lambda');

sgtitle('Section 4.2: Reduced-order Luenberger Observer - Convergence along a Non-equilibrium Trajectory', ...
'FontSize', 14, 'FontWeight', 'bold');