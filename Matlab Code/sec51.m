%% section 5.1: open-loop steering of CLTI (minimum-energy) to the origin
clear; clc; close all;

% reduced-order CLTI system (observable subsystem)
Ao = [0.1883  1.4362;
      0.3178  2.7380];
Bo = [0;
      0.0146];

% operating point
Vx_star     = 20;
lambda_star = 0.2;
x_star      = [Vx_star; lambda_star];
u_star = 725.398;

% time horizon
T  = 2.0;
dt = 1e-3;
t  = 0:dt:T;
N  = numel(t);

% initial condition (absolute to deviation)
x0_abs = [21; 0.35];
x0     = x0_abs - x_star;
xf = [0; 0];

% finite-horizon controllability gramian
W = zeros(2,2);
for k = 1:N
    s  = t(k);
    E  = expm(Ao*s);
    W  = W + (E*Bo*Bo'*E')*dt;
end

Winv = pinv(W);

% minimum-energy open-loop input
delta = (xf - expm(Ao*T)*x0);
u_tilde = zeros(1,N);
for k = 1:N
    tau = t(k);
    u_tilde(k) = (Bo' * expm(Ao'*(T - tau)) * (Winv*delta));
end

u = u_star + u_tilde;

% simulate deviation dynamics
x_tilde = zeros(2,N);
x_tilde(:,1) = x0;
for k = 1:N-1
    xdot = Ao*x_tilde(:,k) + Bo*u_tilde(k);
    x_tilde(:,k+1) = x_tilde(:,k) + xdot*dt;
end

% convert to absolute states for plotting
x_abs = x_tilde + x_star;
Vx    = x_abs(1,:);
lam   = x_abs(2,:);

% terminal check
fprintf("=== section 5.1 (open-loop min-energy) ===\n");
fprintf("T = %.3f s\n", T);
fprintf("x0_abs = [%.4f; %.4f]\n", x0_abs(1), x0_abs(2));
fprintf("x_star = [%.4f; %.4f]\n", x_star(1), x_star(2));
fprintf("x0_tilde = x0_abs - x_star = [%.4f; %.4f]\n", x0(1), x0(2));
fprintf("final deviation x_tilde(T) = [%.6f; %.6f]\n", x_tilde(1,end), x_tilde(2,end));
fprintf("peak |u - u*| = %.3f N*m\n", max(abs(u_tilde)));
fprintf("u range = [%.2f, %.2f] N*m\n", min(u), max(u));

% plotting results
figure('Position', [100 100 900 700]);

subplot(3,1,1);
plot(t, Vx, 'LineWidth', 2); hold on;
yline(Vx_star,'--','LineWidth',1.5);
grid on;
ylabel('Vx (m/s)');
title('Section 5.1: Open-loop Steering of CLTI System (absolute states shown)');
legend('Vx(t)','Vx*','Location','best');

subplot(3,1,2);
plot(t, lam, 'LineWidth', 2); hold on;
yline(lambda_star,'--','LineWidth',1.5);
grid on;
ylabel('lambda');
legend('lambda(t)','lambda*','Location','best');

subplot(3,1,3);
plot(t, u, 'LineWidth', 2); hold on;
yline(u_star,'--','LineWidth',1.5);
grid on;
xlabel('Time (s)');
ylabel('u(t) (Nm)');
legend('u(t)','u*','Location','best');