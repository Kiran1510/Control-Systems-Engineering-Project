%% sec53.m: ABS closed-loop with observer based control
clear; clc; close all;

% system matrices
A = [0.1883  1.4362;
     0.3178  2.7380];
B = [0;
     0.0146];
C = [0 1];

% equilibrium points
Vx_star = 20;
lambda_star = 0.2;
x_star = [Vx_star; lambda_star];
u_star = 725.398;

% controller and observer gains
K = [294.887 474.404];
L = [401.871; 22.926];

% input constraints
u_min = 0; u_max = 1200;

% simulation parameters
dt = 0.002;
Tmax = 5;
t = 0:dt:Tmax;
N = numel(t);

% initialize states
x = zeros(2,N);
xhat = zeros(2,N);
uabs = zeros(1,N);
x(:,1) = [19.5;0.30] - x_star;
xhat(:,1) = [2;-0.1];

% run closed-loop simulation
for k = 1:N-1
    y = C*x(:,k);
    u = u_star - K*xhat(:,k);
    u = min(max(u,u_min),u_max);
    uabs(k) = u;
    u_tilde = u - u_star;
    x(:,k+1) = x(:,k) + (A*x(:,k) + B*u_tilde)*dt;
    xhat(:,k+1) = xhat(:,k) + ...
        (A*xhat(:,k) + B*u_tilde + L*(y - C*xhat(:,k)))*dt;
end

% convert to absolute coordinates
x_abs = x + x_star;
xhat_abs = xhat + x_star;
Vx = x_abs(1,:);
lam = x_abs(2,:);
Vxhat = xhat_abs(1,:);
lamhat = xhat_abs(2,:);

% detect settling time and trim data
tol = 0.02;
settle_idx = find(abs(Vx - Vx_star) < tol*Vx_star,1,'first');
settle_idx = min(settle_idx + round(0.5/dt), N);
t = t(1:settle_idx);
Vx = Vx(1:settle_idx);
Vxhat = Vxhat(1:settle_idx);
lam = lam(1:settle_idx);
lamhat = lamhat(1:settle_idx);
uabs = uabs(1:settle_idx);
N = numel(t);

% setup figure
fig = figure('Color','k','Position',[100 100 1200 600]);
set(fig,'Renderer','opengl');

% left panel: car animation
axA = subplot(2,2,[1 3]);
hold on; grid on;
set(axA,'Color','k','XColor','w','YColor','w');
title(axA,'ABS Animation (Observer-based)','Color','w','FontWeight','bold');
axis(axA,'equal');
xlim(axA,[-6 26]);
ylim(axA,[-1.5 1.5]);
xlabel(axA,'Road position (m)','Color','w');
plot(axA,[-100 100],[0 0],'w','LineWidth',1.5);

% car dimensions
carL=2.2; carH=0.55; carY=0.25;
wheelR=0.22;

car = rectangle(axA,'Position',[0 carY carL carH],...
    'FaceColor',[0.2 0.55 1],'EdgeColor','none');
wheel = rectangle(axA,'Position',[0.45 0 2*wheelR 2*wheelR],...
    'Curvature',[1 1],'EdgeColor','w','LineWidth',2);

% info display at top
hud = annotation(fig,'textbox',[0.05 0.88 0.4 0.05],...
    'String','','Color','w','EdgeColor','w','FontSize',9,...
    'FontName','Courier','HorizontalAlignment','left',...
    'BackgroundColor',[0.1 0.1 0.1 0.8],'FitBoxToText','off');

% right top: velocity plot
axV = subplot(2,2,2); hold on; grid on;
set(axV,'Color','k','XColor','w','YColor','w');
title(axV,'Velocity: true vs estimated','Color','w');
ylabel(axV,'Vx (m/s)','Color','w');
xlim(axV,[0 t(end)]);
ylim(axV,[min([Vx Vxhat])-1, max([Vx Vxhat])+1]);
lineVx = plot(axV, nan, nan, 'b-', 'LineWidth', 2, 'DisplayName', 'Vx');
lineVxhat = plot(axV, nan, nan, 'r--', 'LineWidth', 2, 'DisplayName', 'Vx hat');
yline(axV,Vx_star,'--','Color',[0.7 0.7 0.7],'LineWidth',1,'DisplayName','Vx*');
legend(axV, 'TextColor', 'w', 'Location', 'best');

% right bottom: slip ratio plot
axL = subplot(2,2,4); hold on; grid on;
set(axL,'Color','k','XColor','w','YColor','w');
title(axL,'Slip ratio: true vs estimated','Color','w');
ylabel(axL,'lambda','Color','w');
xlabel(axL,'Time (s)','Color','w');
xlim(axL,[0 t(end)]);
ylim(axL,[min([lam lamhat])-0.2, max([lam lamhat])+0.2]);
lineLam = plot(axL, nan, nan, 'b-', 'LineWidth', 2, 'DisplayName', 'lambda');
lineLamhat = plot(axL, nan, nan, 'r--', 'LineWidth', 2, 'DisplayName', 'lambda hat');
yline(axL,lambda_star,'--','Color',[0.7 0.7 0.7],'LineWidth',1,'DisplayName','lambda*');
legend(axL, 'TextColor', 'w', 'Location', 'best');

% animation and export setup
gifName = 'sec53.gif';
aviName = 'sec53.avi';
vid = VideoWriter(aviName,'Motion JPEG AVI');
vid.FrameRate = 30;
open(vid);
skip = 10;
pos = -5;
scale_factor = 0.5;

for k = 1:skip:N
    % update car position based on velocity
    if k == 1
        pos = -5;
    else
        pos = -5 + scale_factor * trapz(t(1:k), Vx(1:k));
        pos = max(-5, min(pos, 23));
    end
    
    car.Position = [pos carY carL carH];
    wheel.Position = [pos+0.45 0 2*wheelR 2*wheelR];
    hud.String = sprintf('t=%.2fs | Vx=%.2f (est: %.2f) | lam=%.3f (est: %.3f) | u=%.0f Nm',...
        t(k),Vx(k),Vxhat(k),lam(k),lamhat(k),uabs(k));
    
    % update plots with current data
    set(lineVx, 'XData', t(1:k), 'YData', Vx(1:k));
    set(lineVxhat, 'XData', t(1:k), 'YData', Vxhat(1:k));
    set(lineLam, 'XData', t(1:k), 'YData', lam(1:k));
    set(lineLamhat, 'XData', t(1:k), 'YData', lamhat(1:k));
    
    drawnow;
    fr = getframe(fig);
    writeVideo(vid,fr);
    [im,cm] = rgb2ind(fr.cdata,256);
    if k==1
        imwrite(im,cm,gifName,'gif','LoopCount',inf,'DelayTime',1/30);
    else
        imwrite(im,cm,gifName,'gif','WriteMode','append','DelayTime',1/30);
    end
end

close(vid);
fprintf('\nSaved: %s and %s\n',gifName,aviName);