clear; close all; clc; % 一条Lorenz积分曲线--fig1+3+6+12
dbstop if error
[sigma, beta, rho, ~, Q, tspan, sxyz] = setParameters(); % 初始值和参数值
sxyz = [beta,sigma,rho-beta-sigma]; %初始位置
tspan=0:0.001:20;
[~,sols] = ode45(@(t,d) eqLorenz(t, d, rho, beta, sigma), tspan, sxyz); % 数值求解
ntc = 2000+1; %去掉的暂态点数量
data = [sols(ntc:end,1), sols(ntc:end,2), sols(ntc:end,3)]; % 解的三维坐标
clear sols

%% 三维解曲线
figure('color','w')
hold on; grid on;
plot3([Q -Q], [Q -Q], [rho-1 rho-1], 'b.'); % 平衡点Q±
plot3(0,0,0,'k.'); %平衡点O
plot3(data(:,1), data(:,2), data(:,3), 'r');
addAxisName('$x$','$y$','$z$')
text([0 Q -Q], [0 Q -Q], [0 rho-1 rho-1], {'$O$' '$Q_{+}$' '$Q_{-}$'},'interpreter','latex')
view(70,35)

%% 解曲线平面投影
figure('color','w')
% plot(sxyz(1),sxyz(2),'^'); %初始点
hold on; grid on;
plot([Q -Q], [Q -Q], 'b.'); % 平衡点Q±
plot(0,0,'k.'); %平衡点O
plot(data(:,1), data(:,2), 'r');
addAxisName('$x$','$y$')
text([Q -Q], [Q -Q],{'$Q_{+}$' '$Q_{-}$'},'interpreter','latex')
% 边界面y=kx投影

syms Z
n = 60;
xs1 = [linspace(-17,-0.2,n/2), linspace(0.2,17,n/2)]; % x坐标范围
op1 = double(subs(sigma-1-sqrt((1-sigma)^2-4*sigma*(Z-rho)), Z, rho/1.5)./(2*sigma)).*xs1; % 负斜率直线示意图

xs2 = [linspace(-13,-0.2,n/2), linspace(0.2,13,n/2)]; % x坐标范围
op2 = double(subs(sigma-1+sqrt((1-sigma)^2-4*sigma*(Z-rho)), Z, 0)./(2*sigma)).*xs2; % 正斜率直线示意图


