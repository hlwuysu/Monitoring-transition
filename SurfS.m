clear; close all; clc; % 分界面曲面S--fig4
[sigma, beta, rho, ~, Q, tspan, ~] = setParameters(); %初值与其他不同
%% 分界曲面构造
kcn = @(x3) (sigma-1-sqrt((sigma-1)^2-4*sigma*(x3-rho)))./(2*sigma); %临界斜率-负
Vn = @(x1,x2,x3) x2-x1.*kcn(x3); %分界函数表达式-负 坐标范围(x1>=0 & x3<rho)

%% 解
sxyz = [5 10 15];
[~,sols] = ode45(@(t,s) eqLorenz(t,s, rho, beta, sigma), tspan, sxyz);
ntc = 1; %去掉的暂态点数量
data = [sols(ntc:end,1), sols(ntc:end,2), sols(ntc:end,3)];
a=sxyz(1); b=sxyz(2); c=sxyz(3);
clear sols
%% 图窗
set(figure,'color','w', 'Name',['sigma: ',num2str(sigma), ' beta: ',num2str(beta), ' rho: ', num2str(rho)]);
addAxisName('$x$','$y$','$z$');
hold on; grid on;

%% 绘制分界曲面
n = 50;
xn = linspace(-20,-0.15,n); %画图范围x负
yn = linspace(-30,30,n); 
zn = linspace(0,rho,n); 
[X, Y, Z] = meshgrid(xn, yn, zn);
vn = Vn(X,Y,Z);
hn = patch(isosurface(X,Y,Z,vn,0));
isonormals(X,Y,Z,vn,hn)
set(hn,'FaceColor','y', 'EdgeColor','none', 'FaceAlpha',0.5);
xn = linspace(0.15,20,n); %画图范围x正
[X, Y, Z] = meshgrid(xn, yn, zn);
vn = Vn(X,Y,Z);
hn = patch(isosurface(X,Y,Z,vn,0));
isonormals(X,Y,Z,vn,hn)
set(hn,'FaceColor','y', 'EdgeColor','none', 'FaceAlpha',0.5);

light('Position',[1 -1 1],'Style','infinite');
lighting gouraud

%% 绘制最终效果
% 配色表
paper_colors = [
    0,0.2,0.6;   % 1.深蓝
    0.3,0.7,0.9;   % 2.浅蓝
    0,0.4,0.2;   % 3.深绿
    0.5,0.9,0.5;  % 4.浅绿
    0.5,0,0.6;   % 5.深紫
    0.8,0.6,0.9];  % 6.浅紫
% plot3(sxyz(1),sxyz(2),sxyz(3),'k'); %初始点+
plot3(Q, Q, rho-1, 'k.'); %新平衡点Q+
plot3(-Q, -Q, rho-1, 'k.'); %新平衡点Q-
plot3(0,0,0,'k.'); %平衡点O
plot3(data(:,1), data(:,2), data(:,3),'r');
lps = [];
rps = [];
for ii = 2:length(data(:,1))
    % 判断分离点
    if data(ii,3)<rho && data(ii,1)>0 && Vn(data(ii-1,1),data(ii-1,2),data(ii-1,3))<=0 && Vn(data(ii,1),data(ii,2),data(ii,3))>=0 % 负斜率检测仍绕Q+
        rps = [rps ii];
    elseif data(ii,3)<rho && data(ii,1)<0 && Vn(data(ii-1,1),data(ii-1,2),data(ii-1,3))>=0 && Vn(data(ii,1),data(ii,2),data(ii,3))<=0 % 负斜率检测仍绕Q-
        lps = [lps ii];
    end
end
% 绘制分离点
if ~isempty(lps)
    plot3(data(lps,1), data(lps,2), data(lps,3), '.', 'color',paper_colors(1,:));
end
if ~isempty(rps)
    plot3(data(rps,1), data(rps,2), data(rps,3), '.', 'color',paper_colors(2,:));
end

text([-1 Q -Q], [-5 Q -Q], [0 rho-1 rho-1], {'$O$' '$Q_{+}$' '$Q_{-}$'},'interpreter','latex')
view(35,30);
% ass = {'FontSize', 20, 'interpreter','latex', 'HeadLength',5, 'HeadWidth',5};
annotation('textarrow',[0.5,0.6], [0.2,0.3],'String','$S^{+}$', 'FontSize', 24, 'interpreter','latex', 'HeadLength',10, 'HeadWidth',8);
annotation('textarrow',[0.4,0.5], [0.2,0.3],'String','$S^{-}$', 'FontSize', 24, 'interpreter','latex', 'HeadLength',10, 'HeadWidth',8);
zlim([0,50])


