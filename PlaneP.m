clear; close all; clc; % 分界平面P--负斜率--fig5
[sigma, beta, rho, ~, Q, tspan, ~] = setParameters(); %初值与其他不同
%% 分界平面构造
Vpp = @(x1,x2,x3) x2+sqrt(beta*(rho-1))*x3-rho*sqrt(beta*(rho-1));  %分界平面 坐标范围(x1>0 & x2>0)
Vpn = @(x1,x2,x3) -x2+sqrt(beta*(rho-1))*x3-rho*sqrt(beta*(rho-1));  %分界平面 坐标范围(x1<0 & x2<0)
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

% 绘制分界面
n = 50;
xn = linspace(0,20,n); %画图范围
yn = linspace(0,30,n); 
zn = linspace(0,rho,n); 
[X, Y, Z] = meshgrid(xn, yn, zn);
vn = Vpp(X,Y,Z);
hn = patch(isosurface(X,Y,Z,vn,0));
isonormals(X,Y,Z,vn,hn)
set(hn,'FaceColor','y', 'EdgeColor','none', 'FaceAlpha',0.5);
xn = linspace(-20,0,n); %画图范围
yn = linspace(-30,0,n); 
[X, Y, Z] = meshgrid(xn, yn, zn);
vn = Vpn(X,Y,Z);
hn = patch(isosurface(X,Y,Z,vn,0));
isonormals(X,Y,Z,vn,hn)
set(hn,'FaceColor','y', 'EdgeColor','none', 'FaceAlpha',0.5);

light('Position',[1 -1 1],'Style','infinite');
lighting gouraud

% 绘制最终效果
paper_colors = [
    0,0.2,0.6;   % 1.深蓝
    0.3,0.7,0.9;   % 2.浅蓝
    0,0.4,0.2;   % 3.深绿
    0.5,0.9,0.5;  % 4.浅绿
    0.5,0,0.6;   % 5.深紫
    0.8,0.6,0.9];  % 6.浅紫

plot3(Q, Q, rho-1, 'k.'); %新平衡点Q+
plot3(-Q, -Q, rho-1, 'k.'); %新平衡点Q-
plot3(0,0,0,'k.'); %平衡点O
plot3(data(:,1), data(:,2), data(:,3),'r');
lps = [];
rps = [];
for ii = 2:length(data(:,1))
    % 判断分离点
    if data(ii,1)>0 && data(ii,2)>0 && Vpp(data(ii-1,1),data(ii-1,2),data(ii-1,3))<=0 && Vpp(data(ii,1),data(ii,2),data(ii,3))>=0 ...% 检测绕Q+向上
            || data(ii,1)>0 && data(ii,2)>0 && Vpp(data(ii-1,1),data(ii-1,2),data(ii-1,3))>=0 && Vpp(data(ii,1),data(ii,2),data(ii,3))<=0 % 检测绕Q+向下
        rps = [rps ii];
    elseif data(ii,1)<0 && data(ii,2)<0 && Vpn(data(ii-1,1),data(ii-1,2),data(ii-1,3))<=0 && Vpn(data(ii,1),data(ii,2),data(ii,3))>=0 ...% 检测绕Q-向上
            || data(ii,1)<0 && data(ii,2)<0 && Vpn(data(ii-1,1),data(ii-1,2),data(ii-1,3))>=0 && Vpn(data(ii,1),data(ii,2),data(ii,3))<=0 % 检测绕Q-向下
        lps = [lps ii];
    else
        continue
    end
end
% 绘制分离点
if ~isempty(lps)
    plot3(data(lps,1), data(lps,2), data(lps,3), '.', 'color',paper_colors(3,:));
end
if ~isempty(rps)
    plot3(data(rps,1), data(rps,2), data(rps,3), '.', 'color', paper_colors(4,:));
end

% 标注
text([0 Q -Q], [0 Q -Q], [0 rho+4 rho+4], {'$O$' '$Q_{+}$' '$Q_{-}$'},'interpreter','latex')
view(35,30)

annotation('textarrow',[0.2,0.2], [0.4,0.5],'String','$P^{-}$', 'FontSize', 24, 'interpreter','latex', 'HeadLength',10, 'HeadWidth',8)
annotation('textarrow',[0.7,0.7], [0.4,0.5],'String','$P^{+}$', 'FontSize', 24, 'interpreter','latex', 'HeadLength',10, 'HeadWidth',8)
zlim([0,50])

