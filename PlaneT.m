clear; close all; clc; % 分界面T以及投影拐点监测---正斜率---fig7
[sigma, beta, rho, ~, Q, tspan, ~] = setParameters(); %初值与其他不同
%% 分界面构造
kcp = @(x3) (sigma-1+sqrt((sigma-1)^2-4*sigma*(x3-rho)))./(2*sigma); %临界斜率-正
Vp = @(x1,x2,x3) x2-x1.*kcp(0); %分界函数表达式-正 坐标范围(x1>=0 & any x3<rho)

%% 解
sxyz = [5 10 15];
[~,sols] = ode45(@(t,s) eqLorenz(t,s, rho, beta, sigma), tspan, sxyz);
ntc = 1; %去掉的暂态点数量
data = [sols(ntc:end,1), sols(ntc:end,2), sols(ntc:end,3)];
a=sxyz(1); b=sxyz(2); c=sxyz(3);
clear sols

%% 监测解与T的交点
lps = [];
rps = [];
for ii = 2:length(data(:,1))
    % 判断分离点
    if Vp(data(ii-1,1),data(ii-1,2),data(ii-1,3))<=0 && Vp(data(ii,1),data(ii,2),data(ii,3))>=0 % 正斜率检测从Q+到Q-
        rps = [rps ii];
    elseif Vp(data(ii-1,1),data(ii-1,2),data(ii-1,3))>=0 && Vp(data(ii,1),data(ii,2),data(ii,3))<=0 % 正斜率检测从Q-到Q+
        lps = [lps ii];
    end
end

%% 图窗
set(figure,'color','w', 'Name',['sigma: ',num2str(sigma), ' beta: ',num2str(beta), ' rho: ', num2str(rho)]);
addAxisName('$x$','$y$','$z$');
hold on; grid on;

% 绘制分界曲面
n = 50;
xn = linspace(-20,-0.1,n); %画图范围x负
yn = linspace(-30,30,n);
zn = linspace(0,rho,n);
[X, Y, Z] = meshgrid(xn, yn, zn);
vn = Vp(X,Y,Z);
hn = patch(isosurface(X,Y,Z,vn,0));
isonormals(X,Y,Z,vn,hn)
set(hn,'FaceColor','y', 'EdgeColor','none', 'FaceAlpha',0.5);
xn = linspace(0.1,20,n); %画图范围x正
[X, Y, Z] = meshgrid(xn, yn, zn);
vn = Vp(X,Y,Z);
hp = patch(isosurface(X,Y,Z,vn,0)); % 绘制面
isonormals(X,Y,Z,vn,hp)
set(hp,'FaceColor','y', 'EdgeColor','none', 'FaceAlpha',0.5);

light('Position',[1 -1 1],'Style','infinite');
lighting gouraud

% 绘制最终效果
% 论文同图10种高对比配色 - RGB(0-1)，对比鲜明、印刷/黑白打印友好
paper_colors = [
    0,0.2,0.6;   % 1.深蓝
    0.3,0.7,0.9;   % 2.浅蓝
    0,0.4,0.2;   % 3.深绿
    0.5,0.9,0.5;  % 4.浅绿
    0.5,0,0.6;   % 5.深紫
    0.8,0.6,0.9];  % 6.浅紫

% 监测位于准迁移段的拐点以及解与T的交点
[lips,rips, ~,~, ~,~]  = detectips(data(:,1), data(:,2), data(:,3), rho, beta, sigma);
if ~isempty(lips)
    plot3(lips(:,1), lips(:,2), lips(:,3), '*', 'color',paper_colors(6,:));
end
if ~isempty(rips)
    plot3(rips(:,1), rips(:,2), rips(:,3), '*', 'color',paper_colors(5,:));
end

plot3(Q, Q, rho-1, 'k.'); %新平衡点Q+
plot3(-Q, -Q, rho-1, 'k.'); %新平衡点Q-
plot3(0,0,0,'k.'); %平衡点O
plot3(data(:,1), data(:,2), data(:,3),'r');
% 绘制分离点
if ~isempty(lps)
    plot3(data(lps,1), data(lps,2), data(lps,3), '.', 'color',paper_colors(5,:));
end
if ~isempty(rps)
    plot3(data(rps,1), data(rps,2), data(rps,3), '.', 'color',paper_colors(6,:));
end

text([-1 Q -Q], [-5 Q -Q], [0 rho-1 rho-1], {'$O$' '$Q_{+}$' '$Q_{-}$'},'interpreter','latex')
view(35,30);
annotation('textarrow',[0.5,0.6], [0.2,0.3],'String','$T^{+}$', 'FontSize', 20, 'interpreter','latex', 'HeadLength',10, 'HeadWidth',8)
annotation('textarrow',[0.4,0.5], [0.2,0.3],'String','$T^{-}$', 'FontSize', 20, 'interpreter','latex', 'HeadLength',10, 'HeadWidth',8)
zlim([0,50])

%% 
function [lips,rips, LT,RT, LP,RP] =  detectips(x, y, z, r, b, s) 
% 返回【左拐点坐标，右拐点坐标，左终点索引，右终点索引】
% 分界面构造斜率-负-S
kcn = @(x3) (s-1-sqrt((s-1)^2-4*s*(x3-r)))./(2*s); %临界斜率-负
Vn = @(x1,x2,x3) x2-x1.*kcn(x3); %分界函数表达式-负 坐标范围(x1!=0 & any x3<rho)
% 分界面构造斜率-正-T
kcp = @(x3) (s-1+sqrt((s-1)^2-4*s*(x3-r)))./(2*s); %临界斜率-正==k'_c
Vp = @(x1,x2,x3) x2-x1.*kcp(0); %分界函数表达式-正 坐标范围(x1~=0 & 0<=x3<=rho+(s-1)^2/(4s))
% 分界平面构造-P
Vpp = @(x1,x2,x3) x2+sqrt(b*(r-1))*x3-r*sqrt(b*(r-1));  %分界平面 坐标范围(x1>0 & x2>0)
Vpn = @(x1,x2,x3) -x2+sqrt(b*(r-1))*x3-r*sqrt(b*(r-1));  %分界平面 坐标范围(x1<0 & x2<0)
% 监测点
% 斜率正-T
LT = []; RT = [];
for ii = 2:length(x)    % 判断分离点
    if Vp(x(ii-1),y(ii-1),z(ii-1))<=0 && Vp(x(ii),y(ii),z(ii))>=0 % 正斜率检测从Q+到Q-;  data(ii,3)<=rho+(sigma-1)^2/(4*sigma) && 
        LT = [LT ii];
    elseif Vp(x(ii-1),y(ii-1),z(ii-1))>=0 && Vp(x(ii),y(ii),z(ii))<=0 % 正斜率检测从Q-到Q+;  data(ii,3)<=rho+(sigma-1)^2/(4*sigma) && 
        RT = [RT ii];
    end
end
% 平面-P
LP = []; RP = [];
for ii = 2:length(x(:,1))
    % 判断分离点
    if x(ii)>0 && y(ii)>0 && y(ii)<=-double(Vp(x(ii),0,0)) && (Vpp(x(ii-1),y(ii-1),z(ii-1))<=0 && Vpp(x(ii),y(ii),z(ii))>=0 ...% 检测绕Q+向上
            || Vpp(x(ii-1),y(ii-1),z(ii-1))>=0 && Vpp(x(ii),y(ii),z(ii))<=0) % 检测绕Q+向下
        RP = [RP ii];
    elseif x(ii)<0 && y(ii)<0 && y(ii)>=-double(Vp(x(ii),0,0)) && (Vpn(x(ii-1),y(ii-1),z(ii-1))<=0 && Vpn(x(ii),y(ii),z(ii))>=0 ...% 检测绕Q-向上
            || Vpn(x(ii-1),y(ii-1),z(ii-1))>=0 && Vpn(x(ii),y(ii),z(ii))<=0) % 检测绕Q-向下
        LP = [LP ii];
    end
end
% 找去Q-拐点*********
lips = []; lipids = [];
if ~isempty(LT)
    for ii=1:length(LT)
        ij = find(RP<LT(ii), 1, 'last'); %找T交点前的与P的交点
        if isempty(ij)
            if ii==1
                if isempty(find(diff(sign(y(1:LT(ii)))), 1))
                    tp1 = 1;
                else
                    tp1 = 1 +find(diff(sign(y(1:LT(ii)))), 1);
                end
                [lip,lipid] =  DIPS(x(tp1:LT(ii)), y(tp1:LT(ii)), z(tp1:LT(ii)));
            else
                continue;
            end
        else
            tp1 = RP(ij) +find(diff(sign(y(RP(ij):LT(ii)))), 1);
            [lip,lipid] =  DIPS(x(tp1:LT(ii)), y(tp1:LT(ii)), z(tp1:LT(ii)));
        end
        lips = [lips; lip];  %拐点坐标统计
        lipids = [lipids; lipid];  %拐点索引统计
    end
end
% 找去Q+的拐点xxxxxxxxxxx
rips = []; ripids = [];
if ~isempty(RT)
    for ii=1:length(RT)
        ij = find(LP<RT(ii), 1,'last'); %找T交点前的与P的交点
        if isempty(ij)
            if ii==1
                if isempty(find(diff(sign(y(1:RT(ii)))),1))
                    tp1=1;
                else
                    tp1=1 +find(diff(sign(y(1:RT(ii)))),1);
                end
                [rip,ripid] =  DIPS(x(tp1:RT(ii)), y(tp1:RT(ii)), z(tp1:RT(ii)));
            else
                continue;
            end
        else
            tp1 = LP(ij) +find(diff(sign(y(LP(ij):RT(ii)))),1);
            [rip,ripid] =  DIPS(x(tp1:RT(ii)), y(tp1:RT(ii)), z(tp1:RT(ii)));
        end
        rips = [rips; rip];  %拐点坐标统计
        ripids = [ripids; ripid];  % 拐点索引统计
    end
end
end

function [inflection_points,inflection_indices] =  DIPS(x, y, z)
% 输出：inflection_points-拐点坐标, inflection_indices-拐点的索引% 只返回找到的第一个
% 计算一阶导数
dx = x(2:end)-x(1:end-1); %di2gekaishi 
dy = y(2:end)-y(1:end-1);
kap = dy./dx;
Dkap = (kap(2:end)-kap(1:end-1))./(x(3:end)-x(1:end-2));
sign_changes = diff(sign(Dkap)); % ~= 0
id = find(sign_changes);
% 返回拐点坐标
if isempty(id)
    inflection_indices = [];
    inflection_points = [];
else
    inflection_indices = id(1) + 1;
    inflection_points = [x(inflection_indices), y(inflection_indices), z(inflection_indices)];
end
end
