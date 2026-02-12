% 不同rho时的迁移时间统计---fig 14
clear; close all; clc; tic;
dbstop if error
% Lorenz系统参数
sigma = 10;
beta = 8/3;
cid = 1;
% r参数范围（分叉参数）
r_start = 20; % 20
r_end = 160;
r_step = 1; % 1e-1
r_values = r_start:r_step:r_end;

% 初始化存储结果的变量
bps = [];
r_plot = [];
Tall = [];
% 初始条件
x0 = [1,2,5];
% 积分时间设置
tspan = 0:0.001:50;  % 长时间积分以达到稳态
% options = odeset('RelTol',1e-6,'AbsTol',1e-6);
% 进度显示
progress = waitbar(0,'正在计算...');
for ii = 1:length(r_values)% 遍历每个r值
    rho = r_values(ii);
    % 分界面构造斜率-负S
    kcn = @(x3) (sigma-1-sqrt((sigma-1)^2-4*sigma*(x3-rho)))./(2*sigma); %临界斜率-负
    Vn = @(x1,x2,x3) x2-x1.*kcn(x3); %分界函数表达式-负 坐标范围(x1!=0 & any x3<rho)
    % 分界面构造斜率-正T
    kcp = @(x3) (sigma-1+sqrt((sigma-1)^2-4*sigma*(x3-rho)))./(2*sigma); %临界斜率-正==k'_c
    Vp = @(x1,x2,x3) x2-x1.*kcp(0); %分界函数表达式-正 坐标范围(x1~=0 & 0<=x3<=rho+(s-1)^2/(4s))
    % 分界平面构造-P
    Vpp = @(x1,x2,x3) x2+sqrt(beta*(rho-1))*x3-rho*sqrt(beta*(rho-1));  %分界平面 坐标范围(x1>0 & x2>0)
    Vpn = @(x1,x2,x3) -x2+sqrt(beta*(rho-1))*x3-rho*sqrt(beta*(rho-1));  %分界平面 坐标范围(x1<0 & x2<0)

    % 求解Lorenz系统
    [t, sol] = ode45(@(tt,xx) eqLorenz(tt,xx, rho, beta, sigma), tspan, x0);

    % 只保留最后一部分数据（去除暂态）
    transient = 0.8;  % 去除前80%的暂态数据
    idx = round(transient*length(t));
    ts = tspan(idx:end);
    piks = sol(idx:end,:);  % 稳态点

    %% 监测位于准迁移段的拐点以及解与T的交点
[lipids,ripids, LT,RT, LP,RP, LS,RS]  = detectips(piks(:,1), piks(:,2), piks(:,3), rho, beta, sigma);
idlip = find(lipids);
idrip = find(ripids);
if isempty(idlip) && isempty(idrip)
    bps = [bps; 0]; % 平均迁移时间
    Tall = [Tall; 0]; % 迁移时间，总
elseif isempty(idlip) && ~isempty(idrip)
    tp1 = ts(RT(idrip))-ts(ripids(idrip));
    bps = [bps; sum(tp1)/length(tp1)];
    Tall = [Tall; sum(tp1)];
elseif ~isempty(idlip) && isempty(idrip)
    tp1 = ts(LT(idlip))-ts(lipids(idlip));
    bps = [bps; sum(tp1)/length(tp1)];
    Tall = [Tall; sum(tp1)];
else
    tp = min(length(idlip),length(idrip));
    tp1 = ts(LT(idlip))-ts(lipids(idlip));
    tp2 = ts(RT(idrip))-ts(ripids(idrip));
    bps = [bps; (sum(tp1)+sum(tp2))/(length(tp1)+length(tp2))];
    Tall = [Tall; sum(tp1)+sum(tp2)];
end


%%  % 更新进度条
    waitbar(ii/length(r_values), progress, sprintf('计算进度: %.1f%%', ii/length(r_values)*100));
end

% % % save('MeanTime_vsrho.mat', 'bps');
% % % save('r_values4time.mat', 'r_values');
% % % save('TimeRatio_vsrho.mat', 'Tall');

close(progress); toc;

% 绘图
figure('color','w', 'name','Mean duration and time ratio of transitional state versus rho')
hold on; box on;
[hAx, hLine1, hLine2] = plotyy(r_values, bps, r_values, 100*Tall./(ts(end)-ts(1)));% plotyy(共用X轴, 左Y轴数据, 右Y轴数据, 绘图函数1, 绘图函数2)
% 设置左侧Y轴
ylabel(hAx(1), 'Mean duration transitional state (a.u.)', 'color','b', 'interpreter','latex', 'rotation',90);
set(hLine1, 'Color','b', 'LineStyle','-', 'LineWidth',1.5); % 'Marker', '.', 'LineStyle', 'none'
set(hAx(1), 'YColor', 'b');
% 设置右侧Y轴
ylabel(hAx(2), 'Time ratio of transitional state (\%)', 'color','r', 'interpreter','latex', 'rotation',90);
set(hLine2, 'Color','r', 'LineStyle','-', 'LineWidth',1.5); %
set(hAx(2), 'YColor','r')

% 公共设置
xlabel('$\rho$', 'interpreter','latex');
set(hAx, 'XLim',[r_start,r_end]);
set(hLine1, 'Color','b', 'LineStyle','-', 'LineWidth',1.5); %

%% 调用函数
function [lipids,ripids, LT,RT, LP,RP, LS,RS] =  detectips(x, y, z, r, b, s) 
% 返回【左拐点坐标索引，右拐点坐标索引，左终点索引，右终点索引；与左右P交点，与左右S交点】
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
% 斜率负-S
LS = []; RS = [];
for ii = 2:length(x)    % 判断分离点
    if z(ii)<r && x(ii)>0 && Vn(x(ii-1),y(ii-1),z(ii-1))<=0 && Vn(x(ii),y(ii),z(ii))>=0 % 负斜率检测从Q-到Q+
        RS = [RS ii];
    elseif z(ii)<r && x(ii)<0 && Vn(x(ii-1),y(ii-1),z(ii-1))>=0 && Vn(x(ii),y(ii),z(ii))<=0 % 负斜率检测从Q+到Q-
        LS = [LS ii];
    end
end
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
%                 continue;
                lip = [0,0,0];
                lipid = -(tp1-1);
            end
        else
            tp1 = RP(ij) +find(diff(sign(y(RP(ij):LT(ii)))), 1);
            [lip,lipid] =  DIPS(x(tp1:LT(ii)), y(tp1:LT(ii)), z(tp1:LT(ii)));
        end
        lips = [lips; lip];  %拐点坐标统计
        lipids = [lipids; lipid+tp1-1];  %拐点索引统计
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
%                 continue;
                rip = [0,0,0];
                ripid = -(tp1-1);
            end
        else
            tp1 = LP(ij) +find(diff(sign(y(LP(ij):RT(ii)))),1);
            [rip,ripid] =  DIPS(x(tp1:RT(ii)), y(tp1:RT(ii)), z(tp1:RT(ii)));
        end
        rips = [rips; rip];  %拐点坐标统计
        ripids = [ripids; ripid+tp1-1];  % 拐点索引统计
    end
end

end


function [inflection_points,inflection_indices] =  DIPS(x, y, z)
% 输出：inflection_points-拐点坐标, inflection_indices-拐点的索引% 只返回找到的第一个
% 计算一阶导数
dx = x(2:end)-x(1:end-1); %第二个点开始
dy = y(2:end)-y(1:end-1);
kap = dy./dx;
Dkap = (kap(2:end)-kap(1:end-1))./(x(3:end)-x(1:end-2));
sign_changes = diff(sign(Dkap)); % ~= 0
id = find(sign_changes);
% 返回拐点坐标
if isempty(id) %第二个开始无拐点，返回...第一个点 
    inflection_indices = 0; % [] 1
    inflection_points = [0,0,0];
else
    inflection_indices = id(1) + 1;
    inflection_points = [x(inflection_indices), y(inflection_indices), z(inflection_indices)];
end

end
