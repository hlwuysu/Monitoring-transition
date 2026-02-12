clear; close all; clc; % 分界面监测结果---时程图---NEW---fig 10
[sigma, beta, rho, ~, Q, ~, ~] = setParameters(); 
tspan = 0:0.001:50; 
sxyz = [beta,sigma,rho-beta-sigma];
%% 分界曲面构造
kcn = @(x3) (sigma-1-sqrt((sigma-1)^2-4*sigma*(x3-rho)))./(2*sigma); %临界斜率-负
Vn = @(x1,x2,x3) x2-x1.*kcn(x3); %分界曲面表达式-负 坐标范围(x1>=0 & x3<rho)
Vpp = @(x1,x2,x3) x2+sqrt(beta*(rho-1))*(x3-rho); %分界平面表达式-正正 坐标范围(x1>=0 & x2>0)
Vpn = @(x1,x2,x3) -x2+sqrt(beta*(rho-1))*(x3-rho); %分界平面表达式-负负 坐标范围(x1<0 & x2<0)
kcp = @(x3) (sigma-1+sqrt((sigma-1)^2-4*sigma*(x3-rho)))./(2*sigma); %临界斜率-正==k'_c
Vp = @(x1,x2,x3) x2-x1.*kcp(0); %分界函数表达式-正 坐标范围(x1~=0 & 0<=x3<=rho+(s-1)^2/(4s))
%% 解
[~,sols] = ode45(@(t,s) eqLorenz(t,s, rho, beta, sigma), tspan, sxyz);
ntc = 15000; %去掉的暂态点数量
data = [sols(ntc:end,1), sols(ntc:end,2), sols(ntc:end,3)]; % 数值解坐标
t = tspan(ntc:end)';
clear sols

%% 计算监测点
mpid = []; %统计标记点在data中的索引
for ii = 2:length(data(:,1))
    if (data(ii,1)>0 && data(ii,3)<rho && Vn(data(ii-1,1),data(ii-1,2),data(ii-1,3))<=0 && Vn(data(ii,1),data(ii,2),data(ii,3))>=0) ... %负斜率检测Q+附近下穿M
       || (data(ii,1)<0 && data(ii,3)<rho && Vn(data(ii-1,1),data(ii-1,2),data(ii-1,3))>=0 && Vn(data(ii,1),data(ii,2),data(ii,3))<=0) %负斜率检测Q-附近下穿M
        mpid = [mpid; ii-1];
    elseif Vp(data(ii-1,1),data(ii-1,2),data(ii-1,3))<=0 && Vp(data(ii,1),data(ii,2),data(ii,3))>=0 ...% 正斜率检测从Q+到Q-;  data(ii,3)<=rho+(sigma-1)^2/(4*sigma) && 
            || Vp(data(ii-1,1),data(ii-1,2),data(ii-1,3))>=0 && Vp(data(ii,1),data(ii,2),data(ii,3))<=0 % 正斜率检测从Q-到Q+;  data(ii,3)<=rho+(sigma-1)^2/(4*sigma) && 
        mpid = [mpid; ii-1];
    elseif data(ii,1)>0 && data(ii,2)>0 && data(ii,2)<=-double(Vp(data(ii,1),0,0)) && (Vpp(data(ii-1,1),data(ii-1,2),data(ii-1,3))<=0 && Vpp(data(ii,1),data(ii,2),data(ii,3))>=0) ...% 检测绕Q+向上
           || data(ii,1)<0 && data(ii,2)<0 && data(ii,2)>=-double(Vp(data(ii,1),0,0)) && (Vpn(data(ii-1,1),data(ii-1,2),data(ii-1,3))<=0 && Vpn(data(ii,1),data(ii,2),data(ii,3))>=0) ...% 检测绕Q-向上
        mpid = [mpid; ii];
    elseif data(ii,1)>0 && data(ii,2)>0 && data(ii,2)<=-double(Vp(data(ii,1),0,0)) && (Vpp(data(ii-1,1),data(ii-1,2),data(ii-1,3))>=0 && Vpp(data(ii,1),data(ii,2),data(ii,3))<=0) ... %检测Q+附近下穿平面
           || data(ii,1)<0 && data(ii,2)<0 && data(ii,2)>=-double(Vp(data(ii,1),0,0)) && (Vpn(data(ii-1,1),data(ii-1,2),data(ii-1,3))>=0 && Vpn(data(ii,1),data(ii,2),data(ii,3))<=0) %检测Q-附近下穿平面
        mpid = [mpid; ii-1];
    else
    end
end

%% 提取迁移和非迁移段
ii = 1; ij = 1;
ms = {}; nms = {};
ims=1; inms=1;
for ij=1:length(mpid)
    if data(ii,1)*data(mpid(ij),1)>0
        nms{inms} = [data(ii:mpid(ij),:), t(ii:mpid(ij))]; % 非迁移段, 前三列坐标，第四列时间
        inms = inms+1;
    else
        [ipc,ipid]  = detectip(data(ii:mpid(ij),1), data(ii:mpid(ij),2), data(ii:mpid(ij),3));
        if isempty(ipc) % 非迁移段
            nms{inms} = [data(ii:mpid(ij),:), t(ii:mpid(ij))];
            inms = inms+1;
        else
            ipid=ii+ipid;
            nms{inms} = [data(ii:ipid,:), t(ii:ipid)]; % 非迁移段
            inms = inms+1;
            ms{ims} = [data(ipid:mpid(ij),:), t(ipid:mpid(ij))]; % 迁移段
            ims = ims+1;
        end
    end
    ii = mpid(ij);
    if ij==length(mpid)
        if data(ii,1)*data(ii+1,1)<0
            ms{ims-1} = [ms{ims-1}; data(ii:end,:), t(ii:end)]; % 直到时间结束为迁移段
        else
            nms{inms-1} = [nms{inms-1}; data(ii:end,:), t(ii:end)]; % 直到时间结束为非迁移段
        end
    else
        continue;
    end
end

%% 全分量绘图
hh = figure('color','w', 'Name','m and nm');
ccs = {'$x$','$y$','$z$',};
axe = {};
for cid=1:3  % 坐标分量索引---1:x, 2:y, 3:z
    axe{cid} = subplot(3,1,cid);
    addAxisName('$t$',ccs{cid});
    hold on; box on;
    showTH(ms,cid,1) %绘制迁移段
    showTH(nms,cid,2) %绘制非迁移段
    if cid==1
        ys=get(gca,'ylim');
%         xeg1 = [16, -20, 8, 40];
        xeg2 = [25, -5, 10, 10]; % 放大区域[横起点，纵起点，宽度，高度]==[30, -20, 10, 40]
%         rectangle('Position',xeg1, 'EdgeColor','g'); %绘制放大区域
        rectangle('Position',xeg2, 'EdgeColor','g'); %语法：rectangle('Position', [x, y, width, height], 'PropertyName', PropertyValue)
    elseif cid==2
        ylim([-30,30]);
    end
    if cid~=3
        set(axe{cid}, 'XTickLabel', ''); % 隐藏上方子图的横轴刻度标签
        set(axe{cid}.XLabel, 'Visible','off'); % 隐藏上方子图的横轴标签
        line(xlim, [0,0], 'Color','k', 'LineStyle','--'); %绘制x=0，y=0
    else
%         zeg = [t(1), 23, t(end)-t(1), 5]; % 放大区域[横起点，纵起点，宽度，高度]
%         rectangle('Position',zeg, 'EdgeColor','g'); %绘制放大区域
    end
%     line([35.22,35.81;35.22,35.81], ylim, 'Color','k', 'LineStyle',':'); %绘制时间间隔最大的对应点
end

annotation(hh, 'line', [0.13,0.16],[0.96,0.96],'Color', 'b', 'linewidth',2);
annotation(hh, 'textbox',[0.16, 0.89, 0.2, 0.1],'String','Transitional segments', 'FontSize',24, 'FitBoxToText','on', 'EdgeColor','none', 'FaceAlpha',0, 'interpreter','latex');
annotation(hh, 'line', [0.39,0.42],[0.96,0.96],'Color', 'r', 'linewidth',2);
annotation(hh, 'textbox',[0.42, 0.89, 0.2, 0.1],'String','Non-transitional spirals', 'FontSize',24, 'FitBoxToText','on', 'EdgeColor','none', 'FaceAlpha',0, 'interpreter','latex');
annotation(hh, 'textbox', [0.685,0.89,0.2,0.1],'String','$\times $', 'Color','b', 'FontSize',24, 'FitBoxToText','on', 'EdgeColor','none', 'FaceAlpha',0, 'interpreter','latex');
annotation(hh, 'textbox',[0.7, 0.89, 0.2, 0.1],'String','Transitional starts','FontSize',24, 'FitBoxToText','on', 'EdgeColor','none', 'FaceAlpha',0, 'interpreter','latex');

linkaxes([axe{:}], 'x');
xlim([t(1),t(end)])
%% 单分量详细绘图---x
set(figure,'color','w', 'Name','m and nm');
ccs = {'$x$','$y$','$z$',};
cid = 1;  % 坐标分量索引---1:x, 2:y, 3:z
addAxisName('$t$',ccs{cid}); hold on; grid off; box on;
vst = [xeg2(1),xeg2(1)+xeg2(3)]; %放大时段
[ess,ets] = seaAE(nms,vst);
iss = seaAI(ms,vst);
showAETH(nms(ess(1):ess(2)),cid,2,ets)
showTH(ms(iss(1):iss(2)),cid,1)
xlim(vst)
ylim([xeg2(2),xeg2(2)+xeg2(4)])
line(xlim, [0,0], 'Color','k', 'LineStyle','--'); %绘制x=0

%% 调用函数
function [cip,ipid] =  detectip(x, y, z) 
% 返回【拐点坐标，拐点索引】
tp1 = 1 +find(diff(sign(y), 1));
[cip,ipid] =  DIPS(x(tp1:end), y(tp1:end), z(tp1:end));
ipid = tp1+ipid-1; %在X中的索引
end

function [inflection_points,inflection_indices] =  DIPS(x, y, z)
% 输出：inflection_points-拐点坐标, inflection_indices-拐点的索引% 只返回找到的第一个
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


function ret=ctavez(spr,ist) % ---(元胞数组{x,y,z,t}， 标记[迁移为1，否则非1])
 % 迁移标记点的平均Z坐标
 ret = [];
    for ii=1:length(spr)
        tp = spr{ii};
        if ist==1 %是迁移段
            ret = [ret,tp(1,3)];
        end
    end
    ret=sum(ret)/length(ret);
end


function showTH(spr,cid,ist) % ---(元胞数组{x,y,z,t}，坐标索引，标记[迁移为1，否则非1])
 %绘制范围内迁移段和非迁移段所有
    for ii=1:length(spr)
        tp = spr{ii};
        if ist==1 %是迁移段
            plot(tp(:,4),tp(:,cid), 'b');
            plot(tp(1,4),tp(1,cid), 'bx', 'markersize',5);
        else %非迁移段
            plot(tp(:,4),tp(:,cid), 'r');
        end
    end
end

function showAETH(spr,cid,ist,ts) % ---(范围内的元胞数组{x,y,z,t}，坐标索引，标记[迁移为1，否则非1]，t索引)
 %绘制所需范围内，与端点相同性质的区段
    for ii=1:length(spr)
        tp = spr{ii};
        if ii==1 %包含前端点的区段
            if ist==1 %是迁移段
                plot(tp(ts(1):end,4),tp(ts(1):end,cid), 'b');
                plot(tp(1,4),tp(1,cid), 'bx', 'markersize',5);
            else %非迁移段
                plot(tp(ts(1):end,4),tp(ts(1):end,cid), 'r');
            end
        elseif ii==length(spr) %包含后端点的区段
            if ist==1 %是迁移段
                plot(tp(1:ts(2),4),tp(1:ts(2),cid), 'b');
                plot(tp(1,4),tp(1,cid), 'bx', 'markersize',5);
            else %非迁移段
                plot(tp(1:ts(2),4),tp(1:ts(2),cid), 'r');
            end
        else %其他区段
            if ist==1 %是迁移段
                plot(tp(:,4),tp(:,cid), 'b');
                plot(tp(1,4),tp(1,cid), 'bx', 'markersize',5);
            else %非迁移段
                plot(tp(:,4),tp(:,cid), 'r');
            end
        end
    end
end

function [ss,ts]=seaAE(spr,vlt) % ---(元胞数组{x,y,z,t}，t区间)
%寻找t区间所在区段，t两端需同在迁移或非迁移段，且区间内需包含与端点不同性质的段
%ss--[开始区段索引，终止区段索引]；ts--[t前端点在开始区段中的索引，t后端点在终止区段中的索引]
tpt=1;
for ii=1:length(spr)
    tp = spr{ii};
    ish = find(tp(:,4)==vlt(tpt));
    if isempty(ish) %不在这段
        continue;
    else %在这段
        ss(tpt) = ii;
        ts(tpt) = ish;
        if tpt==2
            break;
        else
            tpt=2;
            continue;
        end
    end
end
end

function ss=seaAI(spr,vlt) % ---(元胞数组{x,y,z,t}，t区间)
%寻找包含在t区内的与端点不同性质的段
% ss--[开始区段索引，终止区段索引]
tpt=1;
for ii=1:length(spr)
    tp = spr{ii};
    if tpt==1
        ish = find(tp(:,4)>=vlt(tpt),1);
    else
        ish = find(tp(:,4)<=vlt(tpt),1,'last');
    end
    if isempty(ish) %不在这段
        continue;
    else %在这段
        ss(tpt) = ii;
        if tpt==1
            tpt=2;
        else
            tp1 = tp(:,4);
            if tp1(1)>vlt(2)
                ss(tpt) = ii-1;
                break;
            else
                continue;
            end
        end
    end
end
end

