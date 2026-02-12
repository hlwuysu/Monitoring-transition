clear; close all; clc; % 更多算例---时程图---NEW---fig 13
[sigma, beta, ~, rhos, Q, tspan, sxyz] = setParameters(); % 初值【1 2 5】； %去掉的暂态点数量40000
hh = gobjects(size(rhos));  % 初始化句柄数组
nots = 20;  % 统一的转迁段数 100
for i1=1:length(rhos)
    rho = rhos(i1);
    ts = tspan;
    ms = {}; nms = {}; % 存放迁移和非迁移段
    ini = 0;
    while length(ms)<nots
        if ini~=0
            ts = [ts(1:end-1), ts(end)+tspan];
        end
        ini=ini+1;
        %% 分界曲面构造
        kcn = @(x3) (sigma-1-sqrt((sigma-1)^2-4*sigma*(x3-rho)))./(2*sigma); %临界斜率-负
        Vn = @(x1,x2,x3) x2-x1.*kcn(x3); %分界曲面表达式-负 坐标范围(x1>=0 & x3<rho)
        Vpp = @(x1,x2,x3) x2+sqrt(beta*(rho-1))*(x3-rho); %分界平面表达式-正正 坐标范围(x1>=0 & x2>0)
        Vpn = @(x1,x2,x3) -x2+sqrt(beta*(rho-1))*(x3-rho); %分界平面表达式-负负 坐标范围(x1<0 & x2<0)
        kcp = @(x3) (sigma-1+sqrt((sigma-1)^2-4*sigma*(x3-rho)))./(2*sigma); %临界斜率-正==k'_c
        Vp = @(x1,x2,x3) x2-x1.*kcp(0); %分界函数表达式-正 坐标范围(x1~=0 & 0<=x3<=rho+(s-1)^2/(4s))
        %% 解
        [~,sols] = ode45(@(t,s) eqLorenz(t,s, rho, beta, sigma), ts, sxyz);
        ntc = 20001; %去掉的暂态点数量VVV 20100
        data = [sols(ntc:end,1), sols(ntc:end,2), sols(ntc:end,3)]; % 数值解坐标
        t = ts(ntc:end)';
        clear sols
        coorspan = {[min( min(data(:,1)),min(-data(:,1))), max( max(data(:,1)),max(-data(:,1)))],...
        [min( min(data(:,2)),min(-data(:,2))), max( max(data(:,2)),max(-data(:,2)))]...
        [min( min(data(:,3)),min(-data(:,2))), max( max(data(:,3)),max(-data(:,2)))]}; %{[xmin,xmax], ...}

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
    end
    %% 全分量绘图
    figure('color','w', 'name', ['\rho= ',num2str(rho)], 'position',[100 100 1500 600]);
    ccs = {'$x$','$y$','$z$',};
    cid=1; % 坐标分量索引---1:x, 2:y, 3:z
    addAxisName('$t$',ccs{cid});
    hold on; box on;
    showTH(ms,cid,1) %绘制迁移段
    showTH(nms,cid,2) %绘制非迁移段
    xlim([t(1),t(end)])
    tp = coorspan{cid};
    ylim(tp+0.1*(tp(2)-tp(1))*[-1,1]);
    line(xlim, [0,0], 'Color','k', 'LineStyle','--'); %绘制x=0，y=0
    %% time difference绘图
    if i1==1
        htd = figure('color','w', 'Name','nnn'); % , 'position',[100 100 1500 600]
        xlabel('Number of transitions', 'interpreter','latex','FontSize',12)
%         ylabel('Time difference', 'interpreter','latex', 'rotation',90,'FontSize',12) % 对应函数showTD，纵轴为时间差
        ylabel('$x$-coordinate difference', 'interpreter','latex', 'rotation',90,'FontSize',12) % 对应函数showXD，纵轴为坐标差
        hold on; box on;
    end
%     [hh(i1),hg{i1},mvs{i1}] = showTD(ms,htd,rho,nots); %绘制
    [hh(i1),hg{i1},mvs{i1}] = showXD(ms,htd,rho,nots); %绘制
end

clrs = {'r','b','m','k', 'c','g','y'};
mkrs = {'+','o','*','^', 'x','pentagram','.','hexagram'};
gm = mvs{1};
for ii=1:length(rhos)
    set(hh(ii), 'Color',clrs{ii}, 'Marker',mkrs{ii})
    tpi = mvs{ii};
    gm = [min(gm(1),tpi(1)), max(gm(2),tpi(2))];
%     set(avg(ii), 'Color',clrs{ii})
end
hl1 = line(xlim,[gm(1), gm(1)], 'linestyle','--', 'color','k'); %最小值
hl2 = line(xlim,[gm(2), gm(2)], 'linestyle','-.', 'color','k'); %最大值
%返回当前图窗中的当前坐标区到ax1
ax1 = gca;
%ax2与ax1横纵坐标范围对应
ax2 = axes( 'Position',get(ax1,'Position'),'Visible','off');
ax3 = axes( 'Position',get(ax1,'Position'),'Visible','off');
%画两个legend
legd1 = legend(ax1,[hh(1),hh(2)],{hg{1},hg{2}}, 'interpreter','latex','FontSize',12, 'box','off'); %  'Orientation','horizontal'
legd2 = legend(ax2,[hh(3),hh(4)],{hg{3},hg{4}}, 'interpreter','latex','FontSize',12, 'box','off'); %  'Orientation','horizontal'
legd3 = legend(ax3,[hl2,hl1],{['Maximun:\ ',num2str(gm(2))],['Minimun:\ ',num2str(gm(1))]}, 'interpreter','latex','FontSize',12, 'box','off'); %  'Orientation','horizontal'

% legend([hh,hl2,hl1],[hg,['Maximun:\ ',num2str(gm(2))],['Minimun:\ ',num2str(gm(1))]], 'interpreter','latex','FontSize',12); %  'Orientation','horizontal'
%% 调用函数
function [hh,hlg,mvs] = showTD(spr,fid,rho,nots) % ---(元胞数组{x,y,z,t}， 图窗句柄，rho值,截止转迁段数)
% 至少两段迁移
% 返回图形句柄，图例名，最值
timdif = zeros(1,nots); %绘制迁移段初始点与X=0处的 时间差
for ii=1:nots
    tp = spr{ii};
    x = tp(:,1);
    t = tp(:,4);
        [~,id] = min(abs(x));
        timdif(ii) = t(id)-t(1);
end
figure(fid)
hold on
hh=plot(timdif, '.'); %绘制
mvs = [min(timdif), max(timdif)];
hlg = ['$\rho=\ $',num2str(rho)]; %图例
tp = gca;
tpx = tp.XLim; tpy = tp.YLim;
asx = [0, length(timdif)];
asy = [min(timdif),max(timdif)]+0.1*(max(timdif)-min(timdif)).*[-1,1];
xlim([0, min(tpx(2),asx(2))]); %X坐标范围
ylim([min(tpy(1),asy(1)), max(tpy(2),asy(2))]); %Y坐标范围
% avg = line(xlim,[sum(timdif)/asx(2), sum(timdif)/asx(2)], 'linestyle','--'); %平均值
end

function [hh,hlg,mvs] = showXD(spr,fid,rho,nots) % ---(元胞数组{x,y,z,t}， 图窗句柄，rho值,截止转迁段数)
% 至少两段迁移
% 返回图形句柄，图例名，最值
xdif = zeros(1,nots); %绘制迁移段初始点与X=0处的 X坐标差
for ii=1:nots
    tp = spr{ii};
% %     if tp(1,1)<tp(2,1) %X轨迹的走向由-向+
% %         if tp(1,1)<0
% %             xdif(ii) = -tp(1,1);
% %         end
% %     else
% %         xdif(ii) = tp(1,1);
% %     end
    xdif(ii) = abs(tp(1,1));
end
figure(fid)
hold on
hh=plot(xdif, '.'); %绘制
mvs = [min(xdif), max(xdif)];
hlg = ['$\rho=\ $',num2str(rho)]; %图例
tp = gca;
tpx = tp.XLim; tpy = tp.YLim;
asx = [0, length(xdif)];
asy = [min(xdif),max(xdif)]+0.1*(max(xdif)-min(xdif)).*[-1,1];
xlim([0, min(tpx(2),asx(2))]); %X坐标范围
ylim([min(tpy(1),asy(1)), max(tpy(2),asy(2))]); %Y坐标范围
% avg = line(xlim,[sum(timdif)/asx(2), sum(timdif)/asx(2)], 'linestyle','--'); %平均值
end


function showTH(spr,cid,ist) % ---(元胞数组{x,y,z,t}，坐标索引，标记[迁移为1，否则非1])
 %绘制范围内迁移段和非迁移段所有
if isempty(spr)
    if ist==1
        disp('no transitions')
    else
        disp('no non-transitions')
    end
else
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
end

function [cip,ipid] =  detectip(x, y, z) 
% 返回【拐点坐标，拐点索引】
tp1 = 1 +find(diff(sign(y), 1));
[cip,ipid] =  DIPS(x(tp1:end), y(tp1:end), z(tp1:end));
ipid = tp1+ipid-1; %在X中的索引
end

function [inflection_points,inflection_indices] =  DIPS(x, y, z)
% 输出：inflection_points-拐点坐标, inflection_indices-拐点的索引% 只返回找到的第一个
dx = x(2:end)-x(1:end-1); %第二个点开始
dy = y(2:end)-y(1:end-1);
kap = dy./dx;
Dkap = (kap(2:end)-kap(1:end-1))./(x(3:end)-x(2:end-1)); %倒数第二个点结束
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


