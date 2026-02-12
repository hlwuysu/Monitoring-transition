%% 多方法比较---fig 11
clear; close all; clc;
dbstop if error
[sigma, beta, rho, ~, Q, tspan, sxyz] = setParameters();
%% 解
[t,sols] = ode45(@(t,d) eqLorenz(t, d, rho, beta, sigma), tspan, sxyz); % 数值求解
ntc = 15000; %去掉的暂态点数量
data = [sols(ntc:end,1), sols(ntc:end,2), sols(ntc:end,3)]; % 解的三维坐标
t = t(ntc:end);
clear sols
%% 
idRealTrans = find(data(1:end-1,1).*data(2:end,1)<=0); % x变号位置
tRealT = t(idRealTrans);
%% 本方法
% 分界曲面构造
kcn = @(x3) (sigma-1-sqrt((sigma-1)^2-4*sigma*(x3-rho)))./(2*sigma); %临界斜率-负
Vn = @(x1,x2,x3) x2-x1.*kcn(x3); %分界曲面表达式-负 坐标范围(x1>=0 & x3<rho)
Vpp = @(x1,x2,x3) x2+sqrt(beta*(rho-1))*(x3-rho); %分界平面表达式-正正 坐标范围(x1>=0 & x2>0)
Vpn = @(x1,x2,x3) -x2+sqrt(beta*(rho-1))*(x3-rho); %分界平面表达式-负负 坐标范围(x1<0 & x2<0)
kcp = @(x3) (sigma-1+sqrt((sigma-1)^2-4*sigma*(x3-rho)))./(2*sigma); %临界斜率-正==k'_c
Vp = @(x1,x2,x3) x2-x1.*kcp(0); %分界函数表达式-正 坐标范围(x1~=0 & 0<=x3<=rho+(s-1)^2/(4s))

% 计算监测点
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

% 提取迁移和非迁移段
ii = 1; ij = 1;
ms = {}; nms = {};
ims=1; inms=1;
for ij=1:length(mpid)
    if data(ii,1)*data(mpid(ij),1)>0
        nms{inms} = [data(ii:mpid(ij),:), t(ii:mpid(ij))]; % 非迁移段, 前三列坐标，第四列时间
        inms = inms+1;
    else
        ms{ims} = [data(ii:mpid(ij),:), t(ii:mpid(ij))]; % 迁移段, 前三列坐标，第四列时间
        ims = ims+1;
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
idm = [];
for ii=1:length(ms)
    tp1 = ms{ii};
    idm = [idm;tp1(1,:)];
end
% tthis = idm(:,4); % 点对应的时间
% sum(tRealT-tthis)/length(tRealT); % 平均提前时间
%%  2024_Ma
vlam = 2*Q.*sigma.*(data(:,2)-data(:,1))+2*Q.*(rho.*data(:,1)-data(:,2)-data(:,1).*data(:,3));
id24 = [];
for ii=2:length(vlam)-1
    if (vlam(ii)>0 && vlam(ii)<vlam(ii-1) && vlam(ii)<vlam(ii+1)) || (vlam(ii)<0 && vlam(ii)>vlam(ii-1) && vlam(ii)>vlam(ii+1))
        id24 = [id24,ii];
    end
end
% t24 = t(id24); % 点对应的时间
% sum(tRealT-t24)/length(tRealT); % 平均提前时间
%%  2021_Da
d1 = [data(2:end,:);zeros(1,3)];
d2 = [data(3:end,:);zeros(2,3)];
vs = d1(:,1).*d2(:,2)-data(:,1).*d2(:,2)+data(:,1).*d1(:,2)-d2(:,1).*d1(:,2)+d2(:,1).*data(:,2)+d1(:,1).*data(:,2);
id21 = [];
for ii=2:length(vs)-3
    if (vs(ii)>0 && vs(ii+1)<=0) || (vs(ii)<0 && vs(ii+1)>=0)
        id21 = [id21,ii];
    end
end
% tpt21 = t(id21); % 点对应的时间
% t21 = [tpt21(2),tpt21(4),tpt21(8),tpt21(10),tpt21(12), tpt21(16),tpt21(18),tpt21(20)...
%     tpt21(22),tpt21(24),tpt21(26),tpt21(28),tpt21(32),tpt21(34),tpt21(36),tpt21(38),tpt21(40)...
%     tpt21(44),tpt21(46),tpt21(48),tpt21(50),tpt21(52),tpt21(54),tpt21(56)];
% sum(tRealT-t21')/length(tRealT); % 平均提前时间

%% 分量时程图
ccs = {'$x$','$y$','$z$',};
htd = figure('color','w', 'Name','nnn', 'position',[100 100 1500 600]);
cid=1; % 坐标分量索引---1:x, 2:y, 3:z
addAxisName('$t$',ccs{cid});
hold on; box on;
plot(t,data(:,1), 'b');
hp1 = plot(idm(:,4), idm(:,cid), 'r*', 'markersize',10); %标记
hma = plot(t(id24),data(id24,cid),'k.', 'markersize',16); %Ma_2024标记
hda = plot(t(id21),data(id21,cid),'kx', 'markersize',10); %Da_2021标记

line(xlim, [0,0], 'Color','k', 'LineStyle','--'); %绘制x=0，y=0
line([28.61 28.79; 28.61 28.79], ylim, 'Color','k', 'LineStyle',':'); %绘制差值与Da21
line([idm(21,4) t(id24(21)); idm(21,4) t(id24(21))], ylim, 'Color','k', 'LineStyle',':'); %绘制差值与Ma24

legend([hp1,hma,hda],{'This method','Method in [32]','Method in [31]'}, 'interpreter','latex', 'Orientation','horizontal', 'FontSize',12)
xlim([t(1),t(end)])

% 设置纸张大小和位置（单位：英寸）
% set(gcf, ...
%     'PaperUnits', 'centimeter', ...       % 单位设为英寸
%     'PaperSize', [40 20]);   % 手动控制位置
% % 导出PDF
% print(gcf, 'exact_size.pdf', '-dpdf', '-r300');
annotation('arrow',[.2 .17],[.2 .2],'HeadWidth',8,'HeadLength',10)
annotation('textarrow',[.17 .2],[.2 .2],'String','0.18 ', 'interpreter','latex', 'TextEdgeColor','none','HeadWidth',8,'HeadLength',10,'FontSize',20)
annotation('arrow',[.2 .17],[.1 .1],'HeadWidth',8,'HeadLength',10)
annotation('textarrow',[.17 .2],[.1 .1],'String','0.43 ', 'interpreter','latex', 'TextEdgeColor','none','HeadWidth',8,'HeadLength',10,'FontSize',20)




