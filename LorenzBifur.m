% Lorenz系统分岔图VVV---fig 12
clear; clc; tic;
% Lorenz系统参数
sigma = 10;
beta = 8/3;
cid = 1;
% r参数范围（分叉参数）
r_start = 30;
r_end = 160;
r_step = 1e-2;
r_values = r_start:r_step:r_end;

% 初始化存储结果的变量
bps = [];
r_plot = [];
% 初始条件
x0 = [1,2,5];
% 积分时间设置
tspan = 0:0.001:100;  % 长时间积分以达到稳态
% options = odeset('RelTol',1e-6,'AbsTol',1e-6);
% 进度显示
progress = waitbar(0,'正在计算...');
for ii = 1:length(r_values)% 遍历每个r值
    r = r_values(ii);
    % 求解Lorenz系统
    [t, sol] = ode45(@(tt,xx) eqLorenz(tt,xx, r, beta, sigma), tspan, x0);

    % 只保留最后一部分数据（去除暂态）
    transient = 0.8;  % 去除前80%的暂态数据
    idx = round(transient*length(t));
    ts = tspan(idx:end);
    piks = sol(idx:end,:);  % 稳态点
    % 计算Poincare平面
%     Plane=[1,-1,0,0]; %Poincare平面 x-y=0平面(正方向)
    Plane=[0,0,1,1-r]; %Poincare平面 z=rho-1平面(正方向)
    [tsop,psop]=Solve_Poincare(ts,piks,Plane); %计算与平面交点
    if isempty(tsop) && abs(piks(end,1)-piks(end,2))<1e-3 %如果最后x和y足够接近，则认为收敛了
        tsop=ts(end);
        psop=piks(end,:);
    end
    % 存储结果
    bps = [bps; psop];
    r_plot = [r_plot; r*ones(length(tsop),1)];

    % 更新进度条
    waitbar(ii/length(r_values), progress, sprintf('计算进度: %.1f%%', ii/length(r_values)*100));
end

close(progress); toc;
cdn = {'$x$','$y$','$z$'};
% 绘制分叉图
figure('color','w')
hold on; box on;
scatter(r_plot, bps(:,cid), 1, 'b', 'filled'); % 
addAxisName('$\rho$',cdn{cid});
xlim([r_start r_end]);


% save('r_plot.mat', 'r_plot');  
% save('bpsonp.mat', 'bps');
% saveas(gcf, 'bifur_plot.fig'); 
%% 方程
function [t_List,P_List]=Solve_Poincare(t,sls,Plane)
% Plane 一般情况下是个垂直某个轴的平面
%一般只记录从负到正穿越。如果想反向也记录，可以设置Plane=-Plane。
%第二步，插值得到线与面的交点
P_List=[]; % 交点坐标
t_List=[]; % 交点时间
N=size(sls,1); %计算总共多少个点
tp1=[sls,ones(N,1)];
Dis=dot(tp1,ones(N,1)*Plane,2)./norm(Plane(1:end-1));
for k=1:N-1
    if Dis(k)<=0 && Dis(k+1)>0
        t0=t(k);t1=t(k+1);
        p0=sls(k,:);p1=sls(k+1,:);
        Dis0=Dis(k);Dis1=Dis(k+1);
        %一维线性插值，求Dis=0时的t和y
        %（相比较前面积分的4阶RK，这里用线性插值精度有点低）
        pp=p0+(p1-p0)/(Dis1-Dis0)*(0-Dis0);
        tp=t0+(t1-t0)/(Dis1-Dis0)*(0-Dis0);
        %储存
        P_List=[P_List;pp];
        t_List=[t_List,tp];
    end
end
end


