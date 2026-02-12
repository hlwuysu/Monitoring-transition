clear; close all; clc; % 分界面监测过程 --- S,P,T,Ips ---NEW---fig9
[sigma, beta, rho, ~, Q, ~, ~] = setParameters();
%% 分界曲面构造
kcn = @(x3) (sigma-1-sqrt((sigma-1)^2-4*sigma*(x3-rho)))./(2*sigma); %临界斜率-负
Vn = @(x1,x2,x3) x2-x1.*kcn(x3); %分界曲面表达式-负 坐标范围(x1>=0 & x3<rho)
Vpp = @(x1,x2,x3) x2+sqrt(beta*(rho-1))*x3-rho*sqrt(beta*(rho-1)); %分界平面表达式-正正 坐标范围(x1>=0 & x2>0)
Vpn = @(x1,x2,x3) -x2+sqrt(beta*(rho-1))*x3-rho*sqrt(beta*(rho-1)); %分界平面表达式-负负 坐标范围(x1<0 & x2<0)
kcp = @(x3) (sigma-1+sqrt((sigma-1)^2-4*sigma*(x3-rho)))./(2*sigma); %临界斜率-正==k'_c
Vp = @(x1,x2,x3) x2-x1.*kcp(0); %分界函数表达式-正 坐标范围(x1~=0 & 0<=x3<=rho+(s-1)^2/(4s))
sxyz = [5,10,15];
tspan = 0:0.001:20;
%% 解
[~,sols] = ode45(@(t,s) eqLorenz(t,s, rho, beta, sigma), tspan, sxyz);
ntc = 1; %去掉的暂态点数量，不去暂态点1
data = [sols(ntc:end,1), sols(ntc:end,2), sols(ntc:end,3)];  % 数值解坐标
clear sols

%% 监测点
mps = []; %统计标记点在data中的索引
for ii = 2:length(data(:,1))
    if (data(ii,1)>0 && data(ii,3)<rho && Vn(data(ii-1,1),data(ii-1,2),data(ii-1,3))<=0 && Vn(data(ii,1),data(ii,2),data(ii,3))>=0) ... %检测Q+附近下穿S
       || (data(ii,1)<0 && data(ii,3)<rho && Vn(data(ii-1,1),data(ii-1,2),data(ii-1,3))>=0 && Vn(data(ii,1),data(ii,2),data(ii,3))<=0) %检测Q-附近下穿S
        mps = [mps; ii-1];
    elseif data(ii,3)<=rho+(sigma-1)^2/(4*sigma) && Vp(data(ii-1,1),data(ii-1,2),data(ii-1,3))<=0 && Vp(data(ii,1),data(ii,2),data(ii,3))>=0 ...% 正斜率检测从Q+到Q-
            || data(ii,3)<=rho+(sigma-1)^2/(4*sigma) && Vp(data(ii-1,1),data(ii-1,2),data(ii-1,3))>=0 && Vp(data(ii,1),data(ii,2),data(ii,3))<=0 % 正斜率检测从Q-到Q+
        mps = [mps; ii-1];
    elseif data(ii,1)>0 && data(ii,2)>0 && data(ii,2)<=-double(Vp(data(ii,1),0,0)) && (Vpp(data(ii-1,1),data(ii-1,2),data(ii-1,3))<=0 && Vpp(data(ii,1),data(ii,2),data(ii,3))>=0) ...% 检测绕Q+向上
           || data(ii,1)<0 && data(ii,2)<0 && data(ii,2)>=-double(Vp(data(ii,1),0,0)) && (Vpn(data(ii-1,1),data(ii-1,2),data(ii-1,3))<=0 && Vpn(data(ii,1),data(ii,2),data(ii,3))>=0) ...% 检测绕Q-向上
        mps = [mps; ii];
    elseif data(ii,1)>0 && data(ii,2)>0 && data(ii,2)<=-double(Vp(data(ii,1),0,0)) && (Vpp(data(ii-1,1),data(ii-1,2),data(ii-1,3))>=0 && Vpp(data(ii,1),data(ii,2),data(ii,3))<=0) ... %检测Q+附近下穿P
           || data(ii,1)<0 && data(ii,2)<0 && data(ii,2)>=-double(Vp(data(ii,1),0,0)) && (Vpn(data(ii-1,1),data(ii-1,2),data(ii-1,3))>=0 && Vpn(data(ii,1),data(ii,2),data(ii,3))<=0) %检测Q-附近下穿P
%         if isempty(mps)
%         elseif data(mps(end),1)*data(ii-1,1)<0
%         else
            mps = [mps; ii-1];
%         end
    else
    end
end
%% 创建视频写入对象（指定格式、帧率、分辨率）
videoName = 'test01.mp4';  % 输出文件名
vid = VideoWriter(videoName, 'MPEG-4');  % 格式：MP4（AVI可改'Motion JPEG AVI'）
vid.FrameRate = 24;  % 帧率（每秒20帧，常用15-30）
vid.Quality = 100;   % 视频质量（0-100）
open(vid);  % 打开视频写入流
ptim = 0.5;

%% 图窗
% % set(figure,'color','w', 'Name',['sigma: ',num2str(sigma), ' beta: ',num2str(beta), ' rho: ', num2str(rho)]);
fig = figure('color','w');
addAxisName('$x$','$y$','$z$');
hold on; grid on;
%% 绘制平衡点
plot3(sxyz(1),sxyz(2),sxyz(3),'k+'); %初始点+
%% 绘制分界面
n = 50;
gp = 0.15; %用于显示x不等于0的空白
xn = linspace(-20,-gp,n); %画图范围负斜率---S-
yn = linspace(-30,30,n); 
zn = linspace(0,rho,n); 
[X, Y, Z] = meshgrid(xn, yn, zn);
vn = Vn(X,Y,Z);
hn = patch(isosurface(X,Y,Z,vn,0));
isonormals(X,Y,Z,vn,hn)
set(hn,'FaceColor','y', 'EdgeColor','none', 'FaceAlpha',0.5);
xn = linspace(gp,20,n); %画图范围负斜率---S+
[X, Y, Z] = meshgrid(xn, yn, zn);
vn = Vn(X,Y,Z);
hn = patch(isosurface(X,Y,Z,vn,0));
isonormals(X,Y,Z,vn,hn)
set(hn,'FaceColor','y', 'EdgeColor','none', 'FaceAlpha',0.5);
% *******************
clear X Y Z
xn = linspace(gp,20,n); %画图范围平面---P+
yn = -double(Vp(xn,0,0)); %note sign -
for ii=1:n
    Y(:,ii) = linspace(0,yn(ii),n);
end
X = repmat(xn,n,1);
Z = -Y.*sign(X)./Q+rho;
surf(X,Y,Z, 'EdgeColor','none','FaceColor','y', 'FaceAlpha',0.5)
xn = linspace(-20,-gp,n); %画图范围平面---P-
yn =  -double(Vp(xn,0,0)); %note sign -
for ii=1:n
    Y(:,ii) = linspace(0,yn(ii),n);
end
X = repmat(xn,n,1);
Z = -Y.*sign(X)./Q+rho;
surf(X,Y,Z, 'EdgeColor','none','FaceColor','y', 'FaceAlpha',0.5)
% *******************
clear X Y Z
xn = linspace(-20,-gp,n); %画图范围正斜率---面T-
yn = linspace(-Vp(min(xn),0,0),0,n);
zn = zeros(1,n);
for ii=1:n
    zn(ii) = -sign(xn(1))*yn(ii)./Q+rho;
end
for ii=1:n
    Z(:,ii) = linspace(0,zn(ii),n); % note descending order
end
X = repmat(xn,n,1);
Y = ((sigma-1).*X+X.*sqrt((sigma-1).^2+4*sigma.*rho)) ./(2*sigma);
surf(X,Y,Z, 'EdgeColor','none','FaceColor','y', 'FaceAlpha',0.5)
xn = linspace(gp,20,n); %画图范围正斜率---面T+
yn = linspace(-Vp(max(xn),0,0),0,n);
zn = zeros(1,n);
for ii=1:n
    zn(ii) = -sign(xn(1))*yn(ii)./Q+rho;
end
% *************************
clear X Y Z
for ii=1:n
    Z(:,ii) = linspace(0,zn(n+(ii-1)*(n-1)/(1-n)),n); % note descending order
end
X = repmat(xn,n,1);
Y = ((sigma-1).*X+X.*sqrt((sigma-1).^2+4*sigma.*rho)) ./(2*sigma);
surf(X,Y,Z, 'EdgeColor','none','FaceColor','y', 'FaceAlpha',0.5)

lt1 = light('Position',[0 0 1], 'Style','infinite');
lt1.Color = 'w';
lighting gouraud
lt2 = light('Position',[0 0 -1], 'Style','infinite');
lt2.Color = 'w';
lighting gouraud
% 面的交线*******
xn = linspace(0.1,20,n);
yn = zeros(1,n);
zn = (rho).*ones(1,n);
plot3(xn, yn, zn, 'y'); %S+与P+交线
xn = linspace(-20,-0.1,n);
plot3(xn, yn, zn, 'y'); %S-与P-交线
xn = linspace(0.1,20,n);
yn = kcp(0).*xn;
zn = rho-kcp(0).*abs(xn)./Q;
plot3(xn, yn, zn, 'y'); %T+与P+交线
xn = linspace(-20,-0.1,n);
yn = kcp(0).*xn;
zn = rho-kcp(0).*abs(xn)./Q;
plot3(xn, yn, zn, 'y'); %T-与P-交线
%  % ---*******
% % % zlim([0 40]) % 40 45 50
% % % zticks([0 20 40]) %[0 20 40]  [0 15 30 45] [0 25 50]  
view(40,20); % ---123456;
%% 绘制检测点
% 录制初始帧
frame = getframe(gcf);
writeVideo(vid,frame);
% 填充静止帧（24fps×0.5=12帧）
frame_count = ptim * vid.FrameRate;  % 需填充的帧数=暂停时间×帧率
for iu = 1:frame_count
    writeVideo(vid,frame);  % 重复写入静止帧，填充暂停时间
end

ii = 1; ij = 1;
while ~isempty(mps) && size(mps,1)>1 && ij<=size(mps,1) % 2*sigma*data(mps(ij),2)>(sigma-1+sqrt((sigma-1)^2+4*sigma*rho))*data(mps(ij),1)
    if data(mps(ij),1)>0 && Vp(data(mps(ij),1),data(mps(ij),2),data(mps(ij),3))>0 ...
            || data(mps(ij),1)<0 && Vp(data(mps(ij),1),data(mps(ij),2),data(mps(ij),3))<0
        if ij==1
            plot3(data(ii:mps(ij),1), data(ii:mps(ij),2), data(ii:mps(ij),3),'b'); % tran
            plot3(data(mps(ij),1), data(mps(ij),2), data(mps(ij),3),'k*'); % 拐点
%             pause;
                drawnow;  % 刷新画面
                pause(ptim);
                % 将当前帧写入视频
                frame = getframe(fig);  % 获取当前绘图窗口的帧
                writeVideo(vid, frame);   % 写入视频
                for iu = 1:frame_count
                    writeVideo(vid,frame);  % 重复写入静止帧，填充暂停时间
                end
            continue;
        else
            [ipc,ipid]  = detectip(data(ii:mps(ij),1), data(ii:mps(ij),2), data(ii:mps(ij),3));
            if isempty(ipc)
                plot3(data(ii:mps(ij),1), data(ii:mps(ij),2), data(ii:mps(ij),3),'r'); % non-tran
            else
                ipid=ii+ipid;
                plot3(data(ii:ipid,1), data(ii:ipid,2), data(ii:ipid,3),'r'); % non-tran
                plot3(ipc(1), ipc(2), ipc(3),'k*'); % 拐点
%                 pause;
                    drawnow;  % 刷新画面
                    pause(ptim);
                    % 将当前帧写入视频
                    frame = getframe(fig);  % 获取当前绘图窗口的帧
                    writeVideo(vid, frame);   % 写入视频
                    for iu = 1:frame_count
                        writeVideo(vid,frame);  % 重复写入静止帧，填充暂停时间
                    end
                plot3(data(ipid:mps(ij),1), data(ipid:mps(ij),2), data(ipid:mps(ij),3),'b'); % tran
            end
        end
    else
        plot3(data(ii:mps(ij),1), data(ii:mps(ij),2), data(ii:mps(ij),3),'r'); % 解曲线
    end
    
    if Vn(data(mps(ij),1), data(mps(ij),2), data(mps(ij),3)) *Vn(data(mps(ij)+1,1), data(mps(ij)+1,2), data(mps(ij)+1,3))<=0
        plot3(data(mps(ij),1), data(mps(ij),2), data(mps(ij),3),'ko'); % S交点
    elseif Vp(data(mps(ij),1), data(mps(ij),2), data(mps(ij),3)) *Vp(data(mps(ij)+1,1), data(mps(ij)+1,2), data(mps(ij)+1,3))<=0
        plot3(data(mps(ij),1), data(mps(ij),2), data(mps(ij),3),'kx'); % T交点
    else
        plot3(data(mps(ij),1), data(mps(ij),2), data(mps(ij),3),'k.'); % P交点
    end
    ii = mps(ij);
    ij = ij+1;
%     pause;
        drawnow;  % 刷新画面
        pause(ptim);
        % 将当前帧写入视频
        frame = getframe(fig);  % 获取当前绘图窗口的帧
        writeVideo(vid, frame);   % 写入视频
        for iu = 1:frame_count
            writeVideo(vid,frame);  % 重复写入静止帧，填充暂停时间
        end
end

% 关闭视频写入对象
close(vid);
% close(fig);
% annotation('textbox',[.2 .3 .1 .1],'String','$T^{+}$', 'interpreter','latex', 'EdgeColor','none')
% annotation('textbox',[.2 .1 .1 .1],'String','$S^{-}$', 'interpreter','latex', 'EdgeColor','none')
% annotation('textbox',[.2 .3 .1 .1],'String','$P^{-}$', 'interpreter','latex', 'EdgeColor','none')
% annotation('textbox',[.2 .1 .1 .1],'String','$T^{-}$', 'interpreter','latex', 'EdgeColor','none')

%% 调用函数
function [cip,ipid] =  detectip(x, y, z) 
% 返回【拐点坐标，拐点索引】
tp1 = 1 +find(diff(sign(y), 1));
[cip,ipid] =  DIPS(x(tp1:end), y(tp1:end), z(tp1:end));
ipid = tp1+ipid-1; %在X中的索引
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



