close all
clear
clc
data_origin = xlsread('E:\桌面\经纬度高程.xls');
origin_start = [120.767,52.183,419];  % 起始点
origin_end = [121.413,52.360,675];  % 终点

[X,Y,Z]=griddata(data_origin(:,2),data_origin(:,3),data_origin(:,1), ...
    linspace(min(data_origin(:,2)),max(data_origin(:,2)),100)', ...
    linspace(min(data_origin(:,3)),max(data_origin(:,3)),100));  % 曲面拟合
figure(1)
mesh(X,Y,Z)
xlabel('经度');
ylabel('纬度');
zlabel('海拔');
title('拟合曲面')

Z(ismissing(Z)) = max(max(Z));

%  地图建模
LevelGrid = 500;  % 将高度划分成500份
LevelGrid_h = (max(max(Z)) - min(min(Z))) / LevelGrid;  % （海拔点与点之间相隔距离）

PortGrid = 100;  % 将x于y划分成100份
PortGrid_x = (max(data_origin(:,2)) - min(data_origin(:,2))) / PortGrid;  % （经度点与点之间相隔距离）
PortGrid_y = (max(data_origin(:,3)) - min(data_origin(:,3))) / PortGrid;  % （维度点与点之间相隔距离）

Grid_Z = floor((Z - min(data_origin(:,1))) / LevelGrid_h);  % 建模后相对海拔；
Grid_Z(Grid_Z == 501) = 500;
start_x = round((origin_start(1) - min(data_origin(:,2))) / PortGrid_x);  % 建模后起始点位置
start_y = round((origin_start(2) - min(data_origin(:,3))) / PortGrid_y);  % 建模后起始点位置
start_h = round((origin_start(3) - min(min(Z))) / LevelGrid_h);  % 建模后起始点位置
end_x = round((origin_end(1) - min(data_origin(:,2))) / PortGrid_x);  % 建模后终点位置
end_y = round((origin_end(2) - min(data_origin(:,3))) / PortGrid_y);  % 建模后终点位置
end_h = round((origin_end(3) - min(min(Z))) / LevelGrid_h);  % 建模后终点位置

%算法参数
ant_num = 50;  % 蚂蚁数量
rou = 0.2;
T = ones(PortGrid,PortGrid,LevelGrid);  % 初始化信息素
mingen = 2;  % 迭代次数初始值
maxgen = 300;  % 迭代次数最大值
zuiyou_chang = zeros(1,maxgen);  % 记录每次迭代的最优距离
zuiyou_route = zeros(maxgen, 2 * (end_x-start_x + 1));  % 记录每次迭代的最优距离

% 初始搜索路径
[path,T] = find_path(ant_num,T,PortGrid,start_x,start_y,start_h,end_x,end_y,end_h, ...
    Grid_Z,Z,data_origin,LevelGrid_h,PortGrid_x,PortGrid_y);  % 初始搜索路径
fitness = comput_fit(path,data_origin,PortGrid_x,PortGrid_y,LevelGrid_h,Z);  %适应度计算
[min_chang,min_index] = min(fitness);           %最佳适应度
zuiyou_route(1,:) = path(min_index,:);     %最佳路径
zuiyou_chang(1) = min_chang;               %适应度值记录

% 信息素更新
cfit = 100 / zuiyou_chang(1);
for i = 2:59
    T(i,zuiyou_route(1,i*2-1),zuiyou_route(1,i*2))= ...
        (1-rou)*T(i,zuiyou_route(1,i*2-1),zuiyou_route(1,i*2))+ rou * cfit;
end


while mingen <= maxgen
    % 路径搜索
    [path,T] = find_path(ant_num,T,PortGrid,start_x,start_y,start_h,end_x,end_y,end_h, ...
    Grid_Z,Z,data_origin,LevelGrid_h,PortGrid_x,PortGrid_y); 
    
    % 适应度值计算更新
    fitness = comput_fit(path,data_origin,PortGrid_x,PortGrid_y,LevelGrid_h,Z);                               
    [min_chang,min_index] = min(fitness);
    
    if min_chang < zuiyou_chang(mingen-1)  % 判断此次距离是否比上次距离短
         zuiyou_chang(mingen) = min_chang;  % 是的话记录此次最短距离
         zuiyou_route(mingen,:) = path(min_index,:);  % 是的话记录此次最优路线
    else
         zuiyou_chang(mingen) = zuiyou_chang(mingen-1);  % 否则的话记录上一次最短距离
         zuiyou_route(mingen,:) = zuiyou_route(mingen-1,:);  % 否则的话记录上一次最优路线
    end
    
    % 更新信息素
    cfit = 100 / zuiyou_chang(mingen);
    for i = 2:59
        T(i,zuiyou_route(mingen,i*2-1),zuiyou_route(mingen,i*2)) = (1-rou) * ...
            T(i,zuiyou_route(mingen,i*2-1),zuiyou_route(mingen,i*2)) + rou*cfit;
    end

    figure(2)  % 迭代过程图像
    hold on
    if mingen>=2
        line([mingen-1,mingen],[zuiyou_chang(mingen-1),zuiyou_chang(mingen)]);
    end
    xlabel('迭代次数')
    ylabel('距离/km')
    title('迭代过程')
    hold off
    mingen = mingen+ 1;
end
x_index = start_x:end_x;
true_x = x_index * PortGrid_x + min(min(data_origin(:,2)));
true_y = zuiyou_route(end,1:2:end) * PortGrid_y + min(min(data_origin(:,3)));
disp('最优路线为：')
for i = 1:length(true_x)
    true_h(i) = Z(zuiyou_route(end,2*i-1),x_index(i));
    disp([true_x(i),true_y(i),true_h(i)])
end
%  坐标转换
true_xyz_x = (6371000 + true_h) .* cos(true_y ./ 180 * pi) .* cos(true_x ./ 180 * pi);
true_xyz_y = (6371000 + true_h) .* cos(true_y ./ 180 * pi) .* sin(true_x ./ 180 * pi);
true_xyz_z = (6371000 + true_h) .* sin(true_y ./ 180 * pi);
% 计算欧氏距离
distance = 0;
for i = 1:length(true_xyz_x)-1
    distance = distance + sqrt((true_xyz_x(i)-true_xyz_x(i+1))^2+(true_xyz_y(i)-true_xyz_y(i+1))^2+(true_xyz_z(i)-true_xyz_z(i+1))^2);
end

disp(['总距离为：',num2str(distance/1000),'千米'])
disp(['总时间为：',num2str(distance/5/3600),'小时'])

figure(1)
hold on
plot3(true_x,true_y,true_h,'r-*','linewidth',2)
plot3(true_x(1),true_y(1),true_h(1),'--o','LineWidth',2,...
                       'MarkerEdgeColor','k',...
                       'MarkerFaceColor','g',...
                       'MarkerSize',10)
plot3(true_x(end),true_y(end),true_h(end),'--o','LineWidth',2,...
                       'MarkerEdgeColor','k',...
                       'MarkerFaceColor','g',...
                       'MarkerSize',10)
text(true_x(1),true_y(1),true_h(1),'起点');
text(true_x(end),true_y(end),true_h(end),'终点');


function [path, T] = find_path(ant_num,T,PortGrid,start_x,start_y,start_h,end_x,end_y,end_h,Grid_Z,Z,data_origin,LevelGrid_h,PortGrid_x,PortGrid_y)
% 搜索参数
y_max = 8;   % 蚂蚁纬度方向可视范围
path = zeros(ant_num, 2 * (end_x-start_x));
for ii = 1:ant_num

    path(ii,1:2) = [start_y, start_h];  % 记录起始坐标
    temp = 1;
    now_point=[start_y, start_h];      % 当前坐标
    
    % 计算点适应度值
    for loc_x = start_x+1:end_x-1  % 可视域为经度方向
        temp = temp + 1;
        %计算所有数据点对应的适应度值
        kk = 1;
        for i = -y_max:y_max
            if (now_point(1)+i<=PortGrid) && (now_point(1)+i>0)  % 保证可视域在曲面内
                to_j = Grid_Z(now_point(1)+i,loc_x);
                next_point(kk,:) = [now_point(1)+i,to_j];  % 沿着经度方向更新路线
                qfz(kk) = comput_qfz(next_point(kk,1),next_point(kk,2),now_point(1),now_point(2),end_y,end_h,loc_x,end_x,Z,data_origin,LevelGrid_h,PortGrid_x,PortGrid_y);  % 计算启发值
                qz(kk) = qfz(kk) * T(loc_x,next_point(kk,1),next_point(kk,2));  % 计算可视域中点的概率
                kk = kk+1;
            else
                qz(kk) = 0;  % 如果可视域在曲面外
                kk = kk + 1;
            end
        end
        
        %轮盘赌选择访问的下一个地点
        p = qz/sum(qz);
        pc = cumsum(p); 
        target_index = find(pc>=rand(1));  
        index = target_index(1);  % 找出蚂蚁下一个最可能去的点
        
        old_point = next_point(index,:);

        %更新信息素
        T(loc_x+1,old_point(1),old_point(2)) = 0.9 * T(loc_x+1,old_point(1),old_point(2));
        
        %路径保存
        path(ii,temp*2-1:temp*2) = [old_point(1),old_point(2)];
        now_point = old_point; 
    end
    path(ii,(temp+1)*2-1:(temp+1)*2) = [end_y,end_h];  % 终点
end
end

function qfz = comput_qfz(next_y,next_h,now_y,now_h,end_y,end_h,loc_x,end_x,Z,data_origin,LevelGrid_h,PortGrid_x,PortGrid_y)

% 经纬度与海拔真实值
true_now_x = (loc_x-1) * PortGrid_x + min(data_origin(:,2));
true_now_y = now_y * PortGrid_y + min(data_origin(:,3));
true_now_h = now_h * LevelGrid_h + min(min(Z));

true_next_x = loc_x * PortGrid_x + min(data_origin(:,2));
true_next_y = next_y * PortGrid_y + min(data_origin(:,3));
true_next_h = next_h * LevelGrid_h + min(min(Z));


true_end_x = end_x * PortGrid_x + min(data_origin(:,2));
true_end_y = end_y * PortGrid_y + min(data_origin(:,3));
true_end_h = end_h * LevelGrid_h + min(min(Z));


% 转换直角坐标
true_now_xyz_x = (6371000 + true_now_h) * cos(true_now_y / 180 * pi) * cos(true_now_x / 180 * pi);
true_now_xyz_y = (6371000 + true_now_h) * cos(true_now_y / 180 * pi) * sin(true_now_x / 180 * pi);
true_now_xyz_z = (6371000 + true_now_h) * sin(true_now_y / 180 * pi);

true_next_xyz_x = (6371000 + true_next_h) * cos(true_next_y / 180 * pi) * cos(true_next_x / 180 * pi);
true_next_xyz_y = (6371000 + true_next_h) * cos(true_next_y / 180 * pi) * sin(true_next_x / 180 * pi);
true_next_xyz_z = (6371000 + true_next_h) * sin(true_next_y / 180 * pi);

true_end_xyz_x = (6371000 + true_end_h) * cos(true_end_y / 180 * pi) * cos(true_end_x / 180 * pi);
true_end_xyz_y = (6371000 + true_end_h) * cos(true_end_y / 180 * pi) * sin(true_end_x / 180 * pi);
true_end_xyz_z = (6371000 + true_end_h) * sin(true_end_y / 180 * pi);

% 判断下个点是否可达

slope = abs(true_now_xyz_z-true_next_xyz_z) / sqrt((true_now_xyz_x-true_next_xyz_x)^2+(true_now_xyz_y-true_next_xyz_y)^2);
if (atan(slope) / pi * 180) <= 30  % 坡度约束
    S = 1;
else
    S = 0;
end

%D距离
D = 10000 / (sqrt((true_next_xyz_x-true_now_xyz_x)^2 + (true_next_xyz_y-true_now_xyz_y)^2+(true_next_xyz_z-true_now_xyz_z)^2) + ...  
    sqrt((true_end_xyz_x-true_next_xyz_x)^2+(true_end_xyz_y-true_next_xyz_y)^2+(true_end_xyz_z-true_next_xyz_z)^2));  % 下一可行点到终点的距离 
%计算启发值
qfz = D * S;
end

function fitness = comput_fit(path,data_origin,PortGrid_x,PortGrid_y,LevelGrid_h,Z)
[n,m] = size(path); 
fitness = zeros(n,1);
for i = 1:n
    fitness(i) = 0;
    for j = 2:m/2
        
        % 经纬度与海拔真实值
        true_path_x0 =  PortGrid_x + min(data_origin(:,2));
        true_path_x1 =  2 * PortGrid_x + min(data_origin(:,2));
        
        true_path_y0 = path(i,j*2-1) * PortGrid_y + min(data_origin(:,3));
        true_path_y1 = path(i,(j-1)*2-1) * PortGrid_y + min(data_origin(:,3));
    
        true_path_h0 = path(i,j*2) * LevelGrid_h + min(min(Z));
        true_path_h1 = path(i,(j-1)*2) * LevelGrid_h + min(min(Z));
        
        % 转换直角坐标
        true_xyz_x0 = (6371000 + true_path_h0) * cos(true_path_y0 / 180 * pi) * cos(true_path_x0 / 180 * pi);
        true_xyz_y0 = (6371000 + true_path_h0) * cos(true_path_y0 / 180 * pi) * sin(true_path_x0 / 180 * pi);
        true_xyz_z0 = (6371000 + true_path_h0) * sin(true_path_y0 / 180 * pi);
        
        true_xyz_x1 = (6371000 + true_path_h1) * cos(true_path_y1 / 180 * pi) * cos(true_path_x1 / 180 * pi);
        true_xyz_y1 = (6371000 + true_path_h1) * cos(true_path_y1 / 180 * pi) * sin(true_path_x1 / 180 * pi);
        true_xyz_z1 = (6371000 + true_path_h1) * sin(true_path_y1 / 180 * pi);
        
        %适应度值为长度加高度
        fitness(i) = fitness(i)+ sqrt((true_xyz_x1-true_xyz_x0)^2 + (true_xyz_y1-true_xyz_y0)^2 + (true_xyz_z1-true_xyz_z0)^2)/1000;  % + abs(path(i,j*2));
    end
end
end

