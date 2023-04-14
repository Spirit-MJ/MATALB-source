close all
clear
clc

figure(1)
plot([1,200],[1,200],'.'); % 规划空间
hold on
xlabel('km','fontsize',12)
ylabel('km','fontsize',12)
title('二维规划空间','fontsize',12)
% 描述起点和终点
Start = [20,180];
End = [180,110];
plot([Start(1),End(1)],[Start(2),End(2)],'.');

% 图形标注
text(Start(1)-15,Start(2),'起点');  % 起点
text(End(1)+2,End(2),'终点');  % 终点

position = [40,120;40,190;80,190;80,120;
            100,120;100,150;140,150;140,120;
            150,100;150,120;170,120;170,100;
            80,80;80,110;120,110;120,80];  % 障碍物
        
% 描绘障碍物图形
fill(position(1:4,1),position(1:4,2),[0,0,0]);
fill(position(5:8,1),position(5:8,2),[0,0,0]);
fill(position(9:12,1),position(9:12,2),[0,0,0]);
fill(position(13:16,1),position(13:16,2),[0,0,0]);

%算法参数
ant_num = 100;  % 蚂蚁数量
rou = 0.2;
view_h = 10;  % 横向可视域大小
view_w_min = 90;  % 纵向可视域最小值
view_w_max = 180;  % 纵向可视域最大值
view_w = view_w_max - view_w_min + 1;  % 纵向可视域大小 
num_gen = floor((End(1)-Start(1)-1)/view_h);  % 需访问的点x的总数

T_fiil0 = ones(200);  % 初始化信息素
for ii = 1:length(position)/4
    for i = position(ii*4-3,1):position(ii*4-1,1)
        for j = position(ii*4-3,2):position(ii*4-1,2)
            T_fiil0(i,j) = 0;
        end
    end
end  
T_fill = zeros(End(1)-Start(1)-1, view_w);
for i0 = 1:size(T_fill,1)
    T_fill(i0,:) = T_fiil0(Start(1)+i0,view_w_min:view_w_max);
end
T = ones(num_gen, view_w);  % 初始化信息素
mingen = 2;  % 迭代次数初始值
maxgen = 500;  % 迭代次数最大值
zuiyou_chang = zeros(1,maxgen);  % 记录每次迭代的最优距离
zuiyou_route = zeros(maxgen,num_gen+2);  % 记录每次迭代的最优路线

% 初始搜索路径
[path,T] = find_path(ant_num, num_gen, Start, End, view_h, T, T_fill, view_w_min, view_w_max);  % 初始搜索路径
fitness = comput_fit(path,view_h);  %适应度计算
[min_chang,min_index] = min(fitness);           %最佳适应度
zuiyou_route(1,:) = path(min_index,:);     %最佳路径
zuiyou_chang(1) = min_chang;               %适应度值记录

% 信息素更新
cfit = 1 / zuiyou_chang(1);
mm = 1;
for i = 2:length(zuiyou_route(1,:))-1
    T(mm,zuiyou_route(1,i)-view_w_min+1) = (1-rou)*T(mm,zuiyou_route(1,i)-view_w_min+1)+ rou * cfit;
    mm = mm + 1;
end

while mingen <= maxgen
    % 路径搜索
    [path,T] = find_path(ant_num, num_gen, Start, End, view_h, T, T_fill, view_w_min, view_w_max);  % 初始搜索路径
    
    % 适应度值计算更新
    fitness = comput_fit(path,view_h);  %适应度计算
    [min_chang,min_index] = min(fitness);           %最佳适应度
    
    if min_chang < zuiyou_chang(mingen-1)  % 判断此次距离是否比上次距离短
         zuiyou_chang(mingen) = min_chang;  % 是的话记录此次最短距离
         zuiyou_route(mingen,:) = path(min_index,:);  % 是的话记录此次最优路线
    else
         zuiyou_chang(mingen) = zuiyou_chang(mingen-1);  % 否则的话记录上一次最短距离
         zuiyou_route(mingen,:) = zuiyou_route(mingen-1,:);  % 否则的话记录上一次最优路线
    end
    
    % 更新信息素
    cfit = 1 / zuiyou_chang(mingen);
    mm = 1;
    for i = 2:length(zuiyou_route(mingen,:))-1
        T(mm,zuiyou_route(mingen,i)-view_w_min+1) = (1-rou)*T(mm,zuiyou_route(mingen,i)-view_w_min+1)+ rou * cfit;
        mm = mm + 1;
    end

    figure(2)  % 迭代过程图像
    hold on
    line([mingen-1,mingen],[zuiyou_chang(mingen-1),zuiyou_chang(mingen)]);
    xlabel('迭代次数')
    ylabel('距离')
    title('迭代过程')
    hold off
    mingen = mingen+ 1;
end
loc_x = Start(1);
for i = 1:length(zuiyou_route(end,:))-2
    tenp = Start(1) + view_h * i;
    loc_x = [loc_x,tenp];
end
loc_x = [loc_x,End(1)];
figure(1)
plot(loc_x,zuiyou_route(end,:),'color','red','LineWidth',3,'LineStyle','-.'); % 规划空间
disp(['最优长为：',num2str(zuiyou_chang(end))])

function [path, T] = find_path(ant_num, num_gen, Start, End, view_h, T, T_fill, view_w_min, view_w_max)
path = zeros(ant_num, num_gen+2);
for ii = 1:ant_num

    path(ii,1) = Start(2);  % 记录起始坐标
    temp = 1;
    now_point = Start(2);      % 当前坐标
    
    % 计算点适应度值
    
    for loc_x = 1:num_gen  % 横向可视域
        temp = temp + 1;
        %计算所有数据点对应的适应度值
        kk = 1;
        for i = view_w_min:view_w_max  % 纵向搜索
            next_point(kk) = i;  % 沿着纵向更新路线
            qfz(kk) = comput_qfz(Start(1)+(loc_x-1)*view_h, now_point, next_point(kk),Start(1),view_w_min,End(1),End(2),view_h,T_fill);  % 计算启发值
            qz(kk) = qfz(kk) * T(loc_x,next_point(kk)-view_w_min+1) * T_fill(loc_x*view_h,next_point(kk)-view_w_min+1);  % 计算可视域中点的概率
            kk = kk+1;
        end
        %轮盘赌选择访问的下一个地点
        p = qz/sum(qz);
        pc = cumsum(p); 
        target_index = find(pc>=rand(1));  
        index = target_index(1);  % 找出蚂蚁下一个最可能去的点
        
        old_point = next_point(index);

        %更新信息素
        T(loc_x,old_point-view_w_min+1) = 0.9 * T(loc_x,old_point-view_w_min+1);  % 信息素挥发
        
        %路径保存
        path(ii,temp) = old_point;
        now_point = old_point; 
    end
    path(ii,temp+1) = End(2);  % 终点
end
end

function qfz = comput_qfz(now_x,now_y,next_y,start_x,start_y,end_x,end_y,view_h,T_fill)
next_x = now_x + view_h;
center_x = now_x + ceil(view_h/2) - start_x;
center_y = min(next_y,now_y) + ceil(abs(next_y - now_y)/2) - start_y +1;
if T_fill(center_x, center_y) == 0  % 避障
    S = 0;
else
    S = 1;
end

%D距离
D = 1 / (sqrt(view_h^2 + (now_y-next_y)^2) + ...   % 当前点与下一点距离
    sqrt((next_x-end_x)^2+(next_y-end_y)^2));  % 下一可行点到终点的距离 
%计算启发值
qfz = D * S;
end

function fitness = comput_fit(path,view_h)
[n,m] = size(path); 
fitness = zeros(n,1);
for i = 1:n
    fitness(i) = 0;
    for j = 2:m
        fitness(i) = fitness(i)+ sqrt(view_h^2 + (path(i,j)-path(i,j-1))^2);  
    end
end
end
