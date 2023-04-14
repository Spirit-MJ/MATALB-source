close all
clear
clc
x=[1304 2312;
   3639 1315;
   4177 2244;
   3712 1399;
   3488 1535;
   3326 1556;
   3238 1229;
   4196 1004;
   4312 790;
   4386 570;
   3007 1970;
   2562 1756;
   2788 1491;
   2381 1676;
   1332 695;
   3715 1678;
   3918 2179;
   4061 2370;
   3780 2212;
   3676 2578;
   4029 2838;
   4263 2931;
   3429 1908;
   3507 2367;
   3394 2643;
   3439 3201;
   2935 3240;
   3140 3550;
   2545 2357;
   2778 2826;
   2370 2975];

for i=1:length(x)
    for j=1:length(x)
        if i==j
        D(i,j)=eps;
        else 
        D(i,j)=sqrt((x(i,1)-x(j,1))^2+(x(i,2)-x(j,2))^2);%距离
        end
    end
end
%初始化参数
m=35;%蚂蚁数量
arfa=1;%信息素重要程度因子
beita=2;%启发函数重要程度因子
rou=0.2;%信息素挥发因子
Q=1;
eta=1./D;%启发函数
T=ones(length(x));%信息素矩阵
Table=zeros(m,length(x));%路径记录表
mingen=1;%迭代次数初始值
maxgen=300;%迭代次数最大值
zuiyou_chang=zeros(1,maxgen);%记录每次迭代的最优距离
zuiyou_route=zeros(maxgen,length(x));%记录每次迭代的最优路线
while mingen<=maxgen
    start=zeros(m,1);
for i=1:m
    start(i)=randperm(length(x),1);%m只蚂蚁的起始位置
end
    Table(:,1)=start;%m只蚂蚁的起始位置放到路径记录表的第一列
    citys=1:length(x);
for i=1:m
    for j=2:length(x)
        tabu=Table(i,1:(j-1));%已访问的城市集合
        allow_index=~ismember(citys,tabu);
        allow=citys(allow_index);%没有走过的地点
        p=allow;
        for k=1:length(allow)
            p(k)=T(tabu(end),allow(k))^arfa * eta(tabu(end),allow(k))^beita; %城市转移概率
        end
        p=p/sum(p);
        pc=cumsum(p);%轮盘赌
        mubiao_index=find(pc>=rand(1));
        mubiao=allow(mubiao_index(1));%找出蚂蚁下一个最可能去的点
        Table(i,j)=mubiao;%记录每只的蚂蚁路线
    end
end
chang=zeros(m,1);
for i=1:m
    route=[Table(i,:) Table(i,1)];
    for j=1:length(x)
        chang(i)=chang(i)+D(route(j),route(j+1));%计算m只蚂蚁的行走距离
    end
end
%更新路线
 if mingen==1
     [min_chang,min_index]=min(chang);%找出第一次m只蚂蚁走的最短路径
     zuiyou_chang(mingen)=min_chang;%记录每一次最短的距离
     zuiyou_route(mingen,:)=Table(min_index,:);%记录每一次最优的路线
 else 
     [min_chang, min_index]=min(chang);
     if min_chang<zuiyou_chang(mingen-1)%判断此次距离是否比上次距离短
         zuiyou_chang(mingen)=min_chang;%是的话记录此次最短距离
         zuiyou_route(mingen,:)=Table(min_index,:);%是的话记录此次最优路线
     else
         zuiyou_chang(mingen)=zuiyou_chang(mingen-1);%否则的话记录上一次最短距离
         zuiyou_route(mingen,:)=zuiyou_route(mingen-1,:);%否则的话记录上一次最优路线
     end
 end
    figure(1)%迭代过程图像
    hold on
    if mingen>=2
        line([mingen-1,mingen],[zuiyou_chang(mingen-1),zuiyou_chang(mingen)]);
    end
    xlabel('迭代次数')
    ylabel('最优值')
    title('迭代过程')
    hold off
%更新信息素
    derta=zeros(length(x));
     for i=1:m
         route=[Table(i,:) Table(i,1)];
         for j=1:length(x)
             derta(route(j),route(j+1))=derta(route(j),route(j+1))+Q/chang(i);%释放的信息素浓度      
         end
     end
     T=(1-rou)*T+derta;%更新信息素浓度
     mingen=mingen+1;
     Table=zeros(m,length(x));
end
zuiduan_route=zuiyou_route(end,:);
zuiduan_chang=zuiyou_chang(end);
disp(['最短的路径为：',num2str(zuiduan_route)]);
disp(['最短的距离为：',num2str(zuiduan_chang)]);
for i=1:length(x)
    cityx(i)=x(zuiduan_route(i),1);
    cityy(i)=x(zuiduan_route(i),2);
end
cityx=[cityx x(zuiduan_route(1),1)];
cityy=[cityy x(zuiduan_route(1),2)];
hold off
figure(2)
plot(cityx,cityy,'r-*');
xlabel('城市位置横坐标')
ylabel('城市位置纵坐标')
title('最短路径')

