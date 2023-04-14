close all
clear all
clc
C=[2 15 13 4
   10 4 14 15
   9 14 16 13
   7 8 11 9]
citys=1:size(C,1)+size(C,2);
D=[1000*ones(size(C,1),size(C,1)) C
   1000*ones(size(C,2),size(C,1)) 1000*ones(size(C,2),size(C,2))];
%初始化参数
m=10;%蚂蚁数量
arfa=1;%信息素重要程度因子
beita=5;%启发函数重要程度因子
rou=0.1;%信息素挥发因子
Q=1;
eta=1./D;%启发函数
T=ones(length(citys),length(citys));%信息素矩阵
Table=zeros(m,length(citys));%路径记录表
mingen=1;%迭代次数初始值
maxgen=100;%迭代次数最大值
zuiyou_Time=zeros(1,maxgen);%记录每次迭代的最优距离
zuiyou_route=zeros(maxgen,length(citys));%记录每次迭代的最优路线
while mingen<=maxgen
      start=zeros(m,1);
for i=1:m
      start(i)=randperm(length(citys),1);%m只蚂蚁的起始位置
end
Table(:,1)=start;%m只蚂蚁的起始位置放到路径记录表的第一列
for i=1:m
    for j=2:length(citys)
        tabu=Table(i,1:(j-1));%已访问的城市集合
        allow_index=~ismember(citys,tabu);
        allow=citys(allow_index);%没有走过的地点
        p=allow;
        for k=1:length(allow)
            p(k)=T(tabu(end),allow(k))^arfa*eta(tabu(end),allow(k))^beita; %城市转移概率
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
    route=Table(i,:);
    for j=2:length(citys)
        chang(i)=chang(i)+D(route(j-1),route(j));%计算m只蚂蚁的行走距离
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
   line([mingen-1,mingen],[zuiyou_chang(mingen-1)-3000,zuiyou_chang(mingen)-3000]);
end
xlabel('迭代次数')
ylabel('最优值')
title('迭代过程')
hold off
%更新信息素
derta=zeros(size(D,1),size(D,2));
 for i=1:m
     route=Table(i,:);
     for j=2:length(citys)
         derta(route(j-1),route(j))=derta(route(j-1),route(j))+Q/chang(i);%释放的信息素浓度      
     end
 end
 T=(1-rou)*T+derta;%更新信息素浓度
 mingen=mingen+1;
 Table=zeros(m,length(citys));
end
zuiduan_route=zuiyou_route(end,:)
Youmoney=0;
for i=2:length(citys)
    if D(zuiduan_route(i-1),zuiduan_route(i))~=1000
        Youmoney=Youmoney+D(zuiduan_route(i-1),zuiduan_route(i));
    end
end
You=zeros(size(C,1),size(C,2));
for i=1:size(C,1)
    You(zuiduan_route(2*i-1),zuiduan_route(2*i)-size(C,1))=1;
end
disp(['最优解为：',num2str(Youmoney)]);
disp('最优的指派方案为：')
You

