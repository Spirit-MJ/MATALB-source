function main
close all
clear all
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
N=length(x);
for i=1:N
    for j=1:N
        if i==j
        D(i,j)=eps;
        else 
        D(i,j)=sqrt((x(i,1)-x(j,1))^2+(x(i,2)-x(j,2))^2);%距离矩阵
        end
    end
end

T0=1e30;%初始温度
T_end=1e-30;%终止温度
L=10;%各温度下的迭代次数
Q=0.95;%降温速率
Time=ceil(double(log(T_end/T0)/log(Q)));%迭代次数
count=0;%计数
obj=zeros(Time,1);%距离记录
track=zeros(Time,N);%路线记录
for i=1:L
    S11(i,:)=randperm(N);%每个温度下产生L个可行解
    Dd(i)=juli(S11(i,:),x);
end
S1=XuanZe(S11,Dd);%选择最优的解
%S1=randperm(N);
figure(1)
for i=1:N
S1x(i)=x(S1(i),1);%随机产生一个初始路线
S1y(i)=x(S1(i),2);%随机产生一个初始路线
end
S1x=[S1x x(S1(1),1)];%把开始地点添加到最后画图
S1y=[S1y x(S1(1),2)];%把开始地点添加到最后画图
plot(S1x,S1y,'r-o');%初始解的图像
for i=1:N
    text(x(i,1),x(i,2),['   ' num2str(i)]);%在图像上把每个点序号标上
end
text(x(S1(1),1),x(S1(1),2),'     起点','FontSize',13,'Color',[0 0 1]);%在图像上标出起点
text(x(S1(end),1),x(S1(end),2),'     终点','FontSize',13,'Color',[0 0 1]);%在图像上标出终点
D1=juli(S1,x);
SS1=num2str(S1(1));
for i=2:N
    SS1=[SS1,'―>',num2str(S1(i))];
end
disp('随机产生的路线为：');
disp(SS1);
disp(['随机产生的一个路线距离为：',num2str(D1)])
%迭代优化
while T0>T_end
    count=count+1;
    for i=1:L
        S11(i,:)=randperm(N);%每个温度下产生L个可行解
        Dd(i)=juli(S11(i,:),x);
    end
    S2=XuanZe(S11,Dd);%选择最优的解 
    a=round(rand(1,2).*(N-1)+1);  %产生2个1—31的随机数
    W=S2(a(1));
    S2(a(1))=S2(a(2));
    S2(a(2))=W;  %交换S2(a(1))与S2(a(2))的值，随机交换一个解中的两个位置
    [S1,R]=Metropolis(S1,S2,D,T0,x);
    if count==1||R<obj(count-1)
       obj(count)=R;
    else
       obj(count)=obj(count-1);
    end
    track(count,:)=S1;
    T0=Q*T0;
    S1=S2;   %将前一次路线迭代给此次
    figure(2)
    hold on
    if count>=2
       line([count-1,count],[obj(count-1),obj(count)]);
    end
end
xlabel('迭代次数');
ylabel('距离');
title('优化过程');
hold off
%plot(1:count,obj,'b')
SS2=num2str(track(end,1));
for i=2:N
    SS2=[SS2,'―>',num2str(track(end,i))];
end
disp('最优路线为：') 
disp(SS2);
D_end=juli(track(end,:),x);
disp(['总距离为：',num2str(D_end)]);
figure(3)
for i=1:N
S2x(i)=x(track(end,i),1);
S2y(i)=x(track(end,i),2);
end
S2x=[S2x x(track(end,1),1)];
S2y=[S2y x(track(end,1),2)];
plot(S2x,S2y,'b-o');
for i=1:size(x,1)
    text(x(i,1),x(i,2),['   ' num2str(i)]);
end
text(x(track(end,1),1),x(track(end,1),2),'     起点','FontSize',13,'Color',[1 0 0]);
text(x(track(end,end),1),x(track(end,end),2),'     终点','FontSize',13,'Color',[1 0 0]);
function D1=juli(S,X) %计算任意解的距离
N=length(X);
for i=1:N
S1x(i)=X(S(i),1);
S1y(i)=X(S(i),2);
end
S1x=[S1x X(S(1),1)];
S1y=[S1y X(S(1),2)];
D1=0;
for i=1:N
    D1=D1+sqrt((S1x(i)-S1x(i+1))^2+(S1y(i)-S1y(i+1))^2);
end
end
function [S,R]=Metropolis(S1,S2,D,T,X)
% S1：  当前解
% S2:   新解
% D:    距离矩阵（两两城市的之间的距离）
% T:    当前温度
% S：   下一个当前解
% R：   下一个当前解的路线距离
R1=juli(S1,X);  %计算前一个解路线长度
N=length(X);         %得到城市的个数
R2=juli(S2,X);  %计算此次解路线长度
dC=R2-R1;   %计算能力之差
if dC<0       %如果能力降低 接受新路线
    S=S2;
    R=R2;
elseif exp(-dC/T)>=rand(1)   %以exp(-dC/T)概率接受新路线
    S=S2;
    R=R2;
else        %不接受新路线
    S=S1;
    R=R1;
end
end
function S1=XuanZe(S11,D)%选择  
       [~,temp]=sort(D);%按照距离降序排列
       S1=S11(temp(1),:);%优胜略汰，将距离最近的个体保留下来
end
end
