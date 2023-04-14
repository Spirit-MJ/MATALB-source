close all
clear
clc
city=[1304 2312
   3639 1315
   4177 2244
   3712 1399
   3488 1535
   3326 1556
   3238 1229
   4196 1004
   4312 790
   4386 570
   3007 1970
   2562 1756
   2788 1491
   2381 1676
   1332 695
   3715 1678
   3918 2179
   4061 2370
   3780 2212
   3676 2578
   4029 2838
   4263 2931
   3429 1908
   3507 2367
   3394 2643
   3439 3201
   2935 3240
   3140 3550
   2545 2357
   2778 2826
   2370 2975];%城市坐标矩阵
figure(1)
plot(city(:,1),city(:,2),'ms','LineWidth',2,'MarkerEdgeColor','k','MarkerFaceColor','g')
legend('城市位置')
title('城市分布图','fontsize',12)
xlabel('km','fontsize',12)
ylabel('km','fontsize',12)
grid on
%计算城市间距离
n=size(city,1);%城市数目
D=zeros(n,n);%城市距离矩阵
D=Juli(city);
nMax=200;                      %进化次数
indiNumber=5000;               %个体数目
individual=zeros(indiNumber,n);
%初始化粒子位置
for i=1:indiNumber
    individual(i,:)=randperm(n);    
end
%计算种群适应度
fitx=fitness(individual,city,D);
[~,index]=max(fitx);
tourPbest=individual;                              %当前个体最优
tourGbest=individual(index,:);                    %当前全局最优
recordPbest=zeros(1,indiNumber);                %个体最优记录
recordGbest=fitx(index);                        %群体最优记录
xnew1=individual;
% 循环寻找最优路径
L_best=zeros(1,nMax);
for N=1:nMax 
    % 交叉操作
    for i=1:indiNumber
        shu=randperm(n-1,2);
        chb1=min(shu);
        chb2=max(shu);
        cros=tourPbest(i,chb1:chb2);
        ncros=size(cros,2);      
        %删除与交叉区域相同元素
        for j=1:ncros
            for k=1:n
                if xnew1(i,k)==cros(j)
                    xnew1(i,k)=0;
                    for t=1:n-k
                        temp=xnew1(i,k+t-1);
                        xnew1(i,k+t-1)=xnew1(i,k+t);
                        xnew1(i,k+t)=temp;
                    end
                end
            end
        end
        %插入交叉区域
        xnew1(i,n-ncros+1:n)=cros;
        %新路径长度变短则接受
        dist=0;
        for j=1:n-1
            dist=dist+D(xnew1(i,j),xnew1(i,j+1));
        end
        dist=dist+D(xnew1(i,1),xnew1(i,n));
        if (1/dist)>fitx(i)
            individual(i,:)=xnew1(i,:);
        end
        %与全体最优进行交叉
        shu=randperm(n-1,2);
        chb1=min(shu);
        chb2=max(shu);
        cros=tourGbest(chb1:chb2); 
        ncros=size(cros,2);      
        %删除与交叉区域相同元素
        for j=1:ncros
            for k=1:n
                if xnew1(i,k)==cros(j)
                    xnew1(i,k)=0;
                    for t=1:n-k
                        temp=xnew1(i,k+t-1);
                        xnew1(i,k+t-1)=xnew1(i,k+t);
                        xnew1(i,k+t)=temp;
                    end
                end
            end
        end
        %插入交叉区域
        xnew1(i,n-ncros+1:n)=cros;
        %新路径长度变短则接受
        dist=0;
        for j=1:n-1
            dist=dist+D(xnew1(i,j),xnew1(i,j+1));
        end
        dist=dist+D(xnew1(i,1),xnew1(i,n));
        if (1/dist)>fitx(i)
            individual(i,:)=xnew1(i,:);
        end
       %变异操作
        shu=randperm(n-1,2);
        c1=min(shu);
        c2=max(shu);
        temp=xnew1(i,c1);
        xnew1(i,c1)=xnew1(i,c2);
        xnew1(i,c2)=temp;
        %新路径长度变短则接受
        dist=0;
        for j=1:n-1
            dist=dist+D(xnew1(i,j),xnew1(i,j+1));
        end
        dist=dist+D(xnew1(i,1),xnew1(i,n));
        if (1/dist)>fitx(i)
            individual(i,:)=xnew1(i,:);
        end
        %进化逆转操作
        r=randperm(n,2);%产生两个随机数
        mininverse=min(r);
        maxinverse=max(r);
        xnew1(i,mininverse:maxinverse)=xnew1(i,maxinverse:-1:mininverse);%逆转操作,基因对换位置
        %新路径长度变短则接受
        dist=0;
        for j=1:n-1
            dist=dist+D(xnew1(i,j),xnew1(i,j+1));
        end
        dist=dist+D(xnew1(i,1),xnew1(i,n));
        if (1/dist)>fitx(i)
            individual(i,:)=xnew1(i,:);
        end
    end 
fitx=fitness(individual,city,D);%计算适应度值
%更新当前最优和历史最优
for i=1:indiNumber
    if fitx(i)>recordPbest(i)
       recordPbest(i)=fitx(i);
       tourPbest(i,:)=individual(i,:);
    end
    if fitx(i)>recordGbest
       recordGbest=fitx(i);
       tourGbest=individual(i,:);
    end
end
L_best(N)=recordGbest;
if N>=2
   figure(2)
   line([N-1,N],[1/L_best(N-1),1/L_best(N)]);% 结果作图
   title('算法训练过程')
   xlabel('迭代次数')
   ylabel('距离')
   grid on
end
end

figure(3)
hold on
plot([city(tourGbest(1),1),city(tourGbest(n),1)],[city(tourGbest(1),2),...
    city(tourGbest(n),2)],'ms-','LineWidth',2,'MarkerEdgeColor','k','MarkerFaceColor','g')
hold on
for i=2:n
    plot([city(tourGbest(i-1),1),city(tourGbest(i),1)],[city(tourGbest(i-1),2),...
        city(tourGbest(i),2)],'ms-','LineWidth',2,'MarkerEdgeColor','k','MarkerFaceColor','g')
    hold on
end
legend('规划路径')
scatter(city(:,1),city(:,2));
title('规划路径','fontsize',10)
xlabel('km','fontsize',10)
ylabel('km','fontsize',10)
grid on

function D=Juli(city)
n=size(city,1);
for i=1:n
    for j=i:n
        if i==j
           D(i,j)=eps;
        else
            D(i,j)=((city(i,1)-city(j,1))^2+(city(i,2)-city(j,2))^2)^0.5;
        end
        D(j,i)=D(i,j);
    end
end
end
function fitx=fitness(x,city,D)
%x           input     个体
%cit         input     城市坐标
%D           input     城市距离
%fitx        output    个体适应度值 
m=size(x,1);
n=size(city,1);
fx=zeros(m,1);
for i=1:m
    for j=1:n-1
        fx(i)=fx(i)+D(x(i,j),x(i,j+1));
    end
    fx(i)=fx(i)+D(x(i,1),x(i,n));
end
fitx=1./fx;
end