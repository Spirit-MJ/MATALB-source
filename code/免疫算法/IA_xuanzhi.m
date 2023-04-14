% 免疫优化算法在物流配送中心选址中的应用
close all
clear all
clc
city_coordinate=[1304,2312;
                 3639,1315;
                 4177,2244;
                 3712,1399;
                 3488,1535;
                 3326,1556;
                 3238,1229;
                 4196,1044;
                 4312,790;
                 4386,570;
                 3007,1970;
                 2562,1756;
                 2788,1491;
                 2381,1676;
                 1332,695;
                 3715,1678;
                 3918,2179;
                 4061,2370;
                 3780,2212;
                 3676,2578;
                 4029,2838;
                 4263,2931;
                 3429,1908;
                 3507,2376;
                 3394,2643;
                 3439,3201;
                 2935,3240;
                 3140,3550;
                 2545,2357;
                 2778,2826;
                 2370,2975];
carge=[20,90,90,60,70,70,40,90,90,70,60,40,40,40,20,80,90,70,100,50,50,50,80,70,80,40,40,60,70,50,30];           
np=100;           % 种群规模
overbest=10;          % 记忆库容量
ng=200;            % 迭代次数
pc=0.95;           % 交叉概率
pm=0.2;        % 变异概率
ps=0.9;              % 多样性评价参数
T=0.7;         % 相似度大于阀值，这里阈值为0.7
CenterNum=6;             % 配送中心数
M=np+overbest;      
%识别抗原,将种群信息定义为一个结构体
individuals=struct('fitx',zeros(1,M), 'C',zeros(1,M),'P',zeros(1,M),'x',[]);
%产生初始抗体群
individuals.x=popinit(M,CenterNum);%初始化种群
trace=[]; %记录每代最个体优适应度和平均适应度
for iii=1:ng
     %抗体群多样性评价
     for i=1:M
         individuals.fitx(i)=fitness(individuals.x(i,:),city_coordinate,carge);%计算抗体与抗原亲和力(适应度值）
         individuals.C(i)=concentration(i,M,individuals,T);%计算抗体浓度
     end
     %综合亲和度和浓度评价抗体优秀程度，得出繁殖概率
     individuals.P=excellence(individuals,M,ps);    
     %记录当代最佳个体和种群平均适应度
     [best,index]=max(individuals.fitx);   % 找出最优适应度 
     bestchrom=individuals.x(index,:);    % 找出最优个体
     average=mean(individuals.fitx);       % 计算平均适应度
     trace=[trace;1/best,1/average];              %记录
     %根据excellence，形成父代群，更新记忆库（加入精英保留策略，可由s控制）
     bestindividuals=member(individuals,M,overbest);   % 更新记忆库
     individuals=member(individuals,M,np);      % 形成父代群
     % 选择，交叉，变异操作，再加入记忆库中抗体，产生新种群
     individuals=XuanZe(individuals,np);% 选择
     individuals.x=JiaoCha(pc,individuals.x,np,CenterNum);% 交叉
     individuals.x=BianYi(pm,individuals.x,np,CenterNum);% 变异
     individuals=Reins(individuals,np,bestindividuals,overbest);% 加入记忆库中抗体      
     if iii>=2
      % 画出免疫算法收敛曲线
        figure(1)
        hold on
        draw1=line([iii-1,iii],[trace(iii-1,1),trace(iii,1)]);
        draw2=line([iii-1,iii],[trace(iii-1,2),trace(iii,2)]);
        grid on
        set(draw1,'color',[0 0 1]);
        set(draw2,'color',[1 0 0]);
        xlabel('迭代次数','fontsize',12)
        ylabel('适应度值','fontsize',12)
        title('免疫算法收敛曲线','fontsize',12)
        legend('最优适应度值','平均适应度值')
     end
end
% 画出配送中心选址图
for i=1:31
    distance(i,:)=dist(city_coordinate(i,:),city_coordinate(bestchrom,:)');
end
[a,b]=min(distance');
index=cell(1,CenterNum);
for i=1:CenterNum
    index{i}=find(b==i);%计算各个派送点的地址
end
cargox=city_coordinate(bestchrom,1);
cargoy=city_coordinate(bestchrom,2);
figure(2)
for j=1:length(index)
    for i=1:length(index{j})
        A=[city_coordinate(index{j}(i),1),city_coordinate(index{j}(i),2)];
        B=[cargox(j),cargoy(j)];
        c=[A;B];
        plot(c(:,1),c(:,2),'b-*')
        hold on
    end
end
title('免疫优化算法选址');
xlabel('x坐标');
ylabel('y坐标');
function pop=popinit(n,len)
for i=1:n
    flag=0;
    while flag==0
        pop(i,:)=randperm(31,len);
        flag=test(pop(i,:));
    end
end
end
function flag=test(code)
city_coordinate=[1304,2312;3639,1315;4177,2244;3712,1399;3488,1535;3326,1556;3238,1229;4196,1044;4312,790;4386,570;
                 3007,1970;2562,1756;2788,1491;2381,1676;1332,695;3715,1678;3918,2179;4061,2370;3780,2212;3676,2578;
                 4029,2838;4263,2931;3429,1908;3507,2376;3394,2643;3439,3201;2935,3240;3140,3550;2545,2357;2778,2826;2370,2975];
flag=1;
if max(max(dist(city_coordinate(code,:)')))>3000%配送中心距离约束
    flag=0;
end
end
function fitx=fitness(individual,city_coordinate,carge)
%找出最近配送点
%dist(A,B)计算A中每个行向量与B中每个列向量之间欧氏距离，A的行向量维数必须等于B的列向量维数
for i=1:31
    distance(i,:)=dist(city_coordinate(i,:),city_coordinate(individual,:)');
end
a=min(distance');%找出每个点到六个中心距离最短的点和索引号
for i=1:31  
    expense(i)=carge(i)*a(i);%计算费用
end
fx=sum(expense)+4.0e+4*length(find(a>3000));%距离大于3000取一个惩罚值
fitx=1/fx;
end
function Cv=concentration(i,M,individuals,T)%计算抗体浓度
% 计算个体浓度值
% i              input      第i个抗体
% M              input      种群规模
% individuals    input     个体
% Concentration  output     浓度值
Cv=0;
for j=1:M
    xsd=similar(individuals.x(i,:),individuals.x(j,:));%第i个体与种群个体间的相似度 
    if xsd>T
        Cv=Cv+1;
    end
end
Cv=Cv/M;
end
function Svs=similar(v,s)%计算抗体与抗体之间的亲和力
k=zeros(1,length(v));
for i=1:length(v)
    if find(v(i)==s)
        k(i)=1;
    end
end
Svs=sum(k)/length(v);
end
function P=excellence(individuals,M,ps)%计算期望繁殖概率
% 计算个体繁殖概率
% individuals    input      种群
% M              input      种群规模
% ps             input      多样性评价参数
% exc            output     繁殖概率
A=individuals.fitx;
sumA=sum(A);
C=individuals.C;
sumC=sum(C);
for i=1:M
    P(i)=ps*A(i)/sumA+(1-ps)*C(i)/sumC; 
end
end
function rets=member(individuals,m,n)%记录期望繁殖概率大的和适应度函数大的抗体
% 初始化记忆库,依据excellence，将群体中高适应度低相似度的overbest个个体存入记忆库
% m                  input          抗体数
% n                  input          记忆库个体数\父代群规模
% individuals        input          抗体群
% bestindividuals    output         记忆库\父代群
% 精英保留策略，将fitness最好的s个个体先存起来，避免因其浓度高而被淘汰
s=3;
rets=struct('fitx',zeros(1,n), 'C',zeros(1,n),'P',zeros(1,n),'x',[]);
[~,index]=sort(individuals.fitx,'descend');%按照适应度从大到小排列
for i=1:s
    rets.fitx(i)=individuals.fitx(index(i));   
    rets.C(i)=individuals.C(index(i));
    rets.P(i)=individuals.P(index(i));
    rets.x(i,:)=individuals.x(index(i),:);
end
% 剩余m-s个个体
leftindividuals=struct('fitness',zeros(1,m-s), 'concentration',zeros(1,m-s),'excellence',zeros(1,m-s),'chrom',[]);
for k=1:m-s
    leftindividuals.fitness(k)=individuals.fitx(index(k+s));   
    leftindividuals.concentration(k)=individuals.C(index(k+s));
    leftindividuals.excellence(k)=individuals.P(index(k+s));
    leftindividuals.chrom(k,:)=individuals.x(index(k+s),:);
end
% 将剩余抗体按excellence值排序,将期望繁殖概率大的个体保存下来
[~,index]=sort(leftindividuals.excellence,'descend');%按照期望繁殖概率从大到小排列
% 在剩余抗体群中按excellence再选n-s个最好的个体
for i=s+1:n
    rets.fitx(i)=leftindividuals.fitness(index(i-s));
    rets.C(i)=leftindividuals.concentration(index(i-s));
    rets.P(i)=leftindividuals.excellence(index(i-s));
    rets.x(i,:)=leftindividuals.chrom(index(i-s),:);
end
end
function ret=XuanZe(individuals,np)%选择
% 轮盘赌选择
Excellence=individuals.P;
pselect=Excellence./sum(Excellence);
for j=1:np
      sita=rand(1);
      while sita==0    
            sita=rand(1);        
      end
      for i=1:np
          if sita<=pselect(i)
             individuals.x(j,:)=individuals.x(i,:);
             individuals.fitx(j)=individuals.fitx(i);
             individuals.C(j)=individuals.C(i);
             individuals.P(j)=individuals.P(i);
             break;
           end
      end
end
ret=individuals;
end
function ret=JiaoCha(pcross,chrom,sizepop,len)%交叉
for i=1:sizepop   
    pick=rand(1);
    while pick==0
        pick=rand(1);
    end
    if pick>pcross
        continue;
    end
    shu=randperm(sizepop,2);% 找出交叉个体
    index(1)=shu(1);
    index(2)=shu(2);
    % 选择交叉位置
    pos=ceil(len*rand);
    while pos==1||pos==len
        pos=ceil(len*rand);
    end
    %找出交叉的两个个体
    chrom1=chrom(index(1),:);
    chrom2=chrom(index(2),:);
    k=chrom1(pos:len);
    chrom1(pos:len)=chrom2(pos:len);
    chrom2(pos:len)=k; 
    
    % 满足约束条件赋予新种群
    flag1=test(chrom(index(1),:));
    flag2=test(chrom(index(2),:));
    if flag1*flag2==1
        chrom(index(1),:)=chrom1;
        chrom(index(2),:)=chrom2;
    end   
end
ret=chrom;
end
function ret=BianYi(pm,chrom,np,len)%变异
for i=1:np   
    pick=rand(1);% 变异概率
    while pick==0
        pick=rand(1);
    end
    if pick>pm   % 判断是否变异
        continue;
    end
    index=randperm(np,1);
    pos=randperm(len,1);
    while pos==1
          pos=randperm(len,1);
    end
    nchrom=chrom(index,:);
    nchrom(pos)=randperm(31,1);
    while length(unique(nchrom))==(len-1)
          nchrom(pos)=randperm(31,1);
    end
    flag=test(nchrom);
    if flag==1
        chrom(index,:)=nchrom;
    end
    chrom(index,:)=nchrom;
end
ret=chrom;
end
function newindividuals=Reins(individuals,np,bestindividuals,overbest)%补全种群
m=np+overbest;
newindividuals=struct('fitx',zeros(1,m), 'C',zeros(1,m),'P',zeros(1,m),'x',[]);
% 遗传操作得到的抗体
for i=1:np
    newindividuals.fitx(i)=individuals.fitx(i);   
    newindividuals.C(i)=individuals.C(i);   
    newindividuals.P(i)=individuals.P(i);   
    newindividuals.x(i,:)=individuals.x(i,:);   
end
% 记忆库中抗体
for i=np+1:m
    newindividuals.fitx(i)=bestindividuals.fitx(i-np);   
    newindividuals.C(i)=bestindividuals.C(i-np);   
    newindividuals.P(i)=bestindividuals.P(i-np);   
    newindividuals.x(i,:)=bestindividuals.x(i-np,:);   
end
end
   