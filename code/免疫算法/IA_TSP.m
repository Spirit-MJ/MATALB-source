close all
clear
clc
city=[1304,2312
      3639,1315
      4177,2244
      3712,1399
      3488,1535
      3326,1556
      3238,1229
      4196,1044
      4312,790
      4386,570
      3007,1970
      2562,1756
      2788,1491
      2381,1676
      1332,695
      3715,1678
      3918,2179
      4061,2370
      3780,2212
      3676,2578
      4029,2838
      4263,2931
      3429,1908
      3507,2376
      3394,2643
      3439,3201
      2935,3240
      3140,3550
      2545,2357
      2778,2826
      2370,2975];
np=250;           % 种群规模
overbest=50;          % 记忆库容量
ng=500;            % 迭代次数
pc=0.95;           % 交叉概率
pm=0.2;        % 变异概率
ps=0.9;              % 多样性评价参数
T=0.7;         % 相似度大于阀值，这里阈值为0.7
len=length(city);      %编码长度
M=np+overbest;
D=distanse(city);
%识别抗原,将种群信息定义为一个结构体
individuals=struct('fitx',zeros(1,M), 'C',zeros(1,M),'P',zeros(1,M),'x',[]);
%产生初始抗体群
for i=1:M
    individuals.x(i,:)=randperm(len);%初始化种群
end
trace=[]; %记录每代最个体优适应度和平均适应度
for iii=1:ng
     %抗体群多样性评价
     for i=1:M
         individuals.fitx(i)=fitness(individuals.x(i,:),D);%计算抗体与抗原亲和力(适应度值）
     end
     for i=1:M
         individuals.C(i)=concentration(i,M,individuals,T);%计算抗体浓度
     end
     %综合亲和度和浓度评价抗体优秀程度，得出繁殖概率
     individuals.P=excellence(individuals,M,ps);    
     %记录当代最佳个体和种群平均适应度
     [best,index]=max(individuals.fitx);   % 找出最优适应度 
     bestchrom=individuals.x(index,:);    % 找出最优个体
     average=mean(individuals.fitx);       % 计算平均适应度
     trace=[trace;1/best,1/average];              %记录
     Drawpath(bestchrom,city);
     %根据P，形成父代群，更新记忆库（加入精英保留策略，可由s控制）
     bestindividuals=member(individuals,M,overbest);   % 更新记忆库
     individuals=member(individuals,M,np);      % 形成父代群
     % 选择，交叉，变异操作，再加入记忆库中抗体，产生新种群
     individuals=XuanZe(individuals,np);% 选择
     individuals.x=JiaoCha(individuals.x,pc,bestchrom);% 交叉
     individuals.x=BianYi(individuals.x,pm);% 变异
     individuals.x=Reverse(individuals.x,D);%进化逆转操作
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
        ylabel('距离','fontsize',12)
        title('免疫算法收敛曲线','fontsize',12)
        legend('最优适应度值','平均适应度值')
     end
end
SS1=num2str(bestchrom(1));
for i=2:len
    SS1=[SS1,'―>',num2str(bestchrom(i))];
end
disp(['最优路线为：'])
disp(SS1);
disp(['总距离为：'])
disp(1/best);
function D=distanse(a)%计算距离矩阵
    row=length(a);
    D=zeros(row);
    for i=1:row-1
        for j=i+1:row
            D(i,j)=((a(i,1)-a(j,1))^2+(a(i,2)-a(j,2))^2)^0.5;
            D(j,i)=D(i,j);
        end
        D(i,i)=eps;
    end
end
function len=PathLength(D,x) %计算每个种群中个体的距离
    [NIND,row]=size(x);%种群大小
    len=zeros(NIND,1);
    for j=1:NIND
        d=0;
        for i=1:row-1
            d=d+D(x(j,i),x(j,i+1));
        end
        d=d+D(x(j,end),x(j,1));
        len(j)=d;
    end
end
function fitx=fitness(x,D)
d=0;
for i=1:length(x)-1
     d=d+D(x(i),x(i+1));
end
d=d+D(x(end),x(1));
fitx=1/d;
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
    if find(v(i)==s(i))
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
s=15;
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
function SelCh=JiaoCha(SelCh,Pc,youx)%部分映射交叉
    NSel=size(SelCh,1);%进行选择操作后种群大小
    for i=1:NSel
        if Pc>=rand(1)
            [SelCh(i,:),youx]=intercross(SelCh(i,:),youx);
        end
    end
function [a,b]=intercross(a,b)
    L=length(a);%个体基因的长度，城市的个数
    r1=randsrc(1,2,1:L);%产生2个1:L的随机数
    if r1(1)~=r1(2)
        a0=a;
        b0=b;
        s=min(r1);
        e=max(r1);
        for i=s:e
            a1=a;
            b1=b;
            a(i)=b0(i);%交换个体中的第i个基因
            b(i)=a0(i);%交换个体中的第i个基因
            x=find(a==a(i));%找出交换基因的位置和与交换基因后重复基因的位置
            y=find(b==b(i));%找出交换基因的位置和与交换基因后重复基因的位置
            i1=x(x~=i);%找出交叉后重复基因的位置
            i2=y(y~=i);%找出交叉后重复基因的位置
            if isempty(i1)==0
                a(i1)=a1(i);%利用部分映射消除冲突
            end
            if isempty(i2)==0
                b(i2)=b1(i);%利用部分映射消除冲突
            end
        end
    end
end    
end
function SelCh=BianYi(SelCh,Pm)%变异
    [NSel,L]=size(SelCh);%交叉操作后种群的大小
    for i=1:NSel
        if Pm>=rand
            R0=randperm(L,2);%随机两个基因位置
            R=fliplr(R0);
            SelCh(i,R0)=SelCh(i,R);%交换两个基因的位置
        end
    end
end
function SelCh=Reverse(SelCh,D)%进化逆转
    [row,col]=size(SelCh);%变异操作后种群的大小
    ObjV=PathLength(D,SelCh);%计算变异操作后各个路线的距离
    SelCh1=SelCh;
    for i=1:row
        r=randsrc(1,2,1:col);%产生两个随机数
        mininverse=min(r);
        maxinverse=max(r);
        SelCh1(i,mininverse:maxinverse)=SelCh1(i,maxinverse:-1:mininverse);%逆转操作,基因对换位置
    end
    ObjV1=PathLength(D,SelCh1); %计算路径长度
    index=ObjV1<ObjV;%找出进化逆转操作后路径距离小于进化逆转操作前的索引
    SelCh(index,:)=SelCh1(index,:);%将进化逆转操作后优于原来的个体替换掉
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
function Drawpath(bestchrom,city)%作图
for i=1:length(city)
    xx1(i)=city(bestchrom(i),1);
    yy1(i)=city(bestchrom(i),2);
end
xx=[xx1,city(bestchrom(1),1)];
yy=[yy1,city(bestchrom(1),2)];
figure(2)
plot(xx,yy,'r-*');
for i=1:length(city)
    text(city(i,1),city(i,2),[' ' num2str(i)]);%在图像上把每个点序号标上
end
text(city(bestchrom(1),1),city(bestchrom(1),2),'  起点','FontSize',13,'Color',[0 0 1]);%在图像上标出起点
text(city(bestchrom(end),1),city(bestchrom(end),2),'  终点','FontSize',13,'Color',[0 0 1]);%在图像上标出终点
box on
title('轨迹图');
xlabel('横坐标');
ylabel('纵坐标');
axis('square','equal');
grid on
end  
