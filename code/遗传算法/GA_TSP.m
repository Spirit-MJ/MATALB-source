close all
clear
clc
x0=[1304 2312
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
   2370 2975];
np=2000;%种群大小
len=length(x0);
pc=0.95;%交叉概率
pm=0.05;%变异概率
ng=500;%进化代数
GGAP=0.9;%种群代沟
D=distanse(x0);%距离矩阵
for i=1:np
    x(i,:)=randperm(len);%初始化种群 实数编码
end

for k=1:ng
ObjV=PathLength(D,x);%计算随机种群中个体的距离
[preObjV(k),index]=min(ObjV);%记录每代中最优距离
youx=x(index,:);
%youx(1,:)=[6,7,13,12,14,15,1,29,31,30,27,28,26,25,20,21,22,18,3,17,19,24,11,23,16,4,8,9,10,2,5];
figure(1)%优化过程图
if k>=2
   line([k-1,k],[preObjV(k-1),preObjV(k)]);
   xlabel('遗传代数');
   ylabel('距离');
   title('进化过程');
   grid on
end
fx=1./ObjV;  %计算种群中每个个体的适应度
Drawpath(fx,x,x0);
nx=XuanZe(x,GGAP,np,fx); %选择操作
nx=JiaoCha(nx,pc,youx);%交叉操作
nx=BianYi(nx,pm);%变异操作
nx=Reverse(nx,D);%进化逆转操作
x=Reins(x,nx,ObjV);%子代
end
a=find(fx==max(fx));%找出最终适应度最大个体的索引号
a3=a(1);%找出最终适应度最大个体的索引号
SS1=num2str(x(a3,1));
for i=2:len
    SS1=[SS1,'―>',num2str(x(a3,i))];
end
disp(['最优路线为：'])
disp(SS1);
disp(['总距离为：'])
disp(ObjV(a3));
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
function len=PathLength(D,Chrom) %计算每个种群中个体的距离
    row=size(D,2);%编码长度
    NIND=size(Chrom,1);%种群大小
    len=zeros(NIND,1);
    for i=1:NIND
        p=[Chrom(i,:) Chrom(i,1)];
        i1=p(1:end-1);
        i2=p(2:end);
        for j=1:row
            len(i)=len(i)+D(i1(j),i2(j));%每个个体距离
        end
    end
end
function nx=XuanZe(x,GGAP,np,fx)%选择
       GGAP1=floor(GGAP*np);  %选择留下的个体数      
       [~,temp]=sort(fx,'descend');%按照适应度降序排列
       for j=1:GGAP1
           nx(j,:)=x(temp(j),:);%优胜略汰，将适应度大的个体保留下来
       end
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
function Chrom=Reins(Chrom,SelCh,ObjV)%补全种群
    NIND=size(Chrom,1);%初始种群大小
    NSel=size(SelCh,1);%选择操作后种群的大小
    [~,index]=sort(ObjV);%将遗传操作前的距离排序
    Chrom=[Chrom(index(1:NIND-NSel),:);SelCh];%将遗传操作前最优的一些个体保留下来
end
function Drawpath(fx,x,x0)%作图
a=find(fx==max(fx));%找出最终适应度最大个体的索引号
a3=a(1);%找出最终适应度最大个体的索引号
for i=1:length(x0)
    xx1(i)=x0(x(a3,i),1);
    yy1(i)=x0(x(a3,i),2);
end
xx=[xx1,x0(x(a3,1),1)];
yy=[yy1,x0(x(a3,1),2)];
figure(2)
plot(xx,yy,'r-*');
for i=1:length(x0)
    text(x0(i,1),x0(i,2),[' ' num2str(i)]);%在图像上把每个点序号标上
end
text(x0(x(a3,1),1),x0(x(a3,1),2),'  起点','FontSize',13,'Color',[0 0 1]);%在图像上标出起点
text(x0(x(a3,end),1),x0(x(a3,end),2),'  终点','FontSize',13,'Color',[0 0 1]);%在图像上标出终点
box on
title('轨迹图');
xlabel('横坐标');
ylabel('纵坐标');
axis('square','equal');
grid on
end
