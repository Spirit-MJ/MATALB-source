close all
clear
clc
%除了修改原始样本，还要修改聚类的次数和边界函数
wh = xlsread('E:\桌面\2.xlsx');
x=wh(:,1);
y=wh(:,2);
np=200;%种群大小
pc=0.9;%交叉概率
pm=0.1;%变异概率
ng=200;%进化代数
ggap=0.95;%种群代沟
variatenum=9*2;%变量个数
ca=variatenum/2;%类别数
len=10;%编码长度
yuanshi(x,y,variatenum);%画出不经优化的聚类图
chrom=round(rand(np,variatenum*len));%初始化种群
Field=[repmat(len,1,variatenum);repmat([min(x) min(y);max(x) max(y)],1,variatenum/2);repmat([1;0;0;0],1,variatenum)];%区域描述器
%   FieldD = [len; lb; ub; code; scale; lbin; ubin]
%   len是包含在Chrom中的每个子串的长度，注意sum(len)=size(Chrom, 2)
%   lb和ub分别是每个变量的下界和上界。
%   code指明子串是怎样编码的，1为标准的二进制编码， 0为格雷编码
%   scale指明每个子串所使用的刻度，0表示算术刻度，1表示对数刻度。
%   lbin和ubin指明表示范围中是否包含边界。0表示不包含，1表示包含边界
for k=1:ng
    X=bs2rv(chrom,Field);
    for i=1:np
        fitx(i)=fitness(x,y,X(i,:),variatenum);%计算适应度值
    end
    [Youfitx,index]=max(fitx);
    Youchrom=X(index,:);
    nx=XuanZe(chrom,ggap,np,fitx);%选择操作
    nx=JiaoCha(nx,pc);%交叉操作
    nx=BianYi(nx,pm);%变异操作
    NX=bs2rv(nx,Field);
    for i=1:size(nx,1)
        newfitx(i)=fitness(x,y,NX(i,:),variatenum);%计算适应度值
    end
    [NewYoufitx,Newindex]=max(newfitx);
    NewYouchrom=NX(Newindex,:);
    chrom=Reins(chrom,nx,fitx);%补全种群
    if Youfitx>=NewYoufitx
       trace(1, k)=1/Youfitx-1;
       trace(2, k)=1/mean(fitx)-1;
       Zuiyouchrom=Youchrom;
    else
       trace(1, k)=1/NewYoufitx-1;
       trace(2, k)=1/mean(newfitx)-1;
       Zuiyouchrom=NewYouchrom;
    end
    figure(2)
    if k>=2
       draw1=line([k-1,k],[trace(1, k-1),trace(1, k)]);
       grid on
       set(draw1,'color',[0 0 1]);
       xlabel('进化代数')
       ylabel('准则函数')
       title('进化过程')
    end
end
x0=zeros(1,variatenum/2);
y0=zeros(1,variatenum/2);
for i=1:variatenum/2
    x0(i)=Zuiyouchrom(2*i-1);
    y0(i)=Zuiyouchrom(2*i);
end
num=zeros(1,ca);
for i=1:length(x)
    a=0;
    for j=1:ca
        T(j) = 1 - compute_IoU(x(i),y(i),x0(j),y0(j));
    end
mn=min(T);
a=find(T==mn);
num(a)=num(a)+1;
Kx(a,num(a))=x(i);
Ky(a,num(a))=y(i);
end
for m=1:ca
    x1(m)=mean(Kx(m,1:num(m)));
    y1(m)=mean(Ky(m,1:num(m)));
end 
z0=[x0==x1];
sm0=sum(z0);
z1=[y0==y1];
sm1=sum(z1);
t=0;
while sm0~=ca&&sm1~=ca
    for m=1:ca
        x0(m)=x1(m);
        y0(m)=y1(m);
    end
     num(:)=0;
     for i=1:length(x)
         a=0;
         for j=1:ca
             T(j) = 1 - compute_IoU(x(i),y(i),x0(j),y0(j));
         end
         mn=min(T);
         a=find(T==mn);
         num(a)=num(a)+1;
         Kx(a,num(a))=x(i);
         Ky(a,num(a))=y(i);
     end
      for m=1:ca
          x1(m)=mean(Kx(m,1:num(m)));
          y1(m)=mean(Ky(m,1:num(m)));
          z0=[x0==x1];
          sm0=sum(z0);
          z1=[y0==y1];
          sm1=sum(z1);
      end
      t=t+1;
end
disp(['优化后迭代了',num2str(t),'次']);
figure(3)
hold on
for i=1:ca
    plot(Kx(i,:),Ky(i,:),'*');
end
plot(x0,y0,'kp');
title('遗传算法优化后的聚类')
xlabel('x坐标')
ylabel('y坐标')
axis([min(x) max(x) min(y) max(y)])
hold off
function fitx=fitness(x,y,X,variatenum)
  Kx=zeros(variatenum/2,length(x));
  Ky=zeros(variatenum/2,length(y));
  num=zeros(1,variatenum/2);
for i=1:variatenum/2
    x0(i)=X(2*i-1);
    y0(i)=X(2*i); 
end
for i=1:length(x)
    for j=1:variatenum/2
        T(j) = 1 - compute_IoU(x(i),y(i),x0(j),y0(j));
    end
mn=min(T);
a=find(T==mn);%记录某一类的个数
num(a)=num(a)+1;
Kx(a,num(a))=x(i);
Ky(a,num(a))=y(i);
end
J=0;
for i=1:variatenum/2
      for j=1:num(i)
          J = j + (1 - compute_IoU(Kx(i,j),Ky(i,j),x0(i),y0(i)));
      end
end
fitx=1/(1+J);
end
function nx=XuanZe(x,GGAP,np,fitx)%选择
       GGAP1=floor(GGAP*np);  %选择留下的个体数
       [~,temp]=sort(fitx,'descend');%按照适应度降序排列
       for j=1:GGAP1
           nx(j,:)=x(temp(j),:);%优胜略汰，将适应度大的个体保留下来
       end
end
function nx=JiaoCha(nx,pc)%交叉
    [hang,lie]=size(nx);%进行选择操作后种群大小
    for i=1:2:hang-mod(hang,2)%两两互相交换，防止选择操作后种群数量为单数。
        if pc>=rand(1)
            Rand=randperm(lie,2);%产生两个随机数作为交叉的片段
            Da=max(Rand);%其中大的数
            Xiao=min(Rand);%其中小的数
            temp=nx(i,Xiao:Da);
            nx(i,Xiao:Da)=nx(i+1,Xiao:Da);
            nx(i+1,Xiao:Da)=temp;
        end
    end   
end
function nx=BianYi(nx,pm)%变异
    [hang,lie]=size(nx);%交叉操作后种群的大小
    for i=1:hang
        if pm>=rand
            Rand=randperm(lie,1);
            nx(i,Rand)=~nx(i,Rand);%随机两个基因位置
        end
    end
end
function nx=Reins(x,nx,fitx)%补全种群
    NIND=size(x,1);%初始种群大小
    NSel=size(nx,1);%选择操作后种群的大小
    [~,index]=sort(fitx,'descend');%将遗传操作前的距离排序
    nx=[x(index(1:NIND-NSel),:);nx];%将遗传操作前最优的一些个体保留下来
end
function yuanshi(x,y,variatenum)
ca=variatenum/2;
x0=zeros(1,ca);
y0=zeros(1,ca);
for i=1:ca
    c=round(length(x)*rand(1));
    while c==0
          c=round(length(x)*rand(1));
    end
    x0(i)=x(c);
    y0(i)=y(c);
    if i>1&&x0(i)==x0(i-1)&&y0(i)==y0(i-1)
        i=i-1;
    end
end
num=zeros(1,ca);
for i=1:length(x)
 a=0;
    for j=1:ca
    T(j) = 1 - compute_IoU(x(i),y(i),x0(j),y0(j));
    end
mn=min(T);
a=find(T==mn);
num(a)=num(a)+1;
Kx(a,num(a))=x(i);
Ky(a,num(a))=y(i);
end
for m=1:ca
    x1(m)=mean(Kx(m,1:num(m)));
    y1(m)=mean(Ky(m,1:num(m)));
end 
z0=[x0==x1];
sm0=sum(z0);
z1=[y0==y1];
sm1=sum(z1);
t=0;
while sm0~=ca&&sm1~=ca
    for m=1:ca
    x0(m)=x1(m);
    y0(m)=y1(m);
    end
      num(:)=0;
     for i=1:length(x)
         a=0;
        for j=1:ca
            T(j) = 1 - compute_IoU(x(i),y(i),x0(j),y0(j));
        end
      mn=min(T);
      a=find(T==mn);
      num(a)=num(a)+1;
      Kx(a,num(a))=x(i);
      Ky(a,num(a))=y(i);
     end
      for m=1:ca
    x1(m)=mean(Kx(m,1:num(m)));
    y1(m)=mean(Ky(m,1:num(m)));
    z0=[x0==x1];
    sm0=sum(z0);
    z1=[y0==y1];
    sm1=sum(z1);
      end
      t=t+1;
end
disp(['不经优化迭代了',num2str(t),'次']);
figure(1)
hold on
for i=1:ca
plot(Kx(i,:),Ky(i,:),'*');
end
plot(x0,y0,'kp');
title('不经优化的聚类')
xlabel('x坐标')
ylabel('y坐标')
axis([min(x) max(x) min(y) max(y)])
hold off
end
function IoU = compute_IoU(W1,H1,W2,H2)
w_min = min(W1, W2);
h_min = min(H1, H2);

cross = w_min * h_min;

area1 = W1 * H1;
area2 = W2 * H2;
IoU = cross / (area1 + area2 - cross);
end