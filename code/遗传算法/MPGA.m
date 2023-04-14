close all
clear all
clc
np=200;%个体数目
mp=200;%种群数目
ng=20;%进化代数
ggap=0.9;%种群代沟
variatenum=2;%变量个数
len=20;%每个变量的编码长度
k0=0;%初始迭代次数
k=0;
variatex=-3:0.01:12.1;
variatey=4.1:0.01:5.8;
for i=1:length(variatex)
    for j=1:length(variatey)
        z(i,j)=21.5+variatex(i)*sin(4*pi*variatex(i))+variatey(j)*sin(20*pi*variatey(j));
    end
end
figure(1)
mesh(variatex',variatey',z')
FieldD=[repmat(len,1,variatenum);[-3,4.1;12.1,5.8];repmat([1;0;1;1],1,variatenum)];%区域描述器
%   FieldD = [len; lb; ub; code; scale; lbin; ubin]
%   len是包含在Chrom中的每个子串的长度，注意sum(len)=size(Chrom, 2)
%   lb和ub分别是每个变量的下界和上界。
%   code指明子串是怎样编码的，1为标准的二进制编码， 0为格雷编码
%   scale指明每个子串所使用的刻度，0表示算术刻度，1表示对数刻度。
%   lbin和ubin指明表示范围中是否包含边界。0表示不包含，1表示包含边界
for i=1:mp
    x{i}=round(rand(np,len*variatenum));%产生初始种群
end
%x{1}(1,:)=zeros(1,len);%初始解
pc=0.7+(0.9-0.7)*rand(mp,1);%在[0.7,0.9]区间内随机产生交叉概率
pm=0.05+(0.25-0.05)*rand(mp,1);%在[0.05,0.25]区间内随机产生变异概率
maxfx=0;
Youfx=zeros(mp,1);%记录精华种群
Youx=zeros(mp,len*variatenum);%记录精华种群的编码
for i=1:mp
    fx{i}=objectfunction(bs2rv(x{i},FieldD));%计算初始种群目标函数值
end
while k0<ng
    k=k+1;
    for j=1:mp
        fitx{j}=fx{j};%计算个种群的适应度
        nx{j}=XuanZe(x{j},ggap,np,fitx{j});%选择操作
        nx{j}=JiaoCha(nx{j},pc(j));%交叉操作
        nx{j}=BianYi(nx{j},pm(j));%变异操作
        objvsel=objectfunction(bs2rv(nx{j},FieldD));%计算子代目标函数值
        [x{j},fx{j}]=Reins(x{j},nx{j},fitx{j},fx{j},objvsel);%补全种群
    end
    [x,fx]=Yimin(x,fx,mp);%移民操作
    [Youfx,Youx]=PeopleSelect(x,fx,Youfx,Youx,mp);%人工选择精华种群
    Zuiyoufx(k)=max(Youfx);%找出精华种群中最优的个体
    if Zuiyoufx(k)>maxfx %判断当前优化值是否比前一次优
        maxfx=Zuiyoufx(k);
        k0=0;
    else
        k0=k0+1;
    end  
    if k>=2
       figure(2)
       line([k-1,k],[Zuiyoufx(k-1),Zuiyoufx(k)]);
       xlabel('进化代数')
       ylabel('最优解变化')
       title('进化过程')
    end
end
function fx=objectfunction(X)
  hang=size(X,1);
  for i=1:hang
      fx(i,1)=21.5+X(i,1)*sin(4*pi*X(i,1))+X(i,2)*sin(20*pi*X(i,2));
  end
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
function [x,objv]=Reins(x,nx,fitx,fx,objvsel)%补全种群
    NIND=size(x,1);%初始种群大小
    NSel=size(nx,1);%选择操作后种群的大小
    [~,index]=sort(fitx,'descend');%按照适应度降序排列
    %[~,index]=sort(fx,'ascend');%升序排列
    x=[nx;x(index(1:NIND-NSel),:)];%将遗传操作前最优的一些个体保留下来
    objv=[fx(index(1:NIND-NSel));objvsel];%子代函数值 
end
function [chrom,objv]=Yimin(chrom,objv,mp)%移民操作
for i=1:mp
    [~,Youi]=max(objv{i});%找出第i个种群中最优的个体
    mubiaoi=i+1;%目标种群(移民操作)
    if mubiaoi>mp
        mubiaoi=1;
    end
    [~,Chai]=min(objv{mubiaoi});%找出目标种群中最劣的个体
    chrom{mubiaoi}(Chai,:)=chrom{i}(Youi,:);%目标种群最劣个体替换为源种群中最劣的个体
    objv{mubiaoi}(Chai)=objv{i}(Youi);
end
end
function [maxobjv,maxchrom]=PeopleSelect(x,fx,maxobjv,maxchrom,mp)%人工选择算子
for i=1:mp
    [Youo,Youi]=max(fx{i});%找出第i种群中最优的个体
    if Youo>maxobjv(i)
       maxobjv(i)=Youo;%记录各种群的精华个体
       maxchrom(i,:)=x{i}(Youi,:);
    end
end
end
