close all
clear all
clc
Weight=[71 34 82 23 1 88 12 57 10 68 5 33 37 69 98 24 26 83 16 26 18 43 52 71 22 65 68 8 40 40 24 72 16 34 10 19 28 13 34 98 29 31 79 33 60 74 44 56 54 17];
Value=[26 59 30 19 66 85 94 8 3 44 5 1 41 82 76 1 12 81 73 32 74 54 62 41 19 10 65 53 56 53 70 66 58 22 72 33 96 88 68 45 44 61 78 78 6 66 11 59 83 48];
WeightLimit=300;%限重
np=1000;%种群大小
pc=0.9;%交叉概率
pm=0.3;%变异概率
ng=300;%进化代数
ggap=0.8;%种群代沟
len=length(Weight);%编码长度
x=round(rand(np,len));%初始化种群
x(1,:)=zeros(1,len);%给定初始解
for k=1:ng
    fx=fitness(Weight,Value,x,WeightLimit);
    Objv=fx;
    preObjV(k)=max(Objv);%记录每代中最优
    figure(1)%优化过程图
    if k>=2
       line([k-1,k],[preObjV(k-1),preObjV(k)]);
       xlabel('迭代次数');
       ylabel('价值');
       title('迭代过程');
    end
    nx=XuanZe(x,ggap,np,fx);%选择操作
    nx=JiaoCha(nx,pc);%交叉操作
    nx=BianYi(nx,pm);%变异操作
    nx=Reins(x,nx,Objv);%补全种群
    x=nx;
end
preObjV(ng)
function fx=fitness(Weight,Value,x,WeightLimit)
    for i=1:size(x,1)
        fx(i)=sum(x(i,:).*Value);
        if sum(x(i,:).*Weight)>=WeightLimit
           fx(i)=0; 
        end
    end
end
function nx=XuanZe(x,GGAP,np,fx)%选择
       GGAP1=floor(GGAP*np);  %选择留下的个体数
%        cumfx=cumsum(fx)./sum(fx);%轮盘赌
%        for j=1:GGAP1
%            sita=rand(1);
%            for i=1:length(fx)
%                if cumfx(i)==0
%                   nx(j,:)=x(randperm(np,1),:);
%                else sita<=cumfx(i)
%                    nx(j,:)=x(i,:);
%                    break;
%                end
%            end
%        end
       [~,temp]=sort(fx,'descend');%按照适应度降序排列
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
function nx=Reins(x,nx,Objv)%补全种群
    NIND=size(x,1);%初始种群大小
    NSel=size(nx,1);%选择操作后种群的大小
    [~,index]=sort(Objv,'descend');%将遗传操作前的距离排序
    nx=[x(index(1:NIND-NSel),:);nx];%将遗传操作前最优的一些个体保留下来
end