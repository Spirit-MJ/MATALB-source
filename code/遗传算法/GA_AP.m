close all
clear all
clc
C=[6 7 11 2
   4 5 9  8
   3 1 10 4
   5 9 8  2]
np=1000;%种群大小
pc=0.9;%交叉概率
pm=0.1;%变异概率
ng=100;%进化代数
ggap=0.8;%种群代沟
len=length(C)*length(C);%编码长度
x=round(rand(np,len));%初始化种群
x(1,:)=[1 zeros(1,4) 1 zeros(1,4) 1 zeros(1,4) 1];%给定初始解
for k=1:ng
    fx=fitness(C,x);
    Objv=fx;
    [preObjV(k),index]=max(Objv);%记录每代中最优
    Youx=x(index,:);
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
You=reshape(Youx,length(C),length(C))'
preObjV(k)
function fx=fitness(C,x)
    for i=1:size(x,1)
        x0=reshape(x(i,:),length(C),length(C))';   
        if sum(x0,1)==ones(1,size(C,2))&sum(x0,2)==ones(size(C,1),1)
           fx(i)=1/(sum(sum(x0.*C)));
        else
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