close all
clear
clc
tic  %记录时间
%参数设置
fishnum=100; %生成100只人工鱼
MAXGEN=50; %最多迭代次数
try_number=100;%最多试探次数
visual=1; %感知距离
delta=0.618; %拥挤度因子
step=0.1; %步长
%初始化鱼群
lb_ub=[-10 -10;10 10];
X=AF_init(fishnum,lb_ub);
gen=1;
BestY=-1*ones(1,MAXGEN); %每步中最优的函数值
besty=-100; %最优函数值
Y=AF_foodconsistence(X);
while gen<=MAXGEN
    for i=1:fishnum
        [Xi1,Yi1]=AF_swarm(X,i,visual,step,delta,try_number,lb_ub,Y); %聚群行为
        [Xi2,Yi2]=AF_follow(X,i,visual,step,delta,try_number,lb_ub,Y); %追尾行为
        if Yi1>Yi2
            X(:,i)=Xi1;
            Y(1,i)=Yi1;
        else
            X(:,i)=Xi2;
            Y(1,i)=Yi2;
        end
    end
    [Ymax,index]=max(Y);
    if Ymax>besty
        besty=Ymax;
        bestx=X(:,index);
        BestY(gen)=Ymax;
        bestX(:,gen)=X(:,index);
    else
        BestY(gen)=BestY(gen-1);
        bestX(:,gen)=bestX(:,gen-1);
    end
    if gen>=2
       figure(1)
       line([gen-1,gen],[BestY(gen-1),BestY(gen)])
       xlabel('迭代次数')
       ylabel('优化值')
       title('鱼群算法迭代过程')
    end
    gen=gen+1;
end
figure(2)
hold on
for i=1:MAXGEN
    plot(bestX(1,i),bestX(2,i),'.')
end
plot(bestx(1),bestx(2),'ro','MarkerSize',20)
xlabel('x')
ylabel('y')
title('鱼群算法迭代过程中最优坐标')


disp(['最优解X：',num2str(bestx','%1.5f')])
disp(['最优解Y：',num2str(besty,'%1.5f')])
toc
function X=AF_init(Nfish,lb_ub)
row=size(lb_ub,2);  %变量个数
    for j=1:row
        lb(j)=lb_ub(1,j);
        ub(j)=lb_ub(2,j);
        X(j,:)=lb(j)+(ub(j)-lb(j))*rand(1,Nfish);
    end
end
function [Y]=AF_foodconsistence(X)
fishnum=size(X,2);
for i=1:fishnum
    Y(i)=sin(X(1,i))/X(1,i)*sin(X(2,i))/X(2,i);   
end
end
function [Xnext,Ynext]=AF_prey(Xi,ii,visual,step,try_number,LBUB,lastY)%觅食行为
Xnext=[];
Yi=lastY(ii);
for i=1:try_number
    Xj=Xi+(2*rand(length(Xi),1)-1)*visual;
    Yj=AF_foodconsistence(Xj);
    if Yi<Yj
        Xnext=Xi+rand*step*(Xj-Xi)/norm(Xj-Xi);
        for i=1:length(Xnext)
            if  Xnext(i)>LBUB(2,i)
                Xnext(i)=LBUB(2,i);
            end
            if  Xnext(i)<LBUB(1,i)
                Xnext(i)=LBUB(1,i);
            end
        end
        Xi=Xnext;
        break;
    end
end
%随机行为
if isempty(Xnext)
    Xj=Xi+(2*rand(length(Xi),1)-1)*visual;
    Xnext=Xj;
    for i=1:length(Xnext)
        if  Xnext(i)>LBUB(2,i)
            Xnext(i)=LBUB(2,i);
        end
        if  Xnext(i)<LBUB(1,i)
            Xnext(i)=LBUB(1,i);
        end
    end
end
Ynext=AF_foodconsistence(Xnext);
end
function [Xnext,Ynext]=AF_swarm(X,i,visual,step,deta,try_number,LBUB,lastY)% 聚群行为
Xi=X(:,i);
D=AF_dist(Xi,X);
index=find(D>0&D<visual);
nf=length(index);
if nf>0
    for j=1:size(X,1)
        Xc(j,1)=mean(X(j,index));
    end
    Yc=AF_foodconsistence(Xc);
    Yi=lastY(i);
    if Yc/nf>deta*Yi
        Xnext=Xi+rand*step*(Xc-Xi)/norm(Xc-Xi);
        for i=1:length(Xnext)
            if  Xnext(i)>LBUB(2,i)
                Xnext(i)=LBUB(2,i);
            end
            if  Xnext(i)<LBUB(1,i)
                Xnext(i)=LBUB(1,i);
            end
        end
        Ynext=AF_foodconsistence(Xnext);
    else
        [Xnext,Ynext]=AF_prey(Xi,i,visual,step,try_number,LBUB,lastY);
    end
else
    [Xnext,Ynext]=AF_prey(Xi,i,visual,step,try_number,LBUB,lastY);
end
end
function D=AF_dist(Xi,X)
col=size(X,2);
D=zeros(1,col);
for j=1:col
    D(j)=norm(Xi-X(:,j));
end
end
function [Xnext,Ynext]=AF_follow(X,i,visual,step,deta,try_number,LBUB,lastY)% 追尾行为
Xi=X(:,i);
D=AF_dist(Xi,X);
index=find(D>0 & D<visual);
nf=length(index);
if nf>0
    XX=X(:,index);
    YY=lastY(index);
    [Ymax,Max_index]=max(YY);
    Xmax=XX(:,Max_index);
    Yi=lastY(i);
    if Ymax/nf>deta*Yi
        Xnext=Xi+rand*step*(Xmax-Xi)/norm(Xmax-Xi);
        for i=1:length(Xnext)
            if  Xnext(i)>LBUB(2,i)
                Xnext(i)=LBUB(2,i);
            end
            if  Xnext(i)<LBUB(1,i)
                Xnext(i)=LBUB(1,i);
            end
        end
        Ynext=AF_foodconsistence(Xnext);
    else
        [Xnext,Ynext]=AF_prey(X(:,i),i,visual,step,try_number,LBUB,lastY);
    end
else
    [Xnext,Ynext]=AF_prey(X(:,i),i,visual,step,try_number,LBUB,lastY);
end
end

