close all
clear
clc
tic
x=-1.5:0.01:1.5;
y=-2:0.01:2;
for i=1:length(y)
    for j=1:length(x)
        z(i,j)=sin(sqrt((y(i))^2+x((j))^2))/sqrt((y(i))^2+(x(j))^2)+exp((cos(2*pi*x(j))+cos(2*pi*y(i)))/2)-2.71289;
    end
end
figure(1)
mesh(x,y,z);
grid on
%初始化参数
c1=1.49445;
c2=1.49445;
maxgen=100;   % 迭代次数  
sizepop=100;   %种群规模
Vmax=0.5;  %函数自变量取值  
Vmin=-0.5;
popmax=2;%函数值范围
popmin=-2;
%产生初始粒子和速度
for i=1:sizepop
    % 随机产生一个种群
    pop(i,:)=2*rands(1,2);  %初始种群，随机产生-2——2的数
    V(i,:)=0.5*rands(1,2);  %初始化速度，随机产生-0.5—0.5之间的数
    %计算适应度
    fitx(i)=fitness(pop(i,:));  %计算函数值 
end
%个体极值和群体极值
[bestfitness,bestindex]=max(fitx);%找出适应度最大的
Pbest=pop;   %个体最优
Gbest=pop(bestindex,:);   %全局最佳，找到最优x
fitnessPbest=fitx;   %个体最佳适应度值，此代所有适应度值(Y值)
fitnessGbest=bestfitness;   %全局最佳适应度值，此代最优的适应度值
% VI. 迭代寻优
for i=1:maxgen
    for j=1:sizepop
        % 速度更新
        V(j,:)=V(j,:)+c1*rand(1)*(Pbest(j,:)-pop(j,:))+c2*rand(1)*(Gbest-pop(j,:));
        V(j,find(V(j,:)>Vmax))=Vmax;
        V(j,find(V(j,:)<Vmin))=Vmin;
        
        % 种群更新
        pop(j,:)=pop(j,:)+V(j,:);
        pop(j,find(pop(j,:)>popmax))=popmax;
        pop(j,find(pop(j,:)<popmin))=popmin;
        % 适应度值更新
        fitx(j)=fitness(pop(j,:)); 
    end
    for j=1:sizepop    
        % 个体最优更新
        if fitx(j)>fitnessPbest(j) %如果更新后的种群中有某个粒子适应度大于前一次对应的适应度值
            Pbest(j,:)=pop(j,:); %就将这个粒子(x的值)赋值给gbest(此代对应的x值)保留下来 
            fitnessPbest(j)=fitx(j);%更新此代对应的适应度值(Y值)
        end
        % 群体最优更新
        if fitx(j)>fitnessGbest %如果更新后的种群中有某个粒子适应度大于前一次最优的适应度值
            Gbest=pop(j,:);%将最优的这个粒子(x的值)记录下来赋值给zbest
            fitnessGbest=fitx(j);%将最优的适应度值(y的值)记录下来赋值给fitnesszbest
        end
        fx(i)=fitnessGbest;
    end  
    figure(2)
    if i>1
        line([i-1,i],[fx(i-1),fx(i)]);
    end
    title('最优个体适应度','fontsize',12);
    xlabel('进化代数','fontsize',12);
    ylabel('适应度','fontsize',12);
end
% VII. 输出结果并绘图
disp(['函数的最大值为：',num2str(fx(end))])
figure(3)
hold on
mesh(x,y,z);
plot3(Gbest(1),Gbest(2),fx(end),'r*')
toc
function y=fitness(x)
y=sin(sqrt((x(1))^2+(x(2))^2))/sqrt((x(1))^2+(x(2))^2)+exp((cos(2*pi*x(1))+cos(2*pi*x(2)))/2)-2.71289;
end


