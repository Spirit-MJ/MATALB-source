close all
clear
clc

figure(1)
plot([1,200],[1,200],'.'); % 规划空间
hold on
xlabel('km','fontsize',12)
ylabel('km','fontsize',12)
title('二维规划空间','fontsize',12)

% 描述起点和终点
Start = [20,180];
End = [160,90];
plot([Start(1),End(1)],[Start(2),End(2)],'.');

% 图形标注
text(Start(1)+2,Start(2),'起点');  % 起点
text(End(1)+2,End(2),'终点');  % 终点

np = 1000;  % 种群大小
pc = 0.9;  % 交叉概率
pm = 0.2;  % 变异概率
ng = 200;  % 进化代数
GGAP = 0.95;  % 种群代沟
view = 9;  % 可视域大小
num_gen = floor((End(1)-Start(1)-1)/view);  % 每个染色体上基因个数
zuiyou_fit = zeros(ng,1);
zuiyou_path = zeros(ng,num_gen); 

x = zeros(np,num_gen);
for i=1:np
    x(i,:) = ceil(90+90*rand(1,num_gen));  % 初始化种群实数编码
end

for k = 1:ng
    fx = fitness(x,Start,End,view);  % 计算适应度
    [zuiyou,index] = max(fx);
    zuiyou_fit(k) = 1/zuiyou;
    zuiyou_path(k,:) = x(index,:);
    nx = XuanZe(x,GGAP,np,fx);  % 选择操作
    nx = JiaoCha(nx,pc);  % 交叉操作
    nx = BianYi(nx,pm);  % 变异
    x = Reins(x,nx,fx);  % 补全种群
    figure(2)
    if k>=2
        line([k-1,k],[zuiyou_fit(k-1),zuiyou_fit(k)]);
        xlabel('遗传代数');
        ylabel('距离');
        title('进化过程');
        grid on
    end
end
figure(1)
hold on
loc_y = [Start(2),zuiyou_path(end,:),End(2)];
loc_x = Start(1);
for i = 1:length(zuiyou_path(end,:))
    tenp = Start(1) + view * i;
    loc_x = [loc_x,tenp];
end
loc_x = [loc_x,End(1)];
plot(loc_x,loc_y,'r-.'); % 规划空间 
disp(['最优长为',num2str(zuiyou_fit(end)),'km'])

function fit = fitness(x0,Start,End,view)  % 计算适应度
fit = zeros(size(x0,1),1);
for j = 1:size(x0,1)
    fit(j) = sqrt(view^2 + (Start(2) - x0(j,1))^2);
    for i = 1:size(x0,2)-1
        fit(j) = fit(j) + sqrt(view^2 + (x0(j,i) - x0(j,i+1))^2);
    end
    fit(j) = 1/(fit(j) + sqrt(view^2 + (End(2) - x0(j,end))^2));
end

end

function nx = XuanZe(x0,GGAP,np,fx)  % 选择算子
       GGAP1 = floor(GGAP*np);  %选择留下的个体数      
       [~,temp] = sort(fx,'descend');  % 按照适应度降序排列
       for j = 1:GGAP1
           nx(j,:) = x0(temp(j),:);  % 优胜略汰，将适应度大的个体保留下来
       end
end

function x0 = JiaoCha(x0,Pc)
 [NSel,L] = size(x0);  % 进行选择操作后种群大小
    for i = 1:NSel/2
        if Pc>=rand(1)
            rand_2 = randperm(L,2);  % 随机选择交叉的位置
            temp = x0(2*i-1,min(rand_2):max(rand_2));
            x0(2*i-1,min(rand_2):max(rand_2)) = x0(2*i,min(rand_2):max(rand_2));  % 交叉
            x0(2*i,min(rand_2):max(rand_2)) = temp;  % 交叉
        end
    end
end

function x0 = BianYi(x0,Pm)  % 变异
    [NSel,L] = size(x0);  % 交叉操作后种群的大小
    for i = 1:NSel
        if Pm>=rand
            R0 = randperm(L,1);  % 随机选择变异的基因位置
            R1 = ceil(90+90*rand(1));  % 随机选择变异的基因
            x0(i,R0) = R1;  % 变异
        end
    end
end
 
function Chrom = Reins(Chrom,SelCh,ObjV)  % 补全种群
    NIND = size(Chrom,1);  % 初始种群大小
    NSel = size(SelCh,1);  % 选择操作后种群的大小
    [~,index] = sort(ObjV,'descend');  % 将遗传操作前的距离排序
    size(Chrom(index(1:NIND-NSel),:));
    Chrom = [SelCh;Chrom(index(1:NIND-NSel),:)];  % 将遗传操作前最优的一些个体保留下来
end

