close all
clear 
clc
warning off
% 注意哦，Matlab中的图节点要从1开始编号，所以这里把0全部改为了9
% 编号最好是从1开始连续编号，不要自己随便定义编号
s = [9 9 1 1 2 2 2 7 7 6 6  5  5 4];
t = [1 7 7 2 8 3 5 8 6 8 5  3  4 3];
w = [4 8 3 8 2 7 4 1 6 6 2 14 10 9];
G = graph(s,t,w);
plot(G, 'EdgeLabel', G.Edges.Weight, 'linewidth', 2) 
set( gca, 'XTick', [], 'YTick', [] );  
[P,d] = shortestpath(G, 9, 4)  %注意：该函数matlab2015b之后才有哦

% 在图中高亮我们的最短路径
myplot = plot(G, 'EdgeLabel', G.Edges.Weight, 'linewidth', 2);  %首先将图赋给一个变量
highlight(myplot, P, 'EdgeColor', 'r')   %对这个变量即我们刚刚绘制的图形进行高亮处理（给边加上r红色）

% 求出任意两点的最短路径矩阵
D = distances(G)   %注意：该函数matlab2015b之后才有哦
D(1,2)  % 1 -> 2的最短路径
D(9,4)  % 9 -> 4的最短路径

% 找出给定范围内的所有点  nearest(G,s,d)
% 返回图形 G 中与节点 s 的距离在 d 之内的所有节点
[nodeIDs,dist] = nearest(G, 2, 10)   %注意：该函数matlab2016a之后才有哦

A=[259 7520.74 86.91 2303.82 84.21 13.32 136.74 1073.15 845.04 84.21];
B=[5947.9 173324 720.7 40678.1 1815.0 348.5 2800 22417.73 18254.5 2946.09];
train_P=[A;B];
for i=1:length(A)
    L(i)=(A(i)/sum(A))/(B(i)/sum(B));
end
