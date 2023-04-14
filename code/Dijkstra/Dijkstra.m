close all
clear 
clc
warning off
% ע��Ŷ��Matlab�е�ͼ�ڵ�Ҫ��1��ʼ��ţ����������0ȫ����Ϊ��9
% �������Ǵ�1��ʼ������ţ���Ҫ�Լ���㶨����
s = [9 9 1 1 2 2 2 7 7 6 6  5  5 4];
t = [1 7 7 2 8 3 5 8 6 8 5  3  4 3];
w = [4 8 3 8 2 7 4 1 6 6 2 14 10 9];
G = graph(s,t,w);
plot(G, 'EdgeLabel', G.Edges.Weight, 'linewidth', 2) 
set( gca, 'XTick', [], 'YTick', [] );  
[P,d] = shortestpath(G, 9, 4)  %ע�⣺�ú���matlab2015b֮�����Ŷ

% ��ͼ�и������ǵ����·��
myplot = plot(G, 'EdgeLabel', G.Edges.Weight, 'linewidth', 2);  %���Ƚ�ͼ����һ������
highlight(myplot, P, 'EdgeColor', 'r')   %��������������Ǹոջ��Ƶ�ͼ�ν��и����������߼���r��ɫ��

% ���������������·������
D = distances(G)   %ע�⣺�ú���matlab2015b֮�����Ŷ
D(1,2)  % 1 -> 2�����·��
D(9,4)  % 9 -> 4�����·��

% �ҳ�������Χ�ڵ����е�  nearest(G,s,d)
% ����ͼ�� G ����ڵ� s �ľ����� d ֮�ڵ����нڵ�
[nodeIDs,dist] = nearest(G, 2, 10)   %ע�⣺�ú���matlab2016a֮�����Ŷ

A=[259 7520.74 86.91 2303.82 84.21 13.32 136.74 1073.15 845.04 84.21];
B=[5947.9 173324 720.7 40678.1 1815.0 348.5 2800 22417.73 18254.5 2946.09];
train_P=[A;B];
for i=1:length(A)
    L(i)=(A(i)/sum(A))/(B(i)/sum(B));
end
