close all
clear all
clc
load (['D:\idiot file\各种算法\bp神经网络算法\spectra_data.mat']);%60个样本，每个样本401个变量
% IV. 主成分分析
X0=zscore(NIR);
[P0,Y0,lamda0]=pca(X0);
temp=randperm(size(NIR,1));%获取样本数，并随机产生训练级和测试级
%训练集
train_P=Y0(1:50,1:4)';
train_T=octane(1:50)';
%测试级
test_P=Y0(51:60,1:4)';
test_T=octane(51:end)';
N=size(test_P,2);
%数据归一化
[train_p,ps_input]=mapminmax(train_P,0,1);
test_p=mapminmax('apply',test_P,ps_input);
[train_t,ps_output]=mapminmax(train_T,0,1);
R2=0;
%创建神经网络
while R2<=0.99
net=newff(train_p,train_t,7);
%设置训练参数
net.trainParam.epochs=15000;%迭代次数
net.trainParam.goal=1e-3;%误差
net.trainParam.lr=0.01;%学习率
%训练网络
net=train(net,train_p,train_t);
%仿真测试
sim_t=sim(net,test_p);
%数据反归一化
sim_T=mapminmax('reverse',sim_t,ps_output);
%性能评价
error=abs(sim_T-test_T)./test_T;
%决定系数R^2
R2=(N*sum(sim_T.*test_T)-sum(sim_T)*sum(test_T))^2/((N*sum((sim_T).^2)-(sum(sim_T))^2)*(N*sum((test_T).^2)-(sum(test_T))^2));
end
%结果对比
result=[test_T',sim_T',error']
%绘图
plot(1:N,test_T,'b:*',1:N,sim_T,'r-o')
legend('真实值','预测值')
xlabel('预测样本')
ylabel('预测值')
string={'测试集预测结果对比';['R^2=',num2str(R2)]};
title(string)