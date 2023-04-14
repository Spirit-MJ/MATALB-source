close all
clear
clc
tic
%训练集/测试集产生
load ('E:\idiot file\各种算法\bp神经网络\spectra_data.mat')
% 随机产生训练集和测试集
temp=randperm(size(NIR,1));
% 训练集——50个样本
train_P=NIR(temp(1:50),:)';
train_T=octane(temp(1:50),:)';
% 测试集——10个样本
test_P=NIR(temp(51:end),:)';
test_T=octane(temp(51:end),:)';
N=size(test_P,2);
%数据归一化
% [train_p,ps_input]=mapminmax(train_P,0,1);%训练集输入数据归一化
% test_p=mapminmax('apply',test_P,ps_input);%测试集输入数据归一化
% [train_t,ps_output]=mapminmax(train_T,0,1);%训练集输出数据归一化
%RBF神经网络创建及仿真测试
%创建网络
net=newrbe(train_P,train_T,0.9);
% 仿真测试
sim_T=sim(net,test_P);
%数据反归一化
% sim_T=mapminmax('reverse',T_sim_rbf,ps_output);
%性能评价相对误差error
error_rbf=abs(sim_T-test_T)./test_T;
% 决定系数R^2
R2_rbf=(N*sum(sim_T.*test_T)-sum(sim_T)*sum(test_T))^2/((N*sum((sim_T).^2)-(sum(sim_T))^2)*(N*sum((test_T).^2)-(sum(test_T))^2));
% 结果对比
result_bp=[test_T',sim_T',error_rbf']
% 绘图
figure
plot(1:N,test_T,'b:*',1:N,sim_T,'k-.^')
legend('真实值','RBF预测值')
xlabel('预测样本')
ylabel('辛烷值')
string={'测试集辛烷值含量预测结果对比(RBF)';['R^2=' num2str(R2_rbf)]};
title(string)
toc
