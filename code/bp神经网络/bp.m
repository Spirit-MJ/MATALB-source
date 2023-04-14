close all
clear
clc
Xtrain1=csvread('C:\Users\idiots\Desktop\IE_Project\工业工程类专业教指委-课程设计展示题目\题目4-轴承故障诊断\降噪后时域频域特征.csv');

Xtrain=[Xtrain1;Xtrain2];
Xtest=[Xtest1;Xtest2];
X=[Xtrain Xtest]';
X0=zscore(X);%标准化
[P,Y,lamda]=pca(X0);
temp=randperm(size(Xtrain,2));%随机产生训练集和测试集
%训练集——130个样本
train_P=Y(temp(1:130),1:10)';
train_T=Xlabel(temp(1:130))';
%测试集——10个样本
test_P=Y(temp(131:end),1:10)';
test_T=Xlabel(temp(131:end))';
N=size(test_P,2);
[train_p,ps_input]=mapminmax(train_P,0,1);%训练集输入数据归一化
test_p=mapminmax('apply',test_P,ps_input);%测试集输入数据归一化
[train_t,ps_output]=mapminmax(train_T,0,1);%训练集输出数据归一化
net=newff(train_p,train_t,10);%,{ 'logsig' 'logsig'},'traingdx');%创建神经网络
%设置训练参数
net.trainParam.epochs=1000;%迭代次数
net.trainParam.goal=1e-30;%误差
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
%结果对比
result=[test_T',sim_T',error']
%绘图
plot(1:N,test_T,'b:*',1:N,sim_T,'r-o')
legend('真实值','预测值')
xlabel('预测样本')
ylabel('预测值')
string={'测试集预测结果对比';['R^2=',num2str(R2)]};
title(string)
shuru_yinhanQuanzhi=net.iw{1,1};
yinhan_shuchuQuanzhi=net.lw{2,1};
yinhancengYuzhi=net.b{1};
shuchucengYuzhi=net.b{2};