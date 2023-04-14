close all
clear
clc
x_origin=[487 426 344 413 579 651 1290 1722 2089 2981 3106 3236 3542 ...
   3830 4369 3675 3799 4916 6393 5856 6920 7693 7117 8025 9280 ...
   9633 9694 9374];

lag=3;%自回归阶数
x0=zeros(length(x_origin)-lag,lag); 
for i=1:length(x_origin)-lag
    x0(i,:) = x_origin(i:i+lag-1);
end
y0=x_origin(lag+1:end)';

split = 0.8;%训练集比例
temp0=floor(length(y0)*split);%训练集个数
temp=randperm(length(y0));%随机产生训练集和测试集

%训练集
train_P=x0(temp(1:temp0),:)';
train_T=y0(temp(1:temp0),:)';
%测试级
test_P=x0(temp(temp0+1:end),:)';
test_T=y0(temp(temp0+1:end),:)';
N=size(test_P,2);
%数据归一化
[train_p,ps_input]=mapminmax(train_P,0,1);
test_p=mapminmax('apply',test_P,ps_input);
[train_t,ps_output]=mapminmax(train_T,0,1);
%创建神经网络
nero = 2 * lag - 1; %神经元个数
net=newff(train_p,train_t,nero);
%设置训练参数
net.trainParam.epochs=15000;%迭代次数
net.trainParam.goal=1e-5;%误差
net.trainParam.lr=0.03;%学习率
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

% 测试
% xx=[1.75,1.20];
% sim_t=sim(net,xx)
% sim_T=mapminmax('reverse',sim_t,ps_output)