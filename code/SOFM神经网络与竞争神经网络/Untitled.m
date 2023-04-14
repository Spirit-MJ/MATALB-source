close all
clear 
clc
%训练集/测试集产生
%导入数据
load water_data.mat
%数据归一化
attributes=mapminmax(attributes,0,1);
%训练集——35个样本
train_P=attributes(:,1:35);
train_T=classes(:,1:35);
%测试集——4个样本
test_P=attributes(:,36:end);
test_T=classes(:,36:end);

%竞争神经网络创建、训练及仿真测试
%创建网络
net=newc(minmax(train_P),4,0.01,0.001);
%net=newc(PR,S,KLR,CLR)
%S:神经元的数目
%KLR：Kohonen学习速度，默认为0.01
%CLR：Conscience学习速度，默认为0.001
%设置训练参数
net.trainParam.epochs=500;
%训练网络
net=train(net,train_P);
%仿真测试
%训练集
t_sim_compet_1=sim(net,train_P);
T_sim_compet_1=vec2ind(t_sim_compet_1);
%测试集
t_sim_compet_2=sim(net,test_P);
T_sim_compet_2=vec2ind(t_sim_compet_2);

%SOFM神经网络创建、训练及仿真测试
%创建网络
net=newsom(train_P,[4 4]);
%设置训练参数
net.trainParam.epochs=200;
%训练网络
net=train(net,train_P);
%仿真测试
%训练集
t_sim_sofm_1=sim(net,train_P);
T_sim_sofm_1=vec2ind(t_sim_sofm_1);
%测试集
t_sim_sofm_2=sim(net,test_P);
T_sim_sofm_2=vec2ind(t_sim_sofm_2);

%结果对比
% 竞争神经网络
result_compet_1 = [train_T' T_sim_compet_1']
result_compet_2 = [test_T' T_sim_compet_2']
% SOFM神经网络
result_sofm_1 = [train_T' T_sim_sofm_1']
result_sofm_2 = [test_T' T_sim_sofm_2']
