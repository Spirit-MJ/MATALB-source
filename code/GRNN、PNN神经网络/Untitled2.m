%清空环境变量
close all
clear
clc
%导入数据
load iris_data.mat
%随机产生训练集和测试集
train_P = [];
train_T = [];
test_P = [];
test_T = [];
for i = 1:3
    temp_input=features((i-1)*50+1:i*50,:);
    temp_output=classes((i-1)*50+1:i*50,:);
    n=randperm(50);
    %训练集——120个样本
    train_P=[train_P temp_input(n(1:40),:)'];
    train_T=[train_T temp_output(n(1:40),:)'];
    %测试集——30个样本
    test_P=[test_P temp_input(n(41:50),:)'];
    test_T=[test_T temp_output(n(41:50),:)'];
end
%模型建立 
train_p=train_P;
test_p=test_P;
%GRNN创建及仿真测试
%创建网络
net_grnn=newgrnn(train_p,train_T);
%仿真测试
t_sim_grnn=sim(net_grnn,test_p);
T_sim_grnn=round(t_sim_grnn);
result_grnn=T_sim_grnn';
%PNN创建及仿真测试
Tc_train=ind2vec(train_T);
%创建网络
net_pnn=newpnn(train_p,Tc_train);
%仿真测试
Tc_test=ind2vec(test_T);
t_sim_pnn=sim(net_pnn,test_p);
T_sim_pnn=vec2ind(t_sim_pnn);
result_pnn=T_sim_pnn';
%性能评价
% 正确率accuracy
accuracy_grnn=length(find(result_grnn==test_T'))/length(test_T);
accuracy_pnn=length(find(result_pnn==test_T'))/length(test_T);
%结果对比
result=[test_T' result_grnn result_pnn]
%绘图
figure(1)
plot(1:30,test_T,'bo',1:30,result_grnn,'r-*',1:30,result_pnn,'k:^')
grid on
xlabel('测试集样本编号')
ylabel('测试集样本类别')
string={'测试集预测结果对比(GRNN vs PNN)';['正确率:' num2str(accuracy_grnn*100) '%(GRNN) vs ' num2str(accuracy_pnn*100) '%(PNN)']};
title(string)
legend('真实值','GRNN预测值','PNN预测值')
