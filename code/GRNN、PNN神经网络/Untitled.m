%清空环境变量
close all
clear
clc

%导入数据
load iris_data.mat
%随机产生训练集和测试集
train_P=[];
train_T=[];
test_P=[];
test_T=[];
for i=1:3
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
feature=size(train_P,1);
%模型建立 
result_grnn = [];
result_pnn = [];
for i=1:feature
    for j=i:feature
        p_train=train_P(i:j,:);
        p_test=test_P(i:j,:);
        %GRNN创建及仿真测试
        %创建网络
        net_grnn=newgrnn(p_train,train_T);
        %仿真测试
        t_sim_grnn=sim(net_grnn,p_test);
        T_sim_grnn=round(t_sim_grnn);
        result_grnn=[result_grnn T_sim_grnn'];
        %PNN创建及仿真测试
        Tc_train=ind2vec(train_T);
        %创建网络
        net_pnn=newpnn(p_train,Tc_train);
        %仿真测试
        Tc_test=ind2vec(test_T);
        t_sim_pnn=sim(net_pnn,p_test);
        T_sim_pnn=vec2ind(t_sim_pnn);
        result_pnn=[result_pnn T_sim_pnn'];
    end
end
%性能评价
% 正确率accuracy
accuracy_grnn = [];
accuracy_pnn = [];
for i=1:size(result_grnn,2)
    accuracy_1=length(find(result_grnn(:,i) == test_T'))/length(test_T);
    accuracy_2=length(find(result_pnn(:,i) == test_T'))/length(test_T);
    accuracy_grnn=[accuracy_grnn accuracy_1];
    accuracy_pnn=[accuracy_pnn accuracy_2];
end
%结果对比
result=[test_T' result_grnn result_pnn]
accuracy=[accuracy_grnn;accuracy_pnn]
%绘图
figure(1)
plot(1:30,test_T,'bo',1:30,result_grnn(:,4),'r-*',1:30,result_pnn(:,4),'k:^')
grid on
xlabel('测试集样本编号')
ylabel('测试集样本类别')
string = {'测试集预测结果对比(GRNN vs PNN)';['正确率:' num2str(accuracy_grnn(4)*100) '%(GRNN) vs ' num2str(accuracy_pnn(4)*100) '%(PNN)']};
title(string)
legend('真实值','GRNN预测值','PNN预测值')
figure(2)
plot(1:size(result_grnn,2),accuracy(1,:),'r-*',1:size(result_grnn,2),accuracy(2,:),'b:o')
grid on
xlabel('模型编号')
ylabel('测试集正确率')
title('10个模型的测试集正确率对比(GRNN vs PNN)')
legend('GRNN','PNN')

