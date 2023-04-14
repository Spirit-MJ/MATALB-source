close all
clear all
clc
load(['E:\idiot file\各种算法\极限学习机\iris_data.mat'])%导入数据
temp=randperm(size(classes,1));%随机产生训练集和测试集
% 训练集——120个样本
train_P=features(temp(1:120),:)';
train_T=classes(temp(1:120),:)';
% 测试集——30个样本
test_P=features(temp(121:end),:)';
test_T=classes(temp(121:end),:)';
%进行训练  
[IW,B,LW,TF,TYPE]=XunLian(train_P,train_T,20,'sig',1);
%仿真测试
T_sim_1=fenlei(train_P,IW,B,LW,TF,TYPE);
T_sim_2=fenlei(test_P,IW,B,LW,TF,TYPE);

result_1=[train_T' T_sim_1'];
result_2=[test_T' T_sim_2'];
%计算正确率
k1=length(find(train_T==T_sim_1));
n1=length(train_T);
Accuracy_1=k1 / n1 * 100;
disp(['训练集正确率Accuracy = ' num2str(Accuracy_1) '%(' num2str(k1) '/' num2str(n1) ')'])

k2=length(find(test_T == T_sim_2));
n2=length(test_T);
Accuracy_2 = k2 / n2 * 100;
disp(['测试集正确率Accuracy = ' num2str(Accuracy_2) '%(' num2str(k2) '/' num2str(n2) ')'])
%绘图
figure(2)
plot(1:30,test_T,'bo',1:30,T_sim_2,'r-*')
grid on
xlabel('测试集样本编号')
ylabel('测试集样本类别')
string = {'测试集预测结果对比(ELM)';['(正确率Accuracy = ' num2str(Accuracy_2) '%)' ]};
title(string)
legend('真实值','ELM预测值')
function [IW,B,LW,TF,TYPE]=XunLian(P,T,N,TF,TYPE)
% P  输入矩阵
% T  输出矩阵
% N  隐含层神经元个数
% TF 传递函数类型:
%       ‘sig' for Sigmoidal function(默认)
%       ‘sin' for Sine function
%       ‘hardlim' for Hardlim function
% TYPE 解决问题类型(0：回归。1：分类)
% IW  随机产生的输入层到隐含层之间的权值
% B   隐含层神经元的阈值
% LW  隐含层输出层之间的权值
[R,Q]=size(P);
if TYPE== 1
    T=ind2vec(T);
end
[S,Q] = size(T);
IW=rand(N,R)*2-1;%随机产生的输入层到隐含层之间的权值
B=rand(N,1);%随机产生隐含层神经元的阈值
BiasMatrix=repmat(B,1,Q);
tempH=IW*P+BiasMatrix;%得到隐含层神经元的输入
switch TF
    case 'sig'
        H=1./(1+exp(-tempH));
    case 'sin'
        H=sin(tempH);
    case 'hardlim'
        H=hardlim(tempH);
end
LW=pinv(H')*T';%通过计算得出隐含层输出层之间的权值
end
function Y=fenlei(P,IW,B,LW,TF,TYPE)
% P   输入矩阵
% IW  输入层到隐含层之间的权值
% B   隐含层神经元的阈值
% LW  隐含层输出层之间的权值
% TF  传递函数类型:
%      ‘sig' for Sigmoidal function (default)
%      ‘sin' for Sine function
%      ‘hardlim' for Hardlim function
% TYPE 解决问题类型(0：回归。1：分类)
Q=size(P,2);
BiasMatrix=repmat(B,1,Q);
tempH=IW*P+BiasMatrix;
switch TF
    case 'sig'
        H=1./(1+exp(-tempH));
    case 'sin'
        H=sin(tempH);
    case 'hardlim'
        H=hardlim(tempH);
end
Y=(H'*LW)';
if TYPE==1
    temp_Y=zeros(size(Y));
    for i=1:size(Y,2)
        [~,index] = max(Y(:,i));
        temp_Y(index,i)=1;
    end
    Y=vec2ind(temp_Y); 
end
end