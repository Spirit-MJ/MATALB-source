close all
clear all
clc
load (['D:\idiot file\各种算法\bp神经网络算法\spectra_data.mat'])%60个样本，每个样本401个变量
R2=0;
while R2<=0.99
temp=randperm(size(NIR,1));%获取样本数，并随机产生训练级和测试级
%训练集
train_P=NIR(temp(1:50),:)';
train_T=octane(temp(1:50),:)';
%测试级
test_P=NIR(temp(51:end),:)';
test_T=octane(temp(51:end),:)';
N=size(test_P,2);
[train_p,ps_input]=mapminmax(train_P,0,1);%训练集输入数据归一化
test_p=mapminmax('apply',test_P,ps_input);%测试集输入数据归一化
[train_t,ps_output]=mapminmax(train_T,0,1);%训练集输出数据归一化
[IW,B,LW,TF,TYPE]=XunLian(train_p,train_t,20,'sig',0);%创建训练
sim_t=YuCe(test_p,IW,B,LW,TF,TYPE);%仿真测试
sim_T=mapminmax('reverse',sim_t,ps_output);%数据反归一化
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
function Y=YuCe(P,IW,B,LW,TF,TYPE)
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