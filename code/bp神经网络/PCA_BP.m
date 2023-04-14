close all
clear all
clc
load (['D:\idiot file\�����㷨\bp�������㷨\spectra_data.mat']);%60��������ÿ������401������
% IV. ���ɷַ���
X0=zscore(NIR);
[P0,Y0,lamda0]=pca(X0);
temp=randperm(size(NIR,1));%��ȡ�����������������ѵ�����Ͳ��Լ�
%ѵ����
train_P=Y0(1:50,1:4)';
train_T=octane(1:50)';
%���Լ�
test_P=Y0(51:60,1:4)';
test_T=octane(51:end)';
N=size(test_P,2);
%���ݹ�һ��
[train_p,ps_input]=mapminmax(train_P,0,1);
test_p=mapminmax('apply',test_P,ps_input);
[train_t,ps_output]=mapminmax(train_T,0,1);
R2=0;
%����������
while R2<=0.99
net=newff(train_p,train_t,7);
%����ѵ������
net.trainParam.epochs=15000;%��������
net.trainParam.goal=1e-3;%���
net.trainParam.lr=0.01;%ѧϰ��
%ѵ������
net=train(net,train_p,train_t);
%�������
sim_t=sim(net,test_p);
%���ݷ���һ��
sim_T=mapminmax('reverse',sim_t,ps_output);
%��������
error=abs(sim_T-test_T)./test_T;
%����ϵ��R^2
R2=(N*sum(sim_T.*test_T)-sum(sim_T)*sum(test_T))^2/((N*sum((sim_T).^2)-(sum(sim_T))^2)*(N*sum((test_T).^2)-(sum(test_T))^2));
end
%����Ա�
result=[test_T',sim_T',error']
%��ͼ
plot(1:N,test_T,'b:*',1:N,sim_T,'r-o')
legend('��ʵֵ','Ԥ��ֵ')
xlabel('Ԥ������')
ylabel('Ԥ��ֵ')
string={'���Լ�Ԥ�����Ա�';['R^2=',num2str(R2)]};
title(string)