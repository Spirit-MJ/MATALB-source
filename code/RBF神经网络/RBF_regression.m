close all
clear
clc
tic
%ѵ����/���Լ�����
load ('E:\idiot file\�����㷨\bp������\spectra_data.mat')
% �������ѵ�����Ͳ��Լ�
temp=randperm(size(NIR,1));
% ѵ��������50������
train_P=NIR(temp(1:50),:)';
train_T=octane(temp(1:50),:)';
% ���Լ�����10������
test_P=NIR(temp(51:end),:)';
test_T=octane(temp(51:end),:)';
N=size(test_P,2);
%���ݹ�һ��
% [train_p,ps_input]=mapminmax(train_P,0,1);%ѵ�����������ݹ�һ��
% test_p=mapminmax('apply',test_P,ps_input);%���Լ��������ݹ�һ��
% [train_t,ps_output]=mapminmax(train_T,0,1);%ѵ����������ݹ�һ��
%RBF�����紴�����������
%��������
net=newrbe(train_P,train_T,0.9);
% �������
sim_T=sim(net,test_P);
%���ݷ���һ��
% sim_T=mapminmax('reverse',T_sim_rbf,ps_output);
%��������������error
error_rbf=abs(sim_T-test_T)./test_T;
% ����ϵ��R^2
R2_rbf=(N*sum(sim_T.*test_T)-sum(sim_T)*sum(test_T))^2/((N*sum((sim_T).^2)-(sum(sim_T))^2)*(N*sum((test_T).^2)-(sum(test_T))^2));
% ����Ա�
result_bp=[test_T',sim_T',error_rbf']
% ��ͼ
figure
plot(1:N,test_T,'b:*',1:N,sim_T,'k-.^')
legend('��ʵֵ','RBFԤ��ֵ')
xlabel('Ԥ������')
ylabel('����ֵ')
string={'���Լ�����ֵ����Ԥ�����Ա�(RBF)';['R^2=' num2str(R2_rbf)]};
title(string)
toc
