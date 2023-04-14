close all
clear 
clc
%ѵ����/���Լ�����
%��������
load water_data.mat
%���ݹ�һ��
attributes=mapminmax(attributes,0,1);
%ѵ��������35������
train_P=attributes(:,1:35);
train_T=classes(:,1:35);
%���Լ�����4������
test_P=attributes(:,36:end);
test_T=classes(:,36:end);

%���������紴����ѵ�����������
%��������
net=newc(minmax(train_P),4,0.01,0.001);
%net=newc(PR,S,KLR,CLR)
%S:��Ԫ����Ŀ
%KLR��Kohonenѧϰ�ٶȣ�Ĭ��Ϊ0.01
%CLR��Conscienceѧϰ�ٶȣ�Ĭ��Ϊ0.001
%����ѵ������
net.trainParam.epochs=500;
%ѵ������
net=train(net,train_P);
%�������
%ѵ����
t_sim_compet_1=sim(net,train_P);
T_sim_compet_1=vec2ind(t_sim_compet_1);
%���Լ�
t_sim_compet_2=sim(net,test_P);
T_sim_compet_2=vec2ind(t_sim_compet_2);

%SOFM�����紴����ѵ�����������
%��������
net=newsom(train_P,[4 4]);
%����ѵ������
net.trainParam.epochs=200;
%ѵ������
net=train(net,train_P);
%�������
%ѵ����
t_sim_sofm_1=sim(net,train_P);
T_sim_sofm_1=vec2ind(t_sim_sofm_1);
%���Լ�
t_sim_sofm_2=sim(net,test_P);
T_sim_sofm_2=vec2ind(t_sim_sofm_2);

%����Ա�
% ����������
result_compet_1 = [train_T' T_sim_compet_1']
result_compet_2 = [test_T' T_sim_compet_2']
% SOFM������
result_sofm_1 = [train_T' T_sim_sofm_1']
result_sofm_2 = [test_T' T_sim_sofm_2']
