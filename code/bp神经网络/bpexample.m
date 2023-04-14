close all
clear
clc
x_origin=[487 426 344 413 579 651 1290 1722 2089 2981 3106 3236 3542 ...
   3830 4369 3675 3799 4916 6393 5856 6920 7693 7117 8025 9280 ...
   9633 9694 9374];

lag=3;%�Իع����
x0=zeros(length(x_origin)-lag,lag); 
for i=1:length(x_origin)-lag
    x0(i,:) = x_origin(i:i+lag-1);
end
y0=x_origin(lag+1:end)';

split = 0.8;%ѵ��������
temp0=floor(length(y0)*split);%ѵ��������
temp=randperm(length(y0));%�������ѵ�����Ͳ��Լ�

%ѵ����
train_P=x0(temp(1:temp0),:)';
train_T=y0(temp(1:temp0),:)';
%���Լ�
test_P=x0(temp(temp0+1:end),:)';
test_T=y0(temp(temp0+1:end),:)';
N=size(test_P,2);
%���ݹ�һ��
[train_p,ps_input]=mapminmax(train_P,0,1);
test_p=mapminmax('apply',test_P,ps_input);
[train_t,ps_output]=mapminmax(train_T,0,1);
%����������
nero = 2 * lag - 1; %��Ԫ����
net=newff(train_p,train_t,nero);
%����ѵ������
net.trainParam.epochs=15000;%��������
net.trainParam.goal=1e-5;%���
net.trainParam.lr=0.03;%ѧϰ��
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
%����Ա�
result=[test_T',sim_T',error']
%��ͼ
plot(1:N,test_T,'b:*',1:N,sim_T,'r-o')
legend('��ʵֵ','Ԥ��ֵ')
xlabel('Ԥ������')
ylabel('Ԥ��ֵ')
string={'���Լ�Ԥ�����Ա�';['R^2=',num2str(R2)]};
title(string)

% ����
% xx=[1.75,1.20];
% sim_t=sim(net,xx)
% sim_T=mapminmax('reverse',sim_t,ps_output)