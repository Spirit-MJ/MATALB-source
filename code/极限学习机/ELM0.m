close all
clear all
clc
load (['D:\idiot file\�����㷨\bp�������㷨\spectra_data.mat'])%60��������ÿ������401������
R2=0;
while R2<=0.99
temp=randperm(size(NIR,1));%��ȡ�����������������ѵ�����Ͳ��Լ�
%ѵ����
train_P=NIR(temp(1:50),:)';
train_T=octane(temp(1:50),:)';
%���Լ�
test_P=NIR(temp(51:end),:)';
test_T=octane(temp(51:end),:)';
N=size(test_P,2);
[train_p,ps_input]=mapminmax(train_P,0,1);%ѵ�����������ݹ�һ��
test_p=mapminmax('apply',test_P,ps_input);%���Լ��������ݹ�һ��
[train_t,ps_output]=mapminmax(train_T,0,1);%ѵ����������ݹ�һ��
[IW,B,LW,TF,TYPE]=XunLian(train_p,train_t,20,'sig',0);%����ѵ��
sim_t=YuCe(test_p,IW,B,LW,TF,TYPE);%�������
sim_T=mapminmax('reverse',sim_t,ps_output);%���ݷ���һ��
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
function [IW,B,LW,TF,TYPE]=XunLian(P,T,N,TF,TYPE)
% P  �������
% T  �������
% N  ��������Ԫ����
% TF ���ݺ�������:
%       ��sig' for Sigmoidal function(Ĭ��)
%       ��sin' for Sine function
%       ��hardlim' for Hardlim function
% TYPE �����������(0���ع顣1������)
% IW  �������������㵽������֮���Ȩֵ
% B   ��������Ԫ����ֵ
% LW  �����������֮���Ȩֵ
[R,Q]=size(P);
if TYPE== 1
    T=ind2vec(T);
end
[S,Q] = size(T);
IW=rand(N,R)*2-1;%�������������㵽������֮���Ȩֵ
B=rand(N,1);%���������������Ԫ����ֵ
BiasMatrix=repmat(B,1,Q);
tempH=IW*P+BiasMatrix;%�õ���������Ԫ������
switch TF
    case 'sig'
        H=1./(1+exp(-tempH));
    case 'sin'
        H=sin(tempH);
    case 'hardlim'
        H=hardlim(tempH);
end
LW=pinv(H')*T';%ͨ������ó������������֮���Ȩֵ
end
function Y=YuCe(P,IW,B,LW,TF,TYPE)
% P   �������
% IW  ����㵽������֮���Ȩֵ
% B   ��������Ԫ����ֵ
% LW  �����������֮���Ȩֵ
% TF  ���ݺ�������:
%      ��sig' for Sigmoidal function (default)
%      ��sin' for Sine function
%      ��hardlim' for Hardlim function
% TYPE �����������(0���ع顣1������)
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