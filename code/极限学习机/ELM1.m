close all
clear all
clc
load(['E:\idiot file\�����㷨\����ѧϰ��\iris_data.mat'])%��������
temp=randperm(size(classes,1));%�������ѵ�����Ͳ��Լ�
% ѵ��������120������
train_P=features(temp(1:120),:)';
train_T=classes(temp(1:120),:)';
% ���Լ�����30������
test_P=features(temp(121:end),:)';
test_T=classes(temp(121:end),:)';
%����ѵ��  
[IW,B,LW,TF,TYPE]=XunLian(train_P,train_T,20,'sig',1);
%�������
T_sim_1=fenlei(train_P,IW,B,LW,TF,TYPE);
T_sim_2=fenlei(test_P,IW,B,LW,TF,TYPE);

result_1=[train_T' T_sim_1'];
result_2=[test_T' T_sim_2'];
%������ȷ��
k1=length(find(train_T==T_sim_1));
n1=length(train_T);
Accuracy_1=k1 / n1 * 100;
disp(['ѵ������ȷ��Accuracy = ' num2str(Accuracy_1) '%(' num2str(k1) '/' num2str(n1) ')'])

k2=length(find(test_T == T_sim_2));
n2=length(test_T);
Accuracy_2 = k2 / n2 * 100;
disp(['���Լ���ȷ��Accuracy = ' num2str(Accuracy_2) '%(' num2str(k2) '/' num2str(n2) ')'])
%��ͼ
figure(2)
plot(1:30,test_T,'bo',1:30,T_sim_2,'r-*')
grid on
xlabel('���Լ��������')
ylabel('���Լ��������')
string = {'���Լ�Ԥ�����Ա�(ELM)';['(��ȷ��Accuracy = ' num2str(Accuracy_2) '%)' ]};
title(string)
legend('��ʵֵ','ELMԤ��ֵ')
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
function Y=fenlei(P,IW,B,LW,TF,TYPE)
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