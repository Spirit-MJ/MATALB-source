close all
clear all
clc
%��������
load ('E:\idiot file\�����㷨\֧��������\BreastTissue_data.mat')

%�������ѵ�����Ͳ��Լ�
n=randperm(size(matrix,1));
% ������
% n1=[1:16,22:32,37:51,55:67,71:82,85:101];
% n2=[17:21,33:36,52:54,68:70,83:84,102:106];
%ѵ��������80������
train_matrix=matrix(n(1:80),:);
train_label=label(n(1:80),:);

%���Լ�����26������
test_matrix=matrix(n(81:end),:);
test_label=label(n(81:end),:);

%���ݹ�һ��
[Train_matrix,PS]=mapminmax(train_matrix',0,1);
Train_matrix=Train_matrix';
Test_matrix=mapminmax('apply',test_matrix',PS);
Test_matrix=Test_matrix';

%SVM����/ѵ��(RBF�˺���)
%Ѱ�����c/g��������������֤����
[c,g]=meshgrid(-10:0.2:10,-10:0.2:10);
[m,n]=size(c);
cg=zeros(m,n);
eps=10^(-4);
v=5;
bestc=1;
bestg=0.1;
bestacc=0;
for i=1:m
    for j=1:n
        cmd=['-v ',num2str(v),' -t 2',' -c ',num2str(2^c(i,j)),' -g ',num2str(2^g(i,j))];
        cg(i,j)=libsvmtrain(train_label,Train_matrix,cmd);     
        if cg(i,j)>bestacc
            bestacc=cg(i,j);
            bestc=2^c(i,j);
            bestg=2^g(i,j);
        end        
        if abs(cg(i,j)-bestacc)<=eps&&bestc>2^c(i,j) 
            bestacc=cg(i,j);
            bestc=2^c(i,j);
            bestg=2^g(i,j);
        end               
    end
end
cmd=[' -t 2',' -c ',num2str(bestc),' -g ',num2str(bestg)];

%����/ѵ��SVMģ��
model=libsvmtrain(train_label,Train_matrix,cmd);

%SVM�������
[predict_label_1,accuracy_1,~]=libsvmpredict(train_label,Train_matrix,model);
[predict_label_2,accuracy_2,~]=libsvmpredict(test_label,Test_matrix,model);
disp(['ѵ������׼ȷ��Ϊ��',num2str(accuracy_1(1))]);
disp('ѵ����ԭʼ����Ԥ������')
result_1=[train_label predict_label_1]
disp(['���Լ���׼ȷ��Ϊ��',num2str(accuracy_2(1))]);
disp('���Լ�ԭʼ����Ԥ������')
result_2=[test_label predict_label_2]
%��ͼ
figure
plot(1:length(test_label),test_label,'r-*')
hold on
plot(1:length(test_label),predict_label_2,'b:o')
grid on
legend('��ʵ���','Ԥ�����')
xlabel('���Լ��������')
ylabel('���Լ��������')
string={'���Լ�SVMԤ�����Ա�(RBF�˺���)';
          ['accuracy = ' num2str(accuracy_2(1)) '%']};
title(string)
function [p r]=cross(X,y,fold)
    indices=crossvalind('Kfold',100,fold);      %Indices = crossvalind('Kfold', N, K) K�۽���
    p=0;
    r=0;
    % cp=classperf(X);
    for i=1:fold
        test=(indices==i);
        train=~test;
        train_x=X(train,:);
        train_y=y(train,:);
        test_x=X(test,:);
        test_y=y(test,:);
        b=glmfit(train_x,train_y,'binomial','link','logit');%���߼��ع�������ϵ������ 
        logitFit=glmval(b,test_x,'logit');                  %���߼��ع�Ľ��Ԥ����Լ��Ľ��
        logitFit=(logitFit>=0.5);                           %�������ֵ��С��0.5������Ϊ��1��������0
        pp=0;
        rr=0;
        tpfp=length(find(logitFit==1));                     %Ԥ�����������
        tpfn=length(find(test_y==1));                       %������������
        m=length(test_y);
        tp=0;
        for j=1:m
            if((test_y(j)==logitFit(j))&&(test_y(j)==1))
                tp=tp+1;                                    %Ԥ���������ȷ����������
            end;
        end;
        if(tpfp)
             pp=tp/tpfp;                                    %һ����֤�õ��Ĳ�׼��
        end;
        if(tpfn)
             rr=tp/tpfn;                                    %һ����֤�õ��Ĳ�ȫ��
        end;
        p=p+pp;
        r=r+rr;  
    end
    p=p/fold;                                               %��ƽ��ֵ
    r=r/fold;
end