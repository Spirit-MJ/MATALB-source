close all
clear all
clc
%导入数据
load ('E:\idiot file\各种算法\支持向量机\BreastTissue_data.mat')

%随机产生训练集和测试集
n=randperm(size(matrix,1));
% 留出法
% n1=[1:16,22:32,37:51,55:67,71:82,85:101];
% n2=[17:21,33:36,52:54,68:70,83:84,102:106];
%训练集——80个样本
train_matrix=matrix(n(1:80),:);
train_label=label(n(1:80),:);

%测试集——26个样本
test_matrix=matrix(n(81:end),:);
test_label=label(n(81:end),:);

%数据归一化
[Train_matrix,PS]=mapminmax(train_matrix',0,1);
Train_matrix=Train_matrix';
Test_matrix=mapminmax('apply',test_matrix',PS);
Test_matrix=Test_matrix';

%SVM创建/训练(RBF核函数)
%寻找最佳c/g参数——交叉验证方法
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

%创建/训练SVM模型
model=libsvmtrain(train_label,Train_matrix,cmd);

%SVM仿真测试
[predict_label_1,accuracy_1,~]=libsvmpredict(train_label,Train_matrix,model);
[predict_label_2,accuracy_2,~]=libsvmpredict(test_label,Test_matrix,model);
disp(['训练集的准确率为：',num2str(accuracy_1(1))]);
disp('训练集原始类别和预测的类别：')
result_1=[train_label predict_label_1]
disp(['测试集的准确率为：',num2str(accuracy_2(1))]);
disp('测试集原始类别和预测的类别：')
result_2=[test_label predict_label_2]
%绘图
figure
plot(1:length(test_label),test_label,'r-*')
hold on
plot(1:length(test_label),predict_label_2,'b:o')
grid on
legend('真实类别','预测类别')
xlabel('测试集样本编号')
ylabel('测试集样本类别')
string={'测试集SVM预测结果对比(RBF核函数)';
          ['accuracy = ' num2str(accuracy_2(1)) '%']};
title(string)
function [p r]=cross(X,y,fold)
    indices=crossvalind('Kfold',100,fold);      %Indices = crossvalind('Kfold', N, K) K折交叉
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
        b=glmfit(train_x,train_y,'binomial','link','logit');%用逻辑回归来计算系数矩阵 
        logitFit=glmval(b,test_x,'logit');                  %用逻辑回归的结果预测测试集的结果
        logitFit=(logitFit>=0.5);                           %如果概率值不小于0.5，就认为是1，否则是0
        pp=0;
        rr=0;
        tpfp=length(find(logitFit==1));                     %预测出来的正例
        tpfn=length(find(test_y==1));                       %正例真正个数
        m=length(test_y);
        tp=0;
        for j=1:m
            if((test_y(j)==logitFit(j))&&(test_y(j)==1))
                tp=tp+1;                                    %预测出来的正确的正例个数
            end;
        end;
        if(tpfp)
             pp=tp/tpfp;                                    %一次验证得到的查准率
        end;
        if(tpfn)
             rr=tp/tpfn;                                    %一次验证得到的查全率
        end;
        p=p+pp;
        r=r+rr;  
    end
    p=p/fold;                                               %求平均值
    r=r/fold;
end