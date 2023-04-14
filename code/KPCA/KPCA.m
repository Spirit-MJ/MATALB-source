close all
clear all
clc
Xtrain1=csvread('C:\Users\idiots\Desktop\IE_Project\工业工程类专业教指委-课程设计展示题目\题目4-轴承故障诊断\降噪后时域频域特征.csv');
Xtrain2=csvread('C:\Users\idiots\Desktop\IE_Project\工业工程类专业教指委-课程设计展示题目\题目4-轴承故障诊断\时频特征.csv');Xlabel=csvread('C:\Users\idiots\Desktop\IE_Project\工业工程类专业教指委-课程设计展示题目\题目4-轴承故障诊断\4.label_train.csv');
X=[Xtrain1;Xtrain2];
X=X';
%kpca进行数据提取的函数
rbf_var=10000;
%标准化
X0=zscore(X);
m=size(X0,1);
%计算核矩阵k
for i=1:m
    for j=i:m
        K(i,j)=exp(-norm(X0(i,:)-X0(j,:))^2/rbf_var); %高斯核函数
        K(j,i)=K(i,j);
    end
end
[P,Y,lamda]=pca(K);
per=100*lamda/sum(lamda);
cusmper=zeros(length(per),1);
for i=1:length(lamda)
    cusmper(i)=100*sum(lamda(1:i))/sum(lamda);
end
figure(1)
pareto(per)  %主成分贡献率
xlabel('主成分')
ylabel('贡献率(%)')
title('贡献率')
figure(2)
plot(cusmper,'b-*','LineWidth',1.5)  %累计主成分贡献率
xlabel('累计主成分')
ylabel('累计主成分贡献率(%)')
title('累计贡献率')
text(3+0.2,cusmper(3)-1,[num2str(cusmper(3)),'%'],'FontSize',13)
figure(3)
hold on
grid on
Lei1=find(Xlabel==1);
Lei2=find(Xlabel==2);
Lei3=find(Xlabel==3);
Lei4=find(Xlabel==4);
plot3(Y(Lei1,1),Y(Lei1,2),Y(Lei1,3),'k*')
plot3(Y(Lei2,1),Y(Lei2,2),Y(Lei2,3),'r*')
plot3(Y(Lei3,1),Y(Lei3,2),Y(Lei3,3),'g*')
plot3(Y(Lei4,1),Y(Lei4,2),Y(Lei4,3),'y*')
legend('健康','内圈故障','外圈故障','滚动体故障')
