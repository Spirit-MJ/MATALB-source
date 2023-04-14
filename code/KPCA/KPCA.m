close all
clear all
clc
Xtrain1=csvread('C:\Users\idiots\Desktop\IE_Project\��ҵ������רҵ��ָί-�γ����չʾ��Ŀ\��Ŀ4-��й������\�����ʱ��Ƶ������.csv');
Xtrain2=csvread('C:\Users\idiots\Desktop\IE_Project\��ҵ������רҵ��ָί-�γ����չʾ��Ŀ\��Ŀ4-��й������\ʱƵ����.csv');Xlabel=csvread('C:\Users\idiots\Desktop\IE_Project\��ҵ������רҵ��ָί-�γ����չʾ��Ŀ\��Ŀ4-��й������\4.label_train.csv');
X=[Xtrain1;Xtrain2];
X=X';
%kpca����������ȡ�ĺ���
rbf_var=10000;
%��׼��
X0=zscore(X);
m=size(X0,1);
%����˾���k
for i=1:m
    for j=i:m
        K(i,j)=exp(-norm(X0(i,:)-X0(j,:))^2/rbf_var); %��˹�˺���
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
pareto(per)  %���ɷֹ�����
xlabel('���ɷ�')
ylabel('������(%)')
title('������')
figure(2)
plot(cusmper,'b-*','LineWidth',1.5)  %�ۼ����ɷֹ�����
xlabel('�ۼ����ɷ�')
ylabel('�ۼ����ɷֹ�����(%)')
title('�ۼƹ�����')
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
legend('����','��Ȧ����','��Ȧ����','���������')
