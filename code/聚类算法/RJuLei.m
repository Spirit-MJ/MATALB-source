close all
clear 
clc
OriginX1=xlsread('C:\Users\idiots\Desktop\����ȷ��');
OriginX2=xlsread('C:\Users\idiots\Desktop\��������');
OriginX3=xlsread('C:\Users\idiots\Desktop\��������');
OriginX1_1=cumsum(OriginX1);
OriginX1_1_daysure=zeros(1,size(OriginX1_1,2));
for i=1:size(OriginX1_1,2)
    OriginX1_1_daysure(i)=sum(OriginX1_1(:,i)~=0);
end
OriginX2_1=cumsum(OriginX2);
OriginX2_1_daydeath=zeros(1,size(OriginX2_1,2));
for i=1:size(OriginX2_1,2)
    OriginX2_1_daydeath(i)=sum(OriginX2_1(:,i)~=0);
end
X1=sum(OriginX1)./OriginX1_1_daysure;%ƽ��ÿ��ȷ������
X7_1=cumsum(OriginX2)./cumsum(OriginX1);
X7_1(ismissing(X7_1))=0;
X2=sum(X7_1)./OriginX2_1_daydeath;%ƽ��������
X3=[OriginX1_1_daysure(1:4) find(OriginX1(:,5)==max(OriginX1(:,5)))-(size(OriginX1,1)-OriginX1_1_daysure(5))... 
    find(OriginX1(:,6)==max(OriginX1(:,6)))-(size(OriginX1,1)-OriginX1_1_daysure(6)) OriginX1_1_daysure(7)...
    find(OriginX1(:,8)==max(OriginX1(:,8)))-(size(OriginX1,1)-OriginX1_1_daysure(8))...
    find(OriginX1(:,9)==max(OriginX1(:,9)))-(size(OriginX1,1)-OriginX1_1_daysure(9)) OriginX1_1_daysure(10)...
    find(OriginX1(:,11)==max(OriginX1(:,11)))-(size(OriginX1,1)-OriginX1_1_daysure(11))...
    find(OriginX1(:,12)==max(OriginX1(:,12)))-(size(OriginX1,1)-OriginX1_1_daysure(12))...
    find(OriginX1(:,13)==max(OriginX1(:,13)))-(size(OriginX1,1)-OriginX1_1_daysure(13))...
    find(OriginX1(:,14)==max(OriginX1(:,14)))-(size(OriginX1,1)-OriginX1_1_daysure(14))...
    find(OriginX1(:,15)==max(OriginX1(:,15)))-(size(OriginX1,1)-OriginX1_1_daysure(15))]; %�����鷢չ�������ۼ�ȷ�����߳��ֹյ��ʱ�䣬δ���ֹյ㰴��ֹ������
X4=sum(OriginX3)./sum(OriginX1);%��ֹ6��2���ۼ���������ռ�ۼ�ȷ�������ٷֱ�
X5=sum(OriginX1);%��ֹ6��2���ۼ�ȷ������
X6=max(OriginX1);%�����鷢չ��������ȷ��������ֵ
X7=sum(OriginX2)./sum(OriginX1);%��ֹ6��2���ۼ���������ռ�ۼ�ȷ�������ٷֱ�
X8=max(X7_1);%���������
X9=sum(OriginX2)./OriginX2_1_daydeath;%ƽ��ÿ����������
X=[X1;X2;X3;X4]';
XX=zscore(X);
r=zeros(length(XX));
for i=1:length(XX)
    r(i,i)=1;
    for j=i+1:length(XX)
        r(i,j)=sum(XX(i,:).*XX(j,:))/(sqrt(sum(XX(i,:).^2))*sqrt(sum(XX(j,:).^2)));
        r(j,i)=r(i,j);
    end
end
a=r;
l=input('������Ҫ�ֳɼ��ࣺ');
d=1-abs(a);
z=linkage(d,'ward');%������뷨���� ?
% ? ? ? ?�ַ��� ? ? ? ? �� ?��?
% ? ? ��single�� ? ? ? ��̾��루ȱʡ��?
% ? ? ��complete�� ? ? ?������?
% ? ? ��average�� ? ? ?ƽ������?
% ? ? ��centroid�� ? ? ?���ľ���?
% ? ? ��ward�� ? ? ?���ƽ���ͷ�����Ward������
L=cluster(z,'maxclust',l);
for i=1:l
b=find(L==i);
b=reshape(b,1,length(b));
disp(['��',num2str(i),'����У�',num2str(b)]);
end
dendrogram(z);
xlabel('���');
ylabel('�������');
title('R�;���');