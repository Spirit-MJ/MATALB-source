close all
clear 
clc
OriginX1=xlsread('C:\Users\idiots\Desktop\新增确诊');
OriginX2=xlsread('C:\Users\idiots\Desktop\新增死亡');
OriginX3=xlsread('C:\Users\idiots\Desktop\新增治愈');
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
X1=sum(OriginX1)./OriginX1_1_daysure;%平均每天确诊人数
X7_1=cumsum(OriginX2)./cumsum(OriginX1);
X7_1(ismissing(X7_1))=0;
X2=sum(X7_1)./OriginX2_1_daydeath;%平均致死率
X3=[OriginX1_1_daysure(1:4) find(OriginX1(:,5)==max(OriginX1(:,5)))-(size(OriginX1,1)-OriginX1_1_daysure(5))... 
    find(OriginX1(:,6)==max(OriginX1(:,6)))-(size(OriginX1,1)-OriginX1_1_daysure(6)) OriginX1_1_daysure(7)...
    find(OriginX1(:,8)==max(OriginX1(:,8)))-(size(OriginX1,1)-OriginX1_1_daysure(8))...
    find(OriginX1(:,9)==max(OriginX1(:,9)))-(size(OriginX1,1)-OriginX1_1_daysure(9)) OriginX1_1_daysure(10)...
    find(OriginX1(:,11)==max(OriginX1(:,11)))-(size(OriginX1,1)-OriginX1_1_daysure(11))...
    find(OriginX1(:,12)==max(OriginX1(:,12)))-(size(OriginX1,1)-OriginX1_1_daysure(12))...
    find(OriginX1(:,13)==max(OriginX1(:,13)))-(size(OriginX1,1)-OriginX1_1_daysure(13))...
    find(OriginX1(:,14)==max(OriginX1(:,14)))-(size(OriginX1,1)-OriginX1_1_daysure(14))...
    find(OriginX1(:,15)==max(OriginX1(:,15)))-(size(OriginX1,1)-OriginX1_1_daysure(15))]; %从疫情发展以来到累计确诊曲线出现拐点的时间，未出现拐点按截止日期算
X4=sum(OriginX3)./sum(OriginX1);%截止6月2日累计治愈人数占累计确诊人数百分比
X5=sum(OriginX1);%截止6月2日累计确诊人数
X6=max(OriginX1);%从疫情发展以来新增确诊人数极值
X7=sum(OriginX2)./sum(OriginX1);%截止6月2日累计死亡人数占累计确诊人数百分比
X8=max(X7_1);%最高致死率
X9=sum(OriginX2)./OriginX2_1_daydeath;%平均每天死亡人数
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
l=input('请输入要分成几类：');
d=1-abs(a);
z=linkage(d,'ward');%按最长距离法聚类 ?
% ? ? ? ?字符串 ? ? ? ? 含 ?义?
% ? ? ’single’ ? ? ? 最短距离（缺省）?
% ? ? ’complete’ ? ? ?最大距离?
% ? ? ’average’ ? ? ?平均距离?
% ? ? ’centroid’ ? ? ?重心距离?
% ? ? ’ward’ ? ? ?离差平方和方法（Ward方法）
L=cluster(z,'maxclust',l);
for i=1:l
b=find(L==i);
b=reshape(b,1,length(b));
disp(['第',num2str(i),'类的有：',num2str(b)]);
end
dendrogram(z);
xlabel('类别');
ylabel('分类过程');
title('R型聚类');