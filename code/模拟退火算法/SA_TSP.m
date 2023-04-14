function main
close all
clear all
clc 
x=[1304 2312;
   3639 1315;
   4177 2244;
   3712 1399;
   3488 1535;
   3326 1556;
   3238 1229;
   4196 1004;
   4312 790;
   4386 570;
   3007 1970;
   2562 1756;
   2788 1491;
   2381 1676;
   1332 695;
   3715 1678;
   3918 2179;
   4061 2370;
   3780 2212;
   3676 2578;
   4029 2838;
   4263 2931;
   3429 1908;
   3507 2367;
   3394 2643;
   3439 3201;
   2935 3240;
   3140 3550;
   2545 2357;
   2778 2826;
   2370 2975];
N=length(x);
for i=1:N
    for j=1:N
        if i==j
        D(i,j)=eps;
        else 
        D(i,j)=sqrt((x(i,1)-x(j,1))^2+(x(i,2)-x(j,2))^2);%�������
        end
    end
end

T0=1e30;%��ʼ�¶�
T_end=1e-30;%��ֹ�¶�
L=10;%���¶��µĵ�������
Q=0.95;%��������
Time=ceil(double(log(T_end/T0)/log(Q)));%��������
count=0;%����
obj=zeros(Time,1);%�����¼
track=zeros(Time,N);%·�߼�¼
for i=1:L
    S11(i,:)=randperm(N);%ÿ���¶��²���L�����н�
    Dd(i)=juli(S11(i,:),x);
end
S1=XuanZe(S11,Dd);%ѡ�����ŵĽ�
%S1=randperm(N);
figure(1)
for i=1:N
S1x(i)=x(S1(i),1);%�������һ����ʼ·��
S1y(i)=x(S1(i),2);%�������һ����ʼ·��
end
S1x=[S1x x(S1(1),1)];%�ѿ�ʼ�ص���ӵ����ͼ
S1y=[S1y x(S1(1),2)];%�ѿ�ʼ�ص���ӵ����ͼ
plot(S1x,S1y,'r-o');%��ʼ���ͼ��
for i=1:N
    text(x(i,1),x(i,2),['   ' num2str(i)]);%��ͼ���ϰ�ÿ������ű���
end
text(x(S1(1),1),x(S1(1),2),'     ���','FontSize',13,'Color',[0 0 1]);%��ͼ���ϱ�����
text(x(S1(end),1),x(S1(end),2),'     �յ�','FontSize',13,'Color',[0 0 1]);%��ͼ���ϱ���յ�
D1=juli(S1,x);
SS1=num2str(S1(1));
for i=2:N
    SS1=[SS1,'�D>',num2str(S1(i))];
end
disp('���������·��Ϊ��');
disp(SS1);
disp(['���������һ��·�߾���Ϊ��',num2str(D1)])
%�����Ż�
while T0>T_end
    count=count+1;
    for i=1:L
        S11(i,:)=randperm(N);%ÿ���¶��²���L�����н�
        Dd(i)=juli(S11(i,:),x);
    end
    S2=XuanZe(S11,Dd);%ѡ�����ŵĽ� 
    a=round(rand(1,2).*(N-1)+1);  %����2��1��31�������
    W=S2(a(1));
    S2(a(1))=S2(a(2));
    S2(a(2))=W;  %����S2(a(1))��S2(a(2))��ֵ���������һ�����е�����λ��
    [S1,R]=Metropolis(S1,S2,D,T0,x);
    if count==1||R<obj(count-1)
       obj(count)=R;
    else
       obj(count)=obj(count-1);
    end
    track(count,:)=S1;
    T0=Q*T0;
    S1=S2;   %��ǰһ��·�ߵ������˴�
    figure(2)
    hold on
    if count>=2
       line([count-1,count],[obj(count-1),obj(count)]);
    end
end
xlabel('��������');
ylabel('����');
title('�Ż�����');
hold off
%plot(1:count,obj,'b')
SS2=num2str(track(end,1));
for i=2:N
    SS2=[SS2,'�D>',num2str(track(end,i))];
end
disp('����·��Ϊ��') 
disp(SS2);
D_end=juli(track(end,:),x);
disp(['�ܾ���Ϊ��',num2str(D_end)]);
figure(3)
for i=1:N
S2x(i)=x(track(end,i),1);
S2y(i)=x(track(end,i),2);
end
S2x=[S2x x(track(end,1),1)];
S2y=[S2y x(track(end,1),2)];
plot(S2x,S2y,'b-o');
for i=1:size(x,1)
    text(x(i,1),x(i,2),['   ' num2str(i)]);
end
text(x(track(end,1),1),x(track(end,1),2),'     ���','FontSize',13,'Color',[1 0 0]);
text(x(track(end,end),1),x(track(end,end),2),'     �յ�','FontSize',13,'Color',[1 0 0]);
function D1=juli(S,X) %���������ľ���
N=length(X);
for i=1:N
S1x(i)=X(S(i),1);
S1y(i)=X(S(i),2);
end
S1x=[S1x X(S(1),1)];
S1y=[S1y X(S(1),2)];
D1=0;
for i=1:N
    D1=D1+sqrt((S1x(i)-S1x(i+1))^2+(S1y(i)-S1y(i+1))^2);
end
end
function [S,R]=Metropolis(S1,S2,D,T,X)
% S1��  ��ǰ��
% S2:   �½�
% D:    ��������������е�֮��ľ��룩
% T:    ��ǰ�¶�
% S��   ��һ����ǰ��
% R��   ��һ����ǰ���·�߾���
R1=juli(S1,X);  %����ǰһ����·�߳���
N=length(X);         %�õ����еĸ���
R2=juli(S2,X);  %����˴ν�·�߳���
dC=R2-R1;   %��������֮��
if dC<0       %����������� ������·��
    S=S2;
    R=R2;
elseif exp(-dC/T)>=rand(1)   %��exp(-dC/T)���ʽ�����·��
    S=S2;
    R=R2;
else        %��������·��
    S=S1;
    R=R1;
end
end
function S1=XuanZe(S11,D)%ѡ��  
       [~,temp]=sort(D);%���վ��뽵������
       S1=S11(temp(1),:);%��ʤ��̭������������ĸ��屣������
end
end
