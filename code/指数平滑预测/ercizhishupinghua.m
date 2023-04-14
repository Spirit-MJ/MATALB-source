close all
clear
clc
% x=[2981,3106,3236,3542,3830,4369,3675,3799,4916,6393,5856,6920,7117,8025,9280,9633,9694,9374];
%x=[6611.8,6976.6,8108.3,7068.3,8089.6,8145.9,8174.9,8438.5,8233.3,7200.3,7775.9,8398.2,8810.9,10045.9,8727];
x=[594.24,692.89,727.03,752.02,737.82,746.27,690.49,628.59,671.91];
%x=input('������Ԥ����Ҫ��ԭʼ���ݣ��س�������\n');
arf=0.1:0.1:0.9;
if length(x)>=3
    S10=(x(1)+x(2)+x(3))/3;
else
    S10=x(1);
end
for i=1:9
    S1(i,1)=S10;
    S2(i,1)=arf(i)*x(1)+(1-arf(i))*S10;
end
for i=1:9
    for j=1:length(x)
        S1(i,j+1)=arf(i)*x(j)+(1-arf(i))*S1(i,j);
        S2(i,j+1)=arf(i)*S1(i,j+1)+(1-arf(i))*S2(i,j);
    end
end
for i=1:9
    xsum=0;
    for j=1:length(x)
       xsum=xsum+(S2(i,j+1)-x(j))^2;
    end
    B(i)=sqrt(xsum)/length(x);
end
mn=min(B);
m=find(B==mn);
a0=m/10;
disp(['���ʺϵ�ƽ��ϵ��Ϊ��',num2str(a0)]);
disp(['����ֵ���£�',num2str(x)]);
disp(['�ڴ�ƽ��ϵ���µ�Ԥ��ֵ���£�',num2str(S2(m,2:length(x)+1))]);
a=2*S1(m,length(x)+1)-S2(m,length(x)+1);
b=(a0/(1-a0))*(S1(m,length(x)+1)-S2(m,length(x)+1));
t=input('��������ҪԤ�������:');
Yt=a+b*t;
disp(['��',num2str(t),'���ڵ�Ԥ��ֵΪ��',num2str(Yt)]);
x0=1:length(x);
hold on
plot(x0+2010,x,'b--d');
plot(x0+2010,S2(m,2:length(x)+1),'r-*');
hold off
legend('��������','Ԥ������');
xlabel('���')
ylabel('���(��Ԫ)')
title('����ָ��ƽ��Ԥ��');

