close all
clear
clc
x=[12452.5 8812.6 9763.1 12308.3 15479.4 17475.8 9714.5 10302.6 10898.1 11952.2];
%x=input('������Ԥ����Ҫ�����ݣ��س�������\n');
arf=0.1:0.1:0.9;
if length(x)>=3
    S10=(x(1)+x(2)+x(3))/3;
else
    S10=x(1);
end
for i=1:9
    S1(i,1)=S10;
    S2(i,1)=S10;
    S3(i,1)=S10;
end
for i=1:9
    for j=1:length(x)
        S1(i,j+1)=arf(i)*x(j)+(1-arf(i))*S1(i,j);
        S2(i,j+1)=arf(i)*S1(i,j+1)+(1-arf(i))*S2(i,j);
        S3(i,j+1)=arf(i)*S2(i,j+1)+(1-arf(i))*S3(i,j);
    end
end
for i=1:9
    xsum=0;
    for j=1:length(x)
        xsum=xsum+(S3(i,j+1)-x(j))^2;
    end
    B(i)=sqrt(xsum)/length(x);
end
mn=min(B);
m=find(B==mn);
a0=m/10;
disp(['���ʺϵ�ƽ��ϵ��Ϊ��',num2str(a0)]);
disp(['�ڴ�ƽ��ϵ���µ�Ԥ��ֵ���£�',num2str(S3(m,2:length(x)+1))]);
a=3*S1(m,length(x)+1)-3*S2(m,length(x)+1)+S3(m,length(x)+1);
b=a0/(2*(1-a0)^2)*((6-5*a0)*S1(m,length(x)+1)-2*(5-4*a0)*S2(m,length(x)+1)+(4-3*a0)*S3(m,length(x)+1));
c=a0^2/(2*(1-a0)^2)*(S1(m,length(x)+1)-2*S2(m,length(x)+1)+S3(m,length(x)+1));
t=input('��������ҪԤ���¼��꣺');
Yt=a+b*t+c*t^2;
disp(['��',num2str(t),'��Ԥ��ֵΪ��',num2str(Yt)]);
x0=1:length(x);
hold on
plot(x0,x,'b--d');
plot(x0,S3(m,2:length(x)+1),'r-*');
hold off
legend('��������','Ԥ������');
title('����ָ��ƽ��Ԥ��');
