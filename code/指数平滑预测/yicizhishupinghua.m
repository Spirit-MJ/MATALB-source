close all
clear
clc
x=[43608	46824	49664	46658	52696	54725	42747	44502	44691	43230	41732	39819];
%x=input('��������ҪԤ���ԭʼ���ݣ��س���������\n');
arf=0.1:0.1:0.9;
if length(x)>=3
   S0=(x(1)+x(2)+x(3))/3;
else
    S0=x(1);
end
for i=2:10
    S(i,1)=S0;
end
for i=1:9
    for j=1:length(x)
        S(i,j+1)=arf(i)*x(j)+(1-arf(i))*S(i,j);
    end
end
for i=1:9
    xsum=0;
    for j=1:length(x)
        xsum=xsum+(S(i,j+1)-x(j))^2;
    end
    B(i)=sqrt(xsum)/length(x);
end
%disp(B);
mn=min(B);
a=find(B==mn);
disp(['���ʺϵ�ƽ��ϵ��Ϊ��',num2str(a/10)]);
disp(['����ֵ���£�',num2str(x)]);
disp(['�ڴ�ƽ��ϵ���µ�Ԥ��ֵ���£�',num2str(S(a,2:length(x)+1))]);
F=(a/10)*x(length(x))+(1-(a/10))*S(a,length(x)+1);
disp(['��һ�ڵ�Ԥ��ֵΪ��',num2str(F)]);
z=1:length(x);
hold on
plot(z,x,'r-*');
plot(z,S(a,2:length(x)+1),'b-o');
hold off
legend('��������','Ԥ������');
title('һ��ָ��ƽ��Ԥ��');