close all
clear
clc
x=1999:2018;
y=[590.7 622.9 675.1 741.1 768.3 837.3 916.6 968.7 1059.1 1578.66 1590.52 1807.88 2016.17 2238.28 2248.60 2465.13 2387.44 2504.11 2696.17 2946.09];
R=corrcoef(x,y);
if R(1,2)>0.8
    disp('变量x和y强相关，相关系数为：');
else
    disp('变量x和y非强相关，相关系数为：');
end
disp(R(1,2));
xsum=0;ysum=0;
for i=1:length(x)
    xsum=xsum+x(i);
end
x1=xsum/length(x);
for i=1:length(y)
    ysum=ysum+y(i);
end
y1=ysum/length(y);
b1=0;b2=0;
for i=1:length(x)
    b1=b1+(x(i)-x1)*(y(i)-y1);
    b2=b2+(x(i)-x1)^2;
    b=b1/b2;
end
a=y1-b*x1;
z=1999:2018;
F=a+b*z;
t=input('请输入想要预测的年数：');
Y=a+b*t;
disp([num2str(t),'年的预测数据为：',num2str(Y)]);
hold on
plot(x,y,'*');
plot(z,F,'r-');
hold off



