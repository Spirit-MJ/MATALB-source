close all
clear all
clc
x=1:10;
y=[4375 5131 7939 7595 7677 7316 7779 8311 8844 9240];
x0=0:0.1:10;
A1=polyfit(x,y,1);
y1=polyval(A1,x0);
A2=polyfit(x,y,2);
y2=polyval(A2,x0);
A3=polyfit(x,y,3);
y3=polyval(A3,x0);
A6=polyfit(x,y,6);
y6=polyval(A6,x0);
hold on
scatter(x,y,'*');
plot(x0,y1,'--');
plot(x0,y2,'');
plot(x0,y3,'+');
plot(x0,y6,'x');
axis('square')
hold off
legend('原始数据','一维拟合','二维拟合','三维拟合','六维拟合');
P1=poly2str(A1,'x');
P2=poly2str(A2,'x');
P3=poly2str(A3,'x');
P6=poly2str(A6,'x');

