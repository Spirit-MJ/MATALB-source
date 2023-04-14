close all
clear 
clc
x=[31.45,31.509,31.548,31.09,31.207,31.64,31.635,32.409,31.738];
y=[104.66,104.77,104.618,105.07,105.372,105.172,104.425,104.528,104.749];
w=[33.64,6.09,4.61,9.96,1.56,2.95,3.46,4.46,17.5]*10^8;
a=ones(1,length(x));
for i=1:length(x)
    X(i)=63713930*sind(x(i))*cosd(y(i));
    Y(i)=63713930*sind(x(i))*sind(y(i));
    Z(i)=63713930*cosd(x(i));
end
T1=0;
for i=1:length(x)
   xd1(i)=a(i)*w(i)*X(i);
   T1=T1+xd1(i);
end
T2=0;
for i=1:length(x)
   xd2(i)=a(i)*w(i);
   T2=T2+xd2(i);
end
T3=0;
for i=1:length(x)
   yd1(i)=a(i)*w(i)*Y(i);
   T3=T3+yd1(i);
end
T4=0;
for i=1:length(x)
   yd2(i)=a(i)*w(i)*Z(i);
   T4=T4+yd2(i);
end
T0=0;
for i=1:length(x)
   T0=T0+w(i)*a(i)*sqrt((0-X(i))^2+(0-Y(i))^2+(0-Z(i))^2); 
end
xd(1)=T1/T2;
yd(1)=T3/T2;
zd(1)=T4/T2;
K(1)=T0;
TT0=T0;
TT1=0;
j=0;
while TT0>TT1
    j=j+1;
    T=0; t1=0; t2=0; t3=0; t4=0;
    for i=1:length(x)
      d(i)=sqrt((xd(j)-X(i))^2+(yd(j)-Y(i))^2+(zd(j)-Z(i))^2);
      t1=t1+(a(i)*w(i)*X(i))/d(i);
      t2=t2+(a(i)*w(i))/d(i);
      t3=t3+(a(i)*w(i)*Y(i))/d(i);
      t4=t4+(a(i)*w(i)*Z(i))/d(i);
      xd(j+1)=t1/t2;
      yd(j+1)=t3/t2;
      zd(j+1)=t4/t2;
      T=T+a(i)*w(i)*sqrt((xd(j)-X(i))^2+(yd(j)-Y(i))^2+(zd(j)-Z(i))^2);
      K(j+1)=T;
      TT0=K(j);
      TT1=K(j+1);
    end
end
weidu=acosd(zd(end)/63713930);
SIN=sqrt(1-(zd(end)/63713930)^2);
jingdu=acosd(xd(end)/(63713930*SIN));
format short g
disp('最终选定的坐标为：');
disp([weidu,jingdu]);
disp('最优的费用为：');
disp(K(end));
%绘图
figure(1)
hold on
for i=1:length(x)
    A=[x(i),y(i)];B=[weidu,jingdu];
      c=[A;B];
      fig(i)=plot(c(:,2),c(:,1),'g-')
      set(fig(i),'handlevisibility','off');
end
Fig=scatter(y,x,'k.')
set(Fig,'handlevisibility','off');
scatter(jingdu,weidu,'r*')
legend('选址中心')
text(y(1)-0.1,x(1)-0.05,'涪城区','FontSize',13)
text(y(2)-0.05,x(2)+0.05,'游仙区','FontSize',13)
text(y(3)-0.05,x(3)+0.05,'安州区','FontSize',13)
text(y(4)-0.05,x(4)-0.03,'三台县','FontSize',13)
text(y(5)-0.1,x(5)-0.05,'盐亭县','FontSize',13)
text(y(6)-0.05,x(6)+0.05,'梓潼县','FontSize',13)
text(y(7)-0.03,x(7)+0.05,'北川县','FontSize',13)
text(y(8)-0.03,x(8)+0.05,'平武县','FontSize',13)
text(y(9)-0.03,x(9)+0.05,'江油市','FontSize',13)
title('重心法选址');
xlabel('x坐标');
ylabel('y坐标');

