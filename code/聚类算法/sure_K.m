close all
clear all
clc
aa=rand(1,300);
bb=rand(1,300);
cc=rand(1,300)+10;
dd=rand(1,300)+10;
ee=rand(1,400)-10;
ff=rand(1,400)-10;
x=[aa cc ee];
y=[bb dd ff];
% x=rand(1,1000);
% y=rand(1,1000);
% Data=xlsread('C:\Users\idiots\Desktop\1.xlsx');
% x=Data(:,2);
% y=Data(:,1);
for ca=1:9
    [x0,y0,Kx,Ky,num]=yuanshi(x,y,ca);
    Dpath=0;
    %Rpath=0;
    for i=ca
        Dpath(i)=leizuD(Kx(i,1:num(i)),Ky(i,1:num(i)));
        %Rpath(i)=leizuR(x0(i),y0(i),Kx(i,1:num(i)),Ky(i,1:num(i)));
    end
    meansDpath(ca)=mean(Dpath);
    %meansRpath(ca)=mean(Rpath);
end
plot(1:ca,meansDpath)
xlabel('聚类中心的个数')
ylabel('类族平均直径/半径')
%legend('类族平均直径','类族平均半径')
function [x0,y0,Kx,Ky,num]=yuanshi(x,y,ca)
x0=zeros(1,ca);
y0=zeros(1,ca);
for i=1:ca
    c=round(length(x)*rand(1));
    while c==0
          c=round(length(x)*rand(1));
    end
    x0(i)=x(c);
    y0(i)=y(c);
    if i>1&&x0(i)==x0(i-1)&&y0(i)==y0(i-1)
        i=i-1;
    end
end
num=zeros(1,ca);
for i=1:length(x)
 a=0;
    for j=1:ca
    T(j)=(x(i)-x0(j))^2+(y(i)-y0(j))^2;
    end
mn=min(T);
a=find(T==mn);
num(a)=num(a)+1;
Kx(a,num(a))=x(i);
Ky(a,num(a))=y(i);
end
for m=1:ca
    x1(m)=mean(Kx(m,1:num(m)));
    y1(m)=mean(Ky(m,1:num(m)));
end 
z0=[x0==x1];
sm0=sum(z0);
z1=[y0==y1];
sm1=sum(z1);
t=0;
while sm0~=ca&&sm1~=ca
    for m=1:ca
    x0(m)=x1(m);
    y0(m)=y1(m);
    end
      num(:)=0;
     for i=1:length(x)
         a=0;
        for j=1:ca
         T(j)=(x(i)-x0(j))^2+(y(i)-y0(j))^2;
        end
      mn=min(T);
      a=find(T==mn);
      num(a)=num(a)+1;
      Kx(a,num(a))=x(i);
      Ky(a,num(a))=y(i);
     end
      for m=1:ca
    x1(m)=mean(Kx(m,1:num(m)));
    y1(m)=mean(Ky(m,1:num(m)));
    z0=[x0==x1];
    sm0=sum(z0);
    z1=[y0==y1];
    sm1=sum(z1);
      end
      t=t+1;
end
disp(['不经优化迭代了',num2str(t),'次']);
figure(ca)
hold on
for i=1:ca
plot(Kx(i,:),Ky(i,:),'*');
end
plot(x0,y0,'kp');
title([num2str(ca),'个聚类中心'])
xlabel('x坐标')
ylabel('y坐标')
axis([min(x) max(x) min(y) max(y)])
hold off
end
function D=leizuD(x,y)
for i=1:length(x)
    T(i,i)=eps;
    for j=i+1:length(y)
        T(i,j)=sqrt((x(i)-x(j))^2+(y(i)-y(j))^2);
        T(j,i)=T(i,j);
    end
end
D=max(max(T));
end
function R=leizuR(x0,y0,x,y)
      for i=1:length(x)
          D(i)=sqrt((x(i)-x0)^2+(y(i)-y0)^2);
      end
R=max(D);
end
