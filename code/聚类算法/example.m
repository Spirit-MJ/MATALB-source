clear all
close all
clc 
yuan=load('D:\MATLAB\bin\1.txt');
ca=input('请输入想要分成几类：\n');
y=yuan(:,1);
x=yuan(:,2);
for i=1:ca
    c=round(length(x)*rand(1));
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
        T(j)=sqrt((x(i)-x0(j))^2+(y(i)-y0(j))^2);
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
         T(j)=sqrt((x(i)-x0(j))^2+(y(i)-y0(j))^2);
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
disp(['迭代了',num2str(t),'次']);
format short
[h,l]=size(Kx);
temp=1;
for i=1:h
    for j=1:l
        if Kx(i,j)==0
            break;
        else
        d(temp)=152.5*sqrt((Kx(i,j)-x0(i))^2+(Ky(i,j)-y0(i))^2);
        temp=temp+1;
        end
    end
end
dd=mean(d)
hold on
plot(x0,y0,'kd');
for i=1:ca
plot(Kx(i,:),Ky(i,:),'+');
end
