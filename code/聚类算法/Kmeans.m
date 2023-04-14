close all
clear 
clc
ca=input('请输入想要分成几类：\n');
x=rand(1,1000);
y=rand(1,1000);
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
figure(1)
hold on
for i=1:ca
plot(Kx(i,:),Ky(i,:),'*');
end
plot(x0,y0,'kp');
title('不经优化的聚类')
xlabel('x坐标')
ylabel('y坐标')
axis([min(x) max(x) min(y) max(y)])
hold off