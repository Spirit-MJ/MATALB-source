close all
clear
clc
% x0=[3062 3663 4539 4843 5121 5842 6013 6437 6660 6774 7610];
x0=[8812.6 9763.1 9714.5 10302.6 10898.1 11952.2 12431.3 12486.9 12308.3 15479.4 17475.8];
% x0=[8812.6 9763.1 12308.3 15479.4 17475.8 9714.5 10302.6 10898.1 11952.2 12431.3 12486.9];
for i=2:length(x0)
    lamda(i-1)=x0(i-1)/x0(i);
end
e0=exp(-2/(length(x0)+1));
e1=exp(2/(length(x0)+1));
for i=1:length(x0)-1
    if lamda(i)<e0||lamda(i)>e1
        disp('潩潩潩潩????潩');
        break;
    end
end
T=0;
for i=1:length(x0)
   T=T+x0(i);
   x1(i)=T;
end
for i=2:length(x0)
    Y(i-1,1)=x0(i);
end
for i=1:length(x1)-1
    B(i,1)=-0.5*(x1(i)+x1(i+1));
    B(i,2)=1;
end
u=(inv((B')*B))*(B')*Y;
a=u(1);
b=u(2);
f=dsolve('Dy+a*y=b','y(0)=8812.6','x');
for i=1:length(x1)
    x=i-1;
    as=subs(f);
    x11(i)=eval(as);
end
disp(x11);
for i=2:length(x11)
    x01(i)=x11(i)-x11(i-1);
end
x01(1)=x11(1);
for i=1:length(x0)
    cx(i)=x0(i)-x01(i);
end
cxav=mean(cx);
x0av=mean(x0);
s01=0;s02=0;
for i=1:length(x)
   s01=s01+(x0(i)-x0av)^2;
   s02=s02+(cx(i)-cxav)^2;
end
s1=sqrt(s01/length(x0));    
s2=sqrt(s02/length(x01));   
C=s2/s1;
n=0;
for i=1:length(x0)
   if abs(cx(i)-cxav)<0.6745*s1
    n=n+1;
   end
end
p=n/length(x0);
if C<=0.35
    mc=1;
elseif C>0.35&&C<=0.5
    mc=2;
elseif C>0.5&&C<=0.65
    mc=3;
elseif C>0.65
    mc=4;
end
if p>=0.95
    mp=1;
elseif p>=0.8&&p<0.95
    mp=2;
elseif p>=0.7&&p<0.8
    mp=3;
elseif p<0.7
    mp=4;
end
disp(['C:',num2str(C)])
disp(['P:',num2str(p)])
m=max(mc,mp);
switch m
    case 4
        disp('潩�?�');
    case 3
        disp('潩?');
    case 2
        disp('�?�');
    case 1
        disp('潩');
end
x=11;
p0=subs(f);
x11(12)=eval(p0);
x01(12)=x11(12)-x11(11);
disp(['2020�?潩�?潩潩�?',num2str(x01(12)),'潩潩?�']);
hold on
grid
x=1:length(x0);
plot(x+2008,x0,'b-o');
x=1:length(x0)+1;
plot(x+2008,x01,'r-*');
xlabel('潩�')
ylabel('潩潩潩潩潩?�')
legend('?潩?','?潩?')
