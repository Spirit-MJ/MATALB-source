function xuanzhi
x0=[0,0];
format short g
[x,fval]=fmincon(@fun,x0,[],[],[],[],[],[],[])
end
function f=fun(x)
y=[3,8;
    8,2;
    2,5;
    6,4;
    8,8];
w=[2000,3000,2500,1000,1500];
a=[0.4,0.4,0.6,0.6,0.6];
T=0;
for i=1:5
    T=T+w(i)*a(i)*sqrt((x(1)-y(i,1))^2+(x(2)-y(i,2))^2);
end
f=T;
end
