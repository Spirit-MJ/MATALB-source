close all
clear
clc
% np=input('��������Ⱥ��������\n');
% l=input('��������볤�ȣ�\n');
% pc=input('�����뽻����ʣ�\n');
% pm=input('�����������ʣ�\n');
% ng=input('��������Ⱥ����������\n');
np=20;
l=10;
pc=0.9;
pm=0.05;
ng=200;
x=zeros(np,l);
%������ʼ��Ⱥ
for i=1:np
    x00=randperm(7,1);
    x0=dec2bin(x00,3);
    x0=[x0;repmat(' ',1,3)];
    y=str2num(x0(:).');
    x(i,1:3)=y;
     x00=randperm(7,1);
    x0=dec2bin(x00,3);
    x0=[x0;repmat(' ',1,3)];
    y=str2num(x0(:).');
     x(i,4:6)=y;
end
for i=1:np
    sum0=0;
    sum00=0;
    for j=1:3
      sum0=sum0+2^(3-j)*x(i,j);
    end
    for j=4:6
         sum00=sum00+2^(6-j)*x(i,j);
    end
    fx(i)=sum0^2+sum00^2;%������Ӧֵ
end
for n=1:ng
 hold on
 plot(n,1/max(fx),'b.');
for i=1:np
p(i)=fx(i)/sum(fx);%����ѡ�����
end
for i=1:np
    ppx0=0;
    for j=1:i
     ppx0=ppx0+p(j);
    end
    ppx(i)=ppx0;%�ۼƸ���
end
for m=1:np
sita=rand();
for i=1:np
    if sita<=ppx(i)
        father=i;
        break;
    end
end
mother=randperm(np,1);
cut=randperm(l-1,1);
r1=rand();
if r1<=pc%����
    nx(m,1:cut)=x(father,1:cut);
    nx(m,cut+1:l)=x(mother,cut+1:l);
    r2=rand();
if r2<=pm%����
    mut=randperm(l,1);
    nx(m,mut)=~nx(m,mut);
end
else
    nx(m,:)=x(father,:);
end
end
x=nx;
for i=1:np
    sum0=0;
    sum00=0;
    for j=1:3
      sum0=sum0+2^(3-j)*x(i,j);
    end
    for j=4:6
      sum00=sum00+2^(6-j)*x(i,j);
    end
    fx(i)=sum0^2+sum00^2;%�Ӵ���Ӧֵ
end
end
x
a=find(fx==max(fx));
disp(x(a,:));
disp(fx(a(1)));

