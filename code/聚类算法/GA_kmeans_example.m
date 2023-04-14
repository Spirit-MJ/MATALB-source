close all
clear
clc
%�����޸�ԭʼ��������Ҫ�޸ľ���Ĵ����ͱ߽纯��
wh = xlsread('E:\����\2.xlsx');
x=wh(:,1);
y=wh(:,2);
np=200;%��Ⱥ��С
pc=0.9;%�������
pm=0.1;%�������
ng=200;%��������
ggap=0.95;%��Ⱥ����
variatenum=9*2;%��������
ca=variatenum/2;%�����
len=10;%���볤��
yuanshi(x,y,variatenum);%���������Ż��ľ���ͼ
chrom=round(rand(np,variatenum*len));%��ʼ����Ⱥ
Field=[repmat(len,1,variatenum);repmat([min(x) min(y);max(x) max(y)],1,variatenum/2);repmat([1;0;0;0],1,variatenum)];%����������
%   FieldD = [len; lb; ub; code; scale; lbin; ubin]
%   len�ǰ�����Chrom�е�ÿ���Ӵ��ĳ��ȣ�ע��sum(len)=size(Chrom, 2)
%   lb��ub�ֱ���ÿ���������½���Ͻ硣
%   codeָ���Ӵ�����������ģ�1Ϊ��׼�Ķ����Ʊ��룬 0Ϊ���ױ���
%   scaleָ��ÿ���Ӵ���ʹ�õĿ̶ȣ�0��ʾ�����̶ȣ�1��ʾ�����̶ȡ�
%   lbin��ubinָ����ʾ��Χ���Ƿ�����߽硣0��ʾ��������1��ʾ�����߽�
for k=1:ng
    X=bs2rv(chrom,Field);
    for i=1:np
        fitx(i)=fitness(x,y,X(i,:),variatenum);%������Ӧ��ֵ
    end
    [Youfitx,index]=max(fitx);
    Youchrom=X(index,:);
    nx=XuanZe(chrom,ggap,np,fitx);%ѡ�����
    nx=JiaoCha(nx,pc);%�������
    nx=BianYi(nx,pm);%�������
    NX=bs2rv(nx,Field);
    for i=1:size(nx,1)
        newfitx(i)=fitness(x,y,NX(i,:),variatenum);%������Ӧ��ֵ
    end
    [NewYoufitx,Newindex]=max(newfitx);
    NewYouchrom=NX(Newindex,:);
    chrom=Reins(chrom,nx,fitx);%��ȫ��Ⱥ
    if Youfitx>=NewYoufitx
       trace(1, k)=1/Youfitx-1;
       trace(2, k)=1/mean(fitx)-1;
       Zuiyouchrom=Youchrom;
    else
       trace(1, k)=1/NewYoufitx-1;
       trace(2, k)=1/mean(newfitx)-1;
       Zuiyouchrom=NewYouchrom;
    end
    figure(2)
    if k>=2
       draw1=line([k-1,k],[trace(1, k-1),trace(1, k)]);
       grid on
       set(draw1,'color',[0 0 1]);
       xlabel('��������')
       ylabel('׼����')
       title('��������')
    end
end
x0=zeros(1,variatenum/2);
y0=zeros(1,variatenum/2);
for i=1:variatenum/2
    x0(i)=Zuiyouchrom(2*i-1);
    y0(i)=Zuiyouchrom(2*i);
end
num=zeros(1,ca);
for i=1:length(x)
    a=0;
    for j=1:ca
        T(j) = 1 - compute_IoU(x(i),y(i),x0(j),y0(j));
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
             T(j) = 1 - compute_IoU(x(i),y(i),x0(j),y0(j));
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
disp(['�Ż��������',num2str(t),'��']);
figure(3)
hold on
for i=1:ca
    plot(Kx(i,:),Ky(i,:),'*');
end
plot(x0,y0,'kp');
title('�Ŵ��㷨�Ż���ľ���')
xlabel('x����')
ylabel('y����')
axis([min(x) max(x) min(y) max(y)])
hold off
function fitx=fitness(x,y,X,variatenum)
  Kx=zeros(variatenum/2,length(x));
  Ky=zeros(variatenum/2,length(y));
  num=zeros(1,variatenum/2);
for i=1:variatenum/2
    x0(i)=X(2*i-1);
    y0(i)=X(2*i); 
end
for i=1:length(x)
    for j=1:variatenum/2
        T(j) = 1 - compute_IoU(x(i),y(i),x0(j),y0(j));
    end
mn=min(T);
a=find(T==mn);%��¼ĳһ��ĸ���
num(a)=num(a)+1;
Kx(a,num(a))=x(i);
Ky(a,num(a))=y(i);
end
J=0;
for i=1:variatenum/2
      for j=1:num(i)
          J = j + (1 - compute_IoU(Kx(i,j),Ky(i,j),x0(i),y0(i)));
      end
end
fitx=1/(1+J);
end
function nx=XuanZe(x,GGAP,np,fitx)%ѡ��
       GGAP1=floor(GGAP*np);  %ѡ�����µĸ�����
       [~,temp]=sort(fitx,'descend');%������Ӧ�Ƚ�������
       for j=1:GGAP1
           nx(j,:)=x(temp(j),:);%��ʤ��̭������Ӧ�ȴ�ĸ��屣������
       end
end
function nx=JiaoCha(nx,pc)%����
    [hang,lie]=size(nx);%����ѡ���������Ⱥ��С
    for i=1:2:hang-mod(hang,2)%�������ཻ������ֹѡ���������Ⱥ����Ϊ������
        if pc>=rand(1)
            Rand=randperm(lie,2);%���������������Ϊ�����Ƭ��
            Da=max(Rand);%���д����
            Xiao=min(Rand);%����С����
            temp=nx(i,Xiao:Da);
            nx(i,Xiao:Da)=nx(i+1,Xiao:Da);
            nx(i+1,Xiao:Da)=temp;
        end
    end   
end
function nx=BianYi(nx,pm)%����
    [hang,lie]=size(nx);%�����������Ⱥ�Ĵ�С
    for i=1:hang
        if pm>=rand
            Rand=randperm(lie,1);
            nx(i,Rand)=~nx(i,Rand);%�����������λ��
        end
    end
end
function nx=Reins(x,nx,fitx)%��ȫ��Ⱥ
    NIND=size(x,1);%��ʼ��Ⱥ��С
    NSel=size(nx,1);%ѡ���������Ⱥ�Ĵ�С
    [~,index]=sort(fitx,'descend');%���Ŵ�����ǰ�ľ�������
    nx=[x(index(1:NIND-NSel),:);nx];%���Ŵ�����ǰ���ŵ�һЩ���屣������
end
function yuanshi(x,y,variatenum)
ca=variatenum/2;
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
    T(j) = 1 - compute_IoU(x(i),y(i),x0(j),y0(j));
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
            T(j) = 1 - compute_IoU(x(i),y(i),x0(j),y0(j));
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
disp(['�����Ż�������',num2str(t),'��']);
figure(1)
hold on
for i=1:ca
plot(Kx(i,:),Ky(i,:),'*');
end
plot(x0,y0,'kp');
title('�����Ż��ľ���')
xlabel('x����')
ylabel('y����')
axis([min(x) max(x) min(y) max(y)])
hold off
end
function IoU = compute_IoU(W1,H1,W2,H2)
w_min = min(W1, W2);
h_min = min(H1, H2);

cross = w_min * h_min;

area1 = W1 * H1;
area2 = W2 * H2;
IoU = cross / (area1 + area2 - cross);
end