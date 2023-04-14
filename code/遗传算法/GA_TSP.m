close all
clear
clc
x0=[1304 2312
   3639 1315
   4177 2244
   3712 1399
   3488 1535
   3326 1556
   3238 1229
   4196 1004
   4312 790
   4386 570
   3007 1970
   2562 1756
   2788 1491
   2381 1676
   1332 695
   3715 1678
   3918 2179
   4061 2370
   3780 2212
   3676 2578
   4029 2838
   4263 2931
   3429 1908
   3507 2367
   3394 2643
   3439 3201
   2935 3240
   3140 3550
   2545 2357
   2778 2826
   2370 2975];
np=2000;%��Ⱥ��С
len=length(x0);
pc=0.95;%�������
pm=0.05;%�������
ng=500;%��������
GGAP=0.9;%��Ⱥ����
D=distanse(x0);%�������
for i=1:np
    x(i,:)=randperm(len);%��ʼ����Ⱥ ʵ������
end

for k=1:ng
ObjV=PathLength(D,x);%���������Ⱥ�и���ľ���
[preObjV(k),index]=min(ObjV);%��¼ÿ�������ž���
youx=x(index,:);
%youx(1,:)=[6,7,13,12,14,15,1,29,31,30,27,28,26,25,20,21,22,18,3,17,19,24,11,23,16,4,8,9,10,2,5];
figure(1)%�Ż�����ͼ
if k>=2
   line([k-1,k],[preObjV(k-1),preObjV(k)]);
   xlabel('�Ŵ�����');
   ylabel('����');
   title('��������');
   grid on
end
fx=1./ObjV;  %������Ⱥ��ÿ���������Ӧ��
Drawpath(fx,x,x0);
nx=XuanZe(x,GGAP,np,fx); %ѡ�����
nx=JiaoCha(nx,pc,youx);%�������
nx=BianYi(nx,pm);%�������
nx=Reverse(nx,D);%������ת����
x=Reins(x,nx,ObjV);%�Ӵ�
end
a=find(fx==max(fx));%�ҳ�������Ӧ���������������
a3=a(1);%�ҳ�������Ӧ���������������
SS1=num2str(x(a3,1));
for i=2:len
    SS1=[SS1,'�D>',num2str(x(a3,i))];
end
disp(['����·��Ϊ��'])
disp(SS1);
disp(['�ܾ���Ϊ��'])
disp(ObjV(a3));
function D=distanse(a)%����������
    row=length(a);
    D=zeros(row);
    for i=1:row-1
        for j=i+1:row
            D(i,j)=((a(i,1)-a(j,1))^2+(a(i,2)-a(j,2))^2)^0.5;
            D(j,i)=D(i,j);
        end
        D(i,i)=eps;
    end
end
function len=PathLength(D,Chrom) %����ÿ����Ⱥ�и���ľ���
    row=size(D,2);%���볤��
    NIND=size(Chrom,1);%��Ⱥ��С
    len=zeros(NIND,1);
    for i=1:NIND
        p=[Chrom(i,:) Chrom(i,1)];
        i1=p(1:end-1);
        i2=p(2:end);
        for j=1:row
            len(i)=len(i)+D(i1(j),i2(j));%ÿ���������
        end
    end
end
function nx=XuanZe(x,GGAP,np,fx)%ѡ��
       GGAP1=floor(GGAP*np);  %ѡ�����µĸ�����      
       [~,temp]=sort(fx,'descend');%������Ӧ�Ƚ�������
       for j=1:GGAP1
           nx(j,:)=x(temp(j),:);%��ʤ��̭������Ӧ�ȴ�ĸ��屣������
       end
end
function SelCh=JiaoCha(SelCh,Pc,youx)%����ӳ�佻��
    NSel=size(SelCh,1);%����ѡ���������Ⱥ��С
    for i=1:NSel
        if Pc>=rand(1)
            [SelCh(i,:),youx]=intercross(SelCh(i,:),youx);
        end
    end
function [a,b]=intercross(a,b)
    L=length(a);%�������ĳ��ȣ����еĸ���
    r1=randsrc(1,2,1:L);%����2��1:L�������
    if r1(1)~=r1(2)
        a0=a;
        b0=b;
        s=min(r1);
        e=max(r1);
        for i=s:e
            a1=a;
            b1=b;
            a(i)=b0(i);%���������еĵ�i������
            b(i)=a0(i);%���������еĵ�i������
            x=find(a==a(i));%�ҳ����������λ�ú��뽻��������ظ������λ��
            y=find(b==b(i));%�ҳ����������λ�ú��뽻��������ظ������λ��
            i1=x(x~=i);%�ҳ�������ظ������λ��
            i2=y(y~=i);%�ҳ�������ظ������λ��
            if isempty(i1)==0
                a(i1)=a1(i);%���ò���ӳ��������ͻ
            end
            if isempty(i2)==0
                b(i2)=b1(i);%���ò���ӳ��������ͻ
            end
        end
    end
end    
end
function SelCh=BianYi(SelCh,Pm)%����
    [NSel,L]=size(SelCh);%�����������Ⱥ�Ĵ�С
    for i=1:NSel
        if Pm>=rand
            R0=randperm(L,2);%�����������λ��
            R=fliplr(R0);
            SelCh(i,R0)=SelCh(i,R);%�������������λ��
        end
    end
end
function SelCh=Reverse(SelCh,D)%������ת
    [row,col]=size(SelCh);%�����������Ⱥ�Ĵ�С
    ObjV=PathLength(D,SelCh);%���������������·�ߵľ���
    SelCh1=SelCh;
    for i=1:row
        r=randsrc(1,2,1:col);%�������������
        mininverse=min(r);
        maxinverse=max(r);
        SelCh1(i,mininverse:maxinverse)=SelCh1(i,maxinverse:-1:mininverse);%��ת����,����Ի�λ��
    end
    ObjV1=PathLength(D,SelCh1); %����·������
    index=ObjV1<ObjV;%�ҳ�������ת������·������С�ڽ�����ת����ǰ������
    SelCh(index,:)=SelCh1(index,:);%��������ת����������ԭ���ĸ����滻��
end
function Chrom=Reins(Chrom,SelCh,ObjV)%��ȫ��Ⱥ
    NIND=size(Chrom,1);%��ʼ��Ⱥ��С
    NSel=size(SelCh,1);%ѡ���������Ⱥ�Ĵ�С
    [~,index]=sort(ObjV);%���Ŵ�����ǰ�ľ�������
    Chrom=[Chrom(index(1:NIND-NSel),:);SelCh];%���Ŵ�����ǰ���ŵ�һЩ���屣������
end
function Drawpath(fx,x,x0)%��ͼ
a=find(fx==max(fx));%�ҳ�������Ӧ���������������
a3=a(1);%�ҳ�������Ӧ���������������
for i=1:length(x0)
    xx1(i)=x0(x(a3,i),1);
    yy1(i)=x0(x(a3,i),2);
end
xx=[xx1,x0(x(a3,1),1)];
yy=[yy1,x0(x(a3,1),2)];
figure(2)
plot(xx,yy,'r-*');
for i=1:length(x0)
    text(x0(i,1),x0(i,2),[' ' num2str(i)]);%��ͼ���ϰ�ÿ������ű���
end
text(x0(x(a3,1),1),x0(x(a3,1),2),'  ���','FontSize',13,'Color',[0 0 1]);%��ͼ���ϱ�����
text(x0(x(a3,end),1),x0(x(a3,end),2),'  �յ�','FontSize',13,'Color',[0 0 1]);%��ͼ���ϱ���յ�
box on
title('�켣ͼ');
xlabel('������');
ylabel('������');
axis('square','equal');
grid on
end
