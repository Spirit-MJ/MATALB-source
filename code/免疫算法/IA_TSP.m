close all
clear
clc
city=[1304,2312
      3639,1315
      4177,2244
      3712,1399
      3488,1535
      3326,1556
      3238,1229
      4196,1044
      4312,790
      4386,570
      3007,1970
      2562,1756
      2788,1491
      2381,1676
      1332,695
      3715,1678
      3918,2179
      4061,2370
      3780,2212
      3676,2578
      4029,2838
      4263,2931
      3429,1908
      3507,2376
      3394,2643
      3439,3201
      2935,3240
      3140,3550
      2545,2357
      2778,2826
      2370,2975];
np=250;           % ��Ⱥ��ģ
overbest=50;          % ���������
ng=500;            % ��������
pc=0.95;           % �������
pm=0.2;        % �������
ps=0.9;              % ���������۲���
T=0.7;         % ���ƶȴ��ڷ�ֵ��������ֵΪ0.7
len=length(city);      %���볤��
M=np+overbest;
D=distanse(city);
%ʶ��ԭ,����Ⱥ��Ϣ����Ϊһ���ṹ��
individuals=struct('fitx',zeros(1,M), 'C',zeros(1,M),'P',zeros(1,M),'x',[]);
%������ʼ����Ⱥ
for i=1:M
    individuals.x(i,:)=randperm(len);%��ʼ����Ⱥ
end
trace=[]; %��¼ÿ�����������Ӧ�Ⱥ�ƽ����Ӧ��
for iii=1:ng
     %����Ⱥ����������
     for i=1:M
         individuals.fitx(i)=fitness(individuals.x(i,:),D);%���㿹���뿹ԭ�׺���(��Ӧ��ֵ��
     end
     for i=1:M
         individuals.C(i)=concentration(i,M,individuals,T);%���㿹��Ũ��
     end
     %�ۺ��׺ͶȺ�Ũ�����ۿ�������̶ȣ��ó���ֳ����
     individuals.P=excellence(individuals,M,ps);    
     %��¼������Ѹ������Ⱥƽ����Ӧ��
     [best,index]=max(individuals.fitx);   % �ҳ�������Ӧ�� 
     bestchrom=individuals.x(index,:);    % �ҳ����Ÿ���
     average=mean(individuals.fitx);       % ����ƽ����Ӧ��
     trace=[trace;1/best,1/average];              %��¼
     Drawpath(bestchrom,city);
     %����P���γɸ���Ⱥ�����¼���⣨���뾫Ӣ�������ԣ�����s���ƣ�
     bestindividuals=member(individuals,M,overbest);   % ���¼����
     individuals=member(individuals,M,np);      % �γɸ���Ⱥ
     % ѡ�񣬽��棬����������ټ��������п��壬��������Ⱥ
     individuals=XuanZe(individuals,np);% ѡ��
     individuals.x=JiaoCha(individuals.x,pc,bestchrom);% ����
     individuals.x=BianYi(individuals.x,pm);% ����
     individuals.x=Reverse(individuals.x,D);%������ת����
     individuals=Reins(individuals,np,bestindividuals,overbest);% ���������п���      
     if iii>=2
      % ���������㷨��������
        figure(1)
        hold on
        draw1=line([iii-1,iii],[trace(iii-1,1),trace(iii,1)]);
        draw2=line([iii-1,iii],[trace(iii-1,2),trace(iii,2)]);
        grid on
        set(draw1,'color',[0 0 1]);
        set(draw2,'color',[1 0 0]);
        xlabel('��������','fontsize',12)
        ylabel('����','fontsize',12)
        title('�����㷨��������','fontsize',12)
        legend('������Ӧ��ֵ','ƽ����Ӧ��ֵ')
     end
end
SS1=num2str(bestchrom(1));
for i=2:len
    SS1=[SS1,'�D>',num2str(bestchrom(i))];
end
disp(['����·��Ϊ��'])
disp(SS1);
disp(['�ܾ���Ϊ��'])
disp(1/best);
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
function len=PathLength(D,x) %����ÿ����Ⱥ�и���ľ���
    [NIND,row]=size(x);%��Ⱥ��С
    len=zeros(NIND,1);
    for j=1:NIND
        d=0;
        for i=1:row-1
            d=d+D(x(j,i),x(j,i+1));
        end
        d=d+D(x(j,end),x(j,1));
        len(j)=d;
    end
end
function fitx=fitness(x,D)
d=0;
for i=1:length(x)-1
     d=d+D(x(i),x(i+1));
end
d=d+D(x(end),x(1));
fitx=1/d;
end
function Cv=concentration(i,M,individuals,T)%���㿹��Ũ��
% �������Ũ��ֵ
% i              input      ��i������
% M              input      ��Ⱥ��ģ
% individuals    input     ����
% Concentration  output     Ũ��ֵ
Cv=0;
for j=1:M
    xsd=similar(individuals.x(i,:),individuals.x(j,:));%��i��������Ⱥ���������ƶ� 
    if xsd>T
        Cv=Cv+1;
    end
end
Cv=Cv/M;
end
function Svs=similar(v,s)%���㿹���뿹��֮����׺���
k=zeros(1,length(v));
for i=1:length(v)
    if find(v(i)==s(i))
        k(i)=1;
    end
end
Svs=sum(k)/length(v);
end
function P=excellence(individuals,M,ps)%����������ֳ����
% ������己ֳ����
% individuals    input      ��Ⱥ
% M              input      ��Ⱥ��ģ
% ps             input      ���������۲���
% exc            output     ��ֳ����
A=individuals.fitx;
sumA=sum(A);
C=individuals.C;
sumC=sum(C);
for i=1:M
    P(i)=ps*A(i)/sumA+(1-ps)*C(i)/sumC; 
end
end
function rets=member(individuals,m,n)%��¼������ֳ���ʴ�ĺ���Ӧ�Ⱥ�����Ŀ���
% ��ʼ�������,����excellence����Ⱥ���и���Ӧ�ȵ����ƶȵ�overbest�������������
% m                  input          ������
% n                  input          ����������\����Ⱥ��ģ
% individuals        input          ����Ⱥ
% bestindividuals    output         �����\����Ⱥ
% ��Ӣ�������ԣ���fitness��õ�s�������ȴ���������������Ũ�ȸ߶�����̭
s=15;
rets=struct('fitx',zeros(1,n), 'C',zeros(1,n),'P',zeros(1,n),'x',[]);
[~,index]=sort(individuals.fitx,'descend');%������Ӧ�ȴӴ�С����
for i=1:s
    rets.fitx(i)=individuals.fitx(index(i));   
    rets.C(i)=individuals.C(index(i));
    rets.P(i)=individuals.P(index(i));
    rets.x(i,:)=individuals.x(index(i),:);
end
% ʣ��m-s������
leftindividuals=struct('fitness',zeros(1,m-s), 'concentration',zeros(1,m-s),'excellence',zeros(1,m-s),'chrom',[]);
for k=1:m-s
    leftindividuals.fitness(k)=individuals.fitx(index(k+s));   
    leftindividuals.concentration(k)=individuals.C(index(k+s));
    leftindividuals.excellence(k)=individuals.P(index(k+s));
    leftindividuals.chrom(k,:)=individuals.x(index(k+s),:);
end
% ��ʣ�࿹�尴excellenceֵ����,��������ֳ���ʴ�ĸ��屣������
[~,index]=sort(leftindividuals.excellence,'descend');%����������ֳ���ʴӴ�С����
% ��ʣ�࿹��Ⱥ�а�excellence��ѡn-s����õĸ���
for i=s+1:n
    rets.fitx(i)=leftindividuals.fitness(index(i-s));
    rets.C(i)=leftindividuals.concentration(index(i-s));
    rets.P(i)=leftindividuals.excellence(index(i-s));
    rets.x(i,:)=leftindividuals.chrom(index(i-s),:);
end
end
function ret=XuanZe(individuals,np)%ѡ��
% ���̶�ѡ��
Excellence=individuals.P;
pselect=Excellence./sum(Excellence);
for j=1:np
      sita=rand(1);
      while sita==0    
            sita=rand(1);        
      end
      for i=1:np
          if sita<=pselect(i)
             individuals.x(j,:)=individuals.x(i,:);
             individuals.fitx(j)=individuals.fitx(i);
             individuals.C(j)=individuals.C(i);
             individuals.P(j)=individuals.P(i);
             break;
           end
      end
end
ret=individuals;
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
function newindividuals=Reins(individuals,np,bestindividuals,overbest)%��ȫ��Ⱥ
m=np+overbest;
newindividuals=struct('fitx',zeros(1,m), 'C',zeros(1,m),'P',zeros(1,m),'x',[]);
% �Ŵ������õ��Ŀ���
for i=1:np
    newindividuals.fitx(i)=individuals.fitx(i);   
    newindividuals.C(i)=individuals.C(i);   
    newindividuals.P(i)=individuals.P(i);   
    newindividuals.x(i,:)=individuals.x(i,:);   
end
% ������п���
for i=np+1:m
    newindividuals.fitx(i)=bestindividuals.fitx(i-np);   
    newindividuals.C(i)=bestindividuals.C(i-np);   
    newindividuals.P(i)=bestindividuals.P(i-np);   
    newindividuals.x(i,:)=bestindividuals.x(i-np,:);   
end
end
function Drawpath(bestchrom,city)%��ͼ
for i=1:length(city)
    xx1(i)=city(bestchrom(i),1);
    yy1(i)=city(bestchrom(i),2);
end
xx=[xx1,city(bestchrom(1),1)];
yy=[yy1,city(bestchrom(1),2)];
figure(2)
plot(xx,yy,'r-*');
for i=1:length(city)
    text(city(i,1),city(i,2),[' ' num2str(i)]);%��ͼ���ϰ�ÿ������ű���
end
text(city(bestchrom(1),1),city(bestchrom(1),2),'  ���','FontSize',13,'Color',[0 0 1]);%��ͼ���ϱ�����
text(city(bestchrom(end),1),city(bestchrom(end),2),'  �յ�','FontSize',13,'Color',[0 0 1]);%��ͼ���ϱ���յ�
box on
title('�켣ͼ');
xlabel('������');
ylabel('������');
axis('square','equal');
grid on
end  
