% �����Ż��㷨��������������ѡַ�е�Ӧ��
close all
clear all
clc
city_coordinate=[1304,2312;
                 3639,1315;
                 4177,2244;
                 3712,1399;
                 3488,1535;
                 3326,1556;
                 3238,1229;
                 4196,1044;
                 4312,790;
                 4386,570;
                 3007,1970;
                 2562,1756;
                 2788,1491;
                 2381,1676;
                 1332,695;
                 3715,1678;
                 3918,2179;
                 4061,2370;
                 3780,2212;
                 3676,2578;
                 4029,2838;
                 4263,2931;
                 3429,1908;
                 3507,2376;
                 3394,2643;
                 3439,3201;
                 2935,3240;
                 3140,3550;
                 2545,2357;
                 2778,2826;
                 2370,2975];
carge=[20,90,90,60,70,70,40,90,90,70,60,40,40,40,20,80,90,70,100,50,50,50,80,70,80,40,40,60,70,50,30];           
np=100;           % ��Ⱥ��ģ
overbest=10;          % ���������
ng=200;            % ��������
pc=0.95;           % �������
pm=0.2;        % �������
ps=0.9;              % ���������۲���
T=0.7;         % ���ƶȴ��ڷ�ֵ��������ֵΪ0.7
CenterNum=6;             % ����������
M=np+overbest;      
%ʶ��ԭ,����Ⱥ��Ϣ����Ϊһ���ṹ��
individuals=struct('fitx',zeros(1,M), 'C',zeros(1,M),'P',zeros(1,M),'x',[]);
%������ʼ����Ⱥ
individuals.x=popinit(M,CenterNum);%��ʼ����Ⱥ
trace=[]; %��¼ÿ�����������Ӧ�Ⱥ�ƽ����Ӧ��
for iii=1:ng
     %����Ⱥ����������
     for i=1:M
         individuals.fitx(i)=fitness(individuals.x(i,:),city_coordinate,carge);%���㿹���뿹ԭ�׺���(��Ӧ��ֵ��
         individuals.C(i)=concentration(i,M,individuals,T);%���㿹��Ũ��
     end
     %�ۺ��׺ͶȺ�Ũ�����ۿ�������̶ȣ��ó���ֳ����
     individuals.P=excellence(individuals,M,ps);    
     %��¼������Ѹ������Ⱥƽ����Ӧ��
     [best,index]=max(individuals.fitx);   % �ҳ�������Ӧ�� 
     bestchrom=individuals.x(index,:);    % �ҳ����Ÿ���
     average=mean(individuals.fitx);       % ����ƽ����Ӧ��
     trace=[trace;1/best,1/average];              %��¼
     %����excellence���γɸ���Ⱥ�����¼���⣨���뾫Ӣ�������ԣ�����s���ƣ�
     bestindividuals=member(individuals,M,overbest);   % ���¼����
     individuals=member(individuals,M,np);      % �γɸ���Ⱥ
     % ѡ�񣬽��棬����������ټ��������п��壬��������Ⱥ
     individuals=XuanZe(individuals,np);% ѡ��
     individuals.x=JiaoCha(pc,individuals.x,np,CenterNum);% ����
     individuals.x=BianYi(pm,individuals.x,np,CenterNum);% ����
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
        ylabel('��Ӧ��ֵ','fontsize',12)
        title('�����㷨��������','fontsize',12)
        legend('������Ӧ��ֵ','ƽ����Ӧ��ֵ')
     end
end
% ������������ѡַͼ
for i=1:31
    distance(i,:)=dist(city_coordinate(i,:),city_coordinate(bestchrom,:)');
end
[a,b]=min(distance');
index=cell(1,CenterNum);
for i=1:CenterNum
    index{i}=find(b==i);%����������͵�ĵ�ַ
end
cargox=city_coordinate(bestchrom,1);
cargoy=city_coordinate(bestchrom,2);
figure(2)
for j=1:length(index)
    for i=1:length(index{j})
        A=[city_coordinate(index{j}(i),1),city_coordinate(index{j}(i),2)];
        B=[cargox(j),cargoy(j)];
        c=[A;B];
        plot(c(:,1),c(:,2),'b-*')
        hold on
    end
end
title('�����Ż��㷨ѡַ');
xlabel('x����');
ylabel('y����');
function pop=popinit(n,len)
for i=1:n
    flag=0;
    while flag==0
        pop(i,:)=randperm(31,len);
        flag=test(pop(i,:));
    end
end
end
function flag=test(code)
city_coordinate=[1304,2312;3639,1315;4177,2244;3712,1399;3488,1535;3326,1556;3238,1229;4196,1044;4312,790;4386,570;
                 3007,1970;2562,1756;2788,1491;2381,1676;1332,695;3715,1678;3918,2179;4061,2370;3780,2212;3676,2578;
                 4029,2838;4263,2931;3429,1908;3507,2376;3394,2643;3439,3201;2935,3240;3140,3550;2545,2357;2778,2826;2370,2975];
flag=1;
if max(max(dist(city_coordinate(code,:)')))>3000%�������ľ���Լ��
    flag=0;
end
end
function fitx=fitness(individual,city_coordinate,carge)
%�ҳ�������͵�
%dist(A,B)����A��ÿ����������B��ÿ��������֮��ŷ�Ͼ��룬A��������ά���������B��������ά��
for i=1:31
    distance(i,:)=dist(city_coordinate(i,:),city_coordinate(individual,:)');
end
a=min(distance');%�ҳ�ÿ���㵽�������ľ�����̵ĵ��������
for i=1:31  
    expense(i)=carge(i)*a(i);%�������
end
fx=sum(expense)+4.0e+4*length(find(a>3000));%�������3000ȡһ���ͷ�ֵ
fitx=1/fx;
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
    if find(v(i)==s)
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
s=3;
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
function ret=JiaoCha(pcross,chrom,sizepop,len)%����
for i=1:sizepop   
    pick=rand(1);
    while pick==0
        pick=rand(1);
    end
    if pick>pcross
        continue;
    end
    shu=randperm(sizepop,2);% �ҳ��������
    index(1)=shu(1);
    index(2)=shu(2);
    % ѡ�񽻲�λ��
    pos=ceil(len*rand);
    while pos==1||pos==len
        pos=ceil(len*rand);
    end
    %�ҳ��������������
    chrom1=chrom(index(1),:);
    chrom2=chrom(index(2),:);
    k=chrom1(pos:len);
    chrom1(pos:len)=chrom2(pos:len);
    chrom2(pos:len)=k; 
    
    % ����Լ��������������Ⱥ
    flag1=test(chrom(index(1),:));
    flag2=test(chrom(index(2),:));
    if flag1*flag2==1
        chrom(index(1),:)=chrom1;
        chrom(index(2),:)=chrom2;
    end   
end
ret=chrom;
end
function ret=BianYi(pm,chrom,np,len)%����
for i=1:np   
    pick=rand(1);% �������
    while pick==0
        pick=rand(1);
    end
    if pick>pm   % �ж��Ƿ����
        continue;
    end
    index=randperm(np,1);
    pos=randperm(len,1);
    while pos==1
          pos=randperm(len,1);
    end
    nchrom=chrom(index,:);
    nchrom(pos)=randperm(31,1);
    while length(unique(nchrom))==(len-1)
          nchrom(pos)=randperm(31,1);
    end
    flag=test(nchrom);
    if flag==1
        chrom(index,:)=nchrom;
    end
    chrom(index,:)=nchrom;
end
ret=chrom;
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
   