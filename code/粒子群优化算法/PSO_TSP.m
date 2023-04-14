close all
clear
clc
city=[1304 2312
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
   2370 2975];%�����������
figure(1)
plot(city(:,1),city(:,2),'ms','LineWidth',2,'MarkerEdgeColor','k','MarkerFaceColor','g')
legend('����λ��')
title('���зֲ�ͼ','fontsize',12)
xlabel('km','fontsize',12)
ylabel('km','fontsize',12)
grid on
%������м����
n=size(city,1);%������Ŀ
D=zeros(n,n);%���о������
D=Juli(city);
nMax=200;                      %��������
indiNumber=5000;               %������Ŀ
individual=zeros(indiNumber,n);
%��ʼ������λ��
for i=1:indiNumber
    individual(i,:)=randperm(n);    
end
%������Ⱥ��Ӧ��
fitx=fitness(individual,city,D);
[~,index]=max(fitx);
tourPbest=individual;                              %��ǰ��������
tourGbest=individual(index,:);                    %��ǰȫ������
recordPbest=zeros(1,indiNumber);                %�������ż�¼
recordGbest=fitx(index);                        %Ⱥ�����ż�¼
xnew1=individual;
% ѭ��Ѱ������·��
L_best=zeros(1,nMax);
for N=1:nMax 
    % �������
    for i=1:indiNumber
        shu=randperm(n-1,2);
        chb1=min(shu);
        chb2=max(shu);
        cros=tourPbest(i,chb1:chb2);
        ncros=size(cros,2);      
        %ɾ���뽻��������ͬԪ��
        for j=1:ncros
            for k=1:n
                if xnew1(i,k)==cros(j)
                    xnew1(i,k)=0;
                    for t=1:n-k
                        temp=xnew1(i,k+t-1);
                        xnew1(i,k+t-1)=xnew1(i,k+t);
                        xnew1(i,k+t)=temp;
                    end
                end
            end
        end
        %���뽻������
        xnew1(i,n-ncros+1:n)=cros;
        %��·�����ȱ�������
        dist=0;
        for j=1:n-1
            dist=dist+D(xnew1(i,j),xnew1(i,j+1));
        end
        dist=dist+D(xnew1(i,1),xnew1(i,n));
        if (1/dist)>fitx(i)
            individual(i,:)=xnew1(i,:);
        end
        %��ȫ�����Ž��н���
        shu=randperm(n-1,2);
        chb1=min(shu);
        chb2=max(shu);
        cros=tourGbest(chb1:chb2); 
        ncros=size(cros,2);      
        %ɾ���뽻��������ͬԪ��
        for j=1:ncros
            for k=1:n
                if xnew1(i,k)==cros(j)
                    xnew1(i,k)=0;
                    for t=1:n-k
                        temp=xnew1(i,k+t-1);
                        xnew1(i,k+t-1)=xnew1(i,k+t);
                        xnew1(i,k+t)=temp;
                    end
                end
            end
        end
        %���뽻������
        xnew1(i,n-ncros+1:n)=cros;
        %��·�����ȱ�������
        dist=0;
        for j=1:n-1
            dist=dist+D(xnew1(i,j),xnew1(i,j+1));
        end
        dist=dist+D(xnew1(i,1),xnew1(i,n));
        if (1/dist)>fitx(i)
            individual(i,:)=xnew1(i,:);
        end
       %�������
        shu=randperm(n-1,2);
        c1=min(shu);
        c2=max(shu);
        temp=xnew1(i,c1);
        xnew1(i,c1)=xnew1(i,c2);
        xnew1(i,c2)=temp;
        %��·�����ȱ�������
        dist=0;
        for j=1:n-1
            dist=dist+D(xnew1(i,j),xnew1(i,j+1));
        end
        dist=dist+D(xnew1(i,1),xnew1(i,n));
        if (1/dist)>fitx(i)
            individual(i,:)=xnew1(i,:);
        end
        %������ת����
        r=randperm(n,2);%�������������
        mininverse=min(r);
        maxinverse=max(r);
        xnew1(i,mininverse:maxinverse)=xnew1(i,maxinverse:-1:mininverse);%��ת����,����Ի�λ��
        %��·�����ȱ�������
        dist=0;
        for j=1:n-1
            dist=dist+D(xnew1(i,j),xnew1(i,j+1));
        end
        dist=dist+D(xnew1(i,1),xnew1(i,n));
        if (1/dist)>fitx(i)
            individual(i,:)=xnew1(i,:);
        end
    end 
fitx=fitness(individual,city,D);%������Ӧ��ֵ
%���µ�ǰ���ź���ʷ����
for i=1:indiNumber
    if fitx(i)>recordPbest(i)
       recordPbest(i)=fitx(i);
       tourPbest(i,:)=individual(i,:);
    end
    if fitx(i)>recordGbest
       recordGbest=fitx(i);
       tourGbest=individual(i,:);
    end
end
L_best(N)=recordGbest;
if N>=2
   figure(2)
   line([N-1,N],[1/L_best(N-1),1/L_best(N)]);% �����ͼ
   title('�㷨ѵ������')
   xlabel('��������')
   ylabel('����')
   grid on
end
end

figure(3)
hold on
plot([city(tourGbest(1),1),city(tourGbest(n),1)],[city(tourGbest(1),2),...
    city(tourGbest(n),2)],'ms-','LineWidth',2,'MarkerEdgeColor','k','MarkerFaceColor','g')
hold on
for i=2:n
    plot([city(tourGbest(i-1),1),city(tourGbest(i),1)],[city(tourGbest(i-1),2),...
        city(tourGbest(i),2)],'ms-','LineWidth',2,'MarkerEdgeColor','k','MarkerFaceColor','g')
    hold on
end
legend('�滮·��')
scatter(city(:,1),city(:,2));
title('�滮·��','fontsize',10)
xlabel('km','fontsize',10)
ylabel('km','fontsize',10)
grid on

function D=Juli(city)
n=size(city,1);
for i=1:n
    for j=i:n
        if i==j
           D(i,j)=eps;
        else
            D(i,j)=((city(i,1)-city(j,1))^2+(city(i,2)-city(j,2))^2)^0.5;
        end
        D(j,i)=D(i,j);
    end
end
end
function fitx=fitness(x,city,D)
%x           input     ����
%cit         input     ��������
%D           input     ���о���
%fitx        output    ������Ӧ��ֵ 
m=size(x,1);
n=size(city,1);
fx=zeros(m,1);
for i=1:m
    for j=1:n-1
        fx(i)=fx(i)+D(x(i,j),x(i,j+1));
    end
    fx(i)=fx(i)+D(x(i,1),x(i,n));
end
fitx=1./fx;
end