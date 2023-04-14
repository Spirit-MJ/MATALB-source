close all
clear
clc
x=[1304 2312;
   3639 1315;
   4177 2244;
   3712 1399;
   3488 1535;
   3326 1556;
   3238 1229;
   4196 1004;
   4312 790;
   4386 570;
   3007 1970;
   2562 1756;
   2788 1491;
   2381 1676;
   1332 695;
   3715 1678;
   3918 2179;
   4061 2370;
   3780 2212;
   3676 2578;
   4029 2838;
   4263 2931;
   3429 1908;
   3507 2367;
   3394 2643;
   3439 3201;
   2935 3240;
   3140 3550;
   2545 2357;
   2778 2826;
   2370 2975];

for i=1:length(x)
    for j=1:length(x)
        if i==j
        D(i,j)=eps;
        else 
        D(i,j)=sqrt((x(i,1)-x(j,1))^2+(x(i,2)-x(j,2))^2);%����
        end
    end
end
%��ʼ������
m=35;%��������
arfa=1;%��Ϣ����Ҫ�̶�����
beita=2;%����������Ҫ�̶�����
rou=0.2;%��Ϣ�ػӷ�����
Q=1;
eta=1./D;%��������
T=ones(length(x));%��Ϣ�ؾ���
Table=zeros(m,length(x));%·����¼��
mingen=1;%����������ʼֵ
maxgen=300;%�����������ֵ
zuiyou_chang=zeros(1,maxgen);%��¼ÿ�ε��������ž���
zuiyou_route=zeros(maxgen,length(x));%��¼ÿ�ε���������·��
while mingen<=maxgen
    start=zeros(m,1);
for i=1:m
    start(i)=randperm(length(x),1);%mֻ���ϵ���ʼλ��
end
    Table(:,1)=start;%mֻ���ϵ���ʼλ�÷ŵ�·����¼��ĵ�һ��
    citys=1:length(x);
for i=1:m
    for j=2:length(x)
        tabu=Table(i,1:(j-1));%�ѷ��ʵĳ��м���
        allow_index=~ismember(citys,tabu);
        allow=citys(allow_index);%û���߹��ĵص�
        p=allow;
        for k=1:length(allow)
            p(k)=T(tabu(end),allow(k))^arfa * eta(tabu(end),allow(k))^beita; %����ת�Ƹ���
        end
        p=p/sum(p);
        pc=cumsum(p);%���̶�
        mubiao_index=find(pc>=rand(1));
        mubiao=allow(mubiao_index(1));%�ҳ�������һ�������ȥ�ĵ�
        Table(i,j)=mubiao;%��¼ÿֻ������·��
    end
end
chang=zeros(m,1);
for i=1:m
    route=[Table(i,:) Table(i,1)];
    for j=1:length(x)
        chang(i)=chang(i)+D(route(j),route(j+1));%����mֻ���ϵ����߾���
    end
end
%����·��
 if mingen==1
     [min_chang,min_index]=min(chang);%�ҳ���һ��mֻ�����ߵ����·��
     zuiyou_chang(mingen)=min_chang;%��¼ÿһ����̵ľ���
     zuiyou_route(mingen,:)=Table(min_index,:);%��¼ÿһ�����ŵ�·��
 else 
     [min_chang, min_index]=min(chang);
     if min_chang<zuiyou_chang(mingen-1)%�жϴ˴ξ����Ƿ���ϴξ����
         zuiyou_chang(mingen)=min_chang;%�ǵĻ���¼�˴���̾���
         zuiyou_route(mingen,:)=Table(min_index,:);%�ǵĻ���¼�˴�����·��
     else
         zuiyou_chang(mingen)=zuiyou_chang(mingen-1);%����Ļ���¼��һ����̾���
         zuiyou_route(mingen,:)=zuiyou_route(mingen-1,:);%����Ļ���¼��һ������·��
     end
 end
    figure(1)%��������ͼ��
    hold on
    if mingen>=2
        line([mingen-1,mingen],[zuiyou_chang(mingen-1),zuiyou_chang(mingen)]);
    end
    xlabel('��������')
    ylabel('����ֵ')
    title('��������')
    hold off
%������Ϣ��
    derta=zeros(length(x));
     for i=1:m
         route=[Table(i,:) Table(i,1)];
         for j=1:length(x)
             derta(route(j),route(j+1))=derta(route(j),route(j+1))+Q/chang(i);%�ͷŵ���Ϣ��Ũ��      
         end
     end
     T=(1-rou)*T+derta;%������Ϣ��Ũ��
     mingen=mingen+1;
     Table=zeros(m,length(x));
end
zuiduan_route=zuiyou_route(end,:);
zuiduan_chang=zuiyou_chang(end);
disp(['��̵�·��Ϊ��',num2str(zuiduan_route)]);
disp(['��̵ľ���Ϊ��',num2str(zuiduan_chang)]);
for i=1:length(x)
    cityx(i)=x(zuiduan_route(i),1);
    cityy(i)=x(zuiduan_route(i),2);
end
cityx=[cityx x(zuiduan_route(1),1)];
cityy=[cityy x(zuiduan_route(1),2)];
hold off
figure(2)
plot(cityx,cityy,'r-*');
xlabel('����λ�ú�����')
ylabel('����λ��������')
title('���·��')

