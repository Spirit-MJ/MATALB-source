close all
clear all
clc
C=[2 15 13 4
   10 4 14 15
   9 14 16 13
   7 8 11 9]
citys=1:size(C,1)+size(C,2);
D=[1000*ones(size(C,1),size(C,1)) C
   1000*ones(size(C,2),size(C,1)) 1000*ones(size(C,2),size(C,2))];
%��ʼ������
m=10;%��������
arfa=1;%��Ϣ����Ҫ�̶�����
beita=5;%����������Ҫ�̶�����
rou=0.1;%��Ϣ�ػӷ�����
Q=1;
eta=1./D;%��������
T=ones(length(citys),length(citys));%��Ϣ�ؾ���
Table=zeros(m,length(citys));%·����¼��
mingen=1;%����������ʼֵ
maxgen=100;%�����������ֵ
zuiyou_Time=zeros(1,maxgen);%��¼ÿ�ε��������ž���
zuiyou_route=zeros(maxgen,length(citys));%��¼ÿ�ε���������·��
while mingen<=maxgen
      start=zeros(m,1);
for i=1:m
      start(i)=randperm(length(citys),1);%mֻ���ϵ���ʼλ��
end
Table(:,1)=start;%mֻ���ϵ���ʼλ�÷ŵ�·����¼��ĵ�һ��
for i=1:m
    for j=2:length(citys)
        tabu=Table(i,1:(j-1));%�ѷ��ʵĳ��м���
        allow_index=~ismember(citys,tabu);
        allow=citys(allow_index);%û���߹��ĵص�
        p=allow;
        for k=1:length(allow)
            p(k)=T(tabu(end),allow(k))^arfa*eta(tabu(end),allow(k))^beita; %����ת�Ƹ���
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
    route=Table(i,:);
    for j=2:length(citys)
        chang(i)=chang(i)+D(route(j-1),route(j));%����mֻ���ϵ����߾���
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
   line([mingen-1,mingen],[zuiyou_chang(mingen-1)-3000,zuiyou_chang(mingen)-3000]);
end
xlabel('��������')
ylabel('����ֵ')
title('��������')
hold off
%������Ϣ��
derta=zeros(size(D,1),size(D,2));
 for i=1:m
     route=Table(i,:);
     for j=2:length(citys)
         derta(route(j-1),route(j))=derta(route(j-1),route(j))+Q/chang(i);%�ͷŵ���Ϣ��Ũ��      
     end
 end
 T=(1-rou)*T+derta;%������Ϣ��Ũ��
 mingen=mingen+1;
 Table=zeros(m,length(citys));
end
zuiduan_route=zuiyou_route(end,:)
Youmoney=0;
for i=2:length(citys)
    if D(zuiduan_route(i-1),zuiduan_route(i))~=1000
        Youmoney=Youmoney+D(zuiduan_route(i-1),zuiduan_route(i));
    end
end
You=zeros(size(C,1),size(C,2));
for i=1:size(C,1)
    You(zuiduan_route(2*i-1),zuiduan_route(2*i)-size(C,1))=1;
end
disp(['���Ž�Ϊ��',num2str(Youmoney)]);
disp('���ŵ�ָ�ɷ���Ϊ��')
You

