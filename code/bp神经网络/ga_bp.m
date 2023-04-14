close all
clear all
clc
x0=[0.2286 0.1292 0.0720 0.1592 0.1335 0.0733 0.1159 0.0940 0.0522 0.1345 0.0090 0.1260 0.3619 0.0690 0.1828
    0.2090 0.0947 0.1393 0.1387 0.2558 0.0900 0.0771 0.0882 0.0393 0.1430 0.0126 0.1670 0.2450 0.0508 0.1328
    0.0442 0.0880 0.1147 0.0563 0.3347 0.1150 0.1453 0.0429 0.1818 0.0378 0.0092 0.2251 0.1516 0.0858 0.0670
    0.2603 0.1715 0.0702 0.2711 0.1491 0.1330 0.0968 0.1911 0.2545 0.0871 0.0060 0.1793 0.1002 0.0789 0.0909
    0.3690 0.2222 0.0562 0.5157 0.1872 0.1614 0.1425 0.1506 0.1310 0.0500 0.0078 0.0348 0.0451 0.0707 0.0880
    0.0359 0.1149 0.1230 0.5460 0.1977 0.1248 0.0624 0.0832 0.1640 0.1002 0.0059 0.1503 0.1837 0.1295 0.0700
    0.1759 0.2347 0.1829 0.1811 0.2922 0.0655 0.0774 0.0227 0.2056 0.0925 0.0078 0.1852 0.3501 0.1680 0.2668
    0.0724 0.1909 0.1340 0.2409 0.2842 0.0450 0.0824 0.1064 0.1909 0.1586 0.0116 0.1698 0.3644 0.2718 0.2494
    0.2634 0.2258 0.1165 0.1154 0.1074 0.0657 0.0610 0.2623 0.2588 0.1155 0.0050 0.0978 0.1511 0.2273 0.3220
    0.2101 0.0950 0.1298 0.1359 0.2601 0.1001 0.0753 0.0890 0.0389 0.1451 0.0128 0.1590 0.2452 0.0512 0.1319
    0.2593 0.1800 0.0711 0.2801 0.1501 0.1298 0.1001 0.1891 0.2531 0.0875 0.0058 0.1803 0.0992 0.0802 0.1002
    0.2599 0.2235 0.1201 0.0071 0.1102 0.0683 0.0621 0.2597 0.2602 0.1167 0.0048 0.1002 0.1521 0.2881 0.3205];
y0=[1 0 0
    1 0 0
    1 0 0
    0 1 0
    0 1 0
    0 1 0
    0 0 1
    0 0 1
    0 0 1
    1 0 0
    0 1 0
    0 0 1];
S1=size(x0,2);        %������Ԫ����
S2=2*size(x0,2)+1;    %������Ԫ����
S3=size(y0,2);        %�����Ԫ����

w1num=S1*S2;          %����㵽�������Ȩֵ����
w2num=S2*S3;          %�����㵽������Ȩֵ����
S=w1num+S2+w2num+S3;  %�Ż���������

%ѵ������
train_P=x0(1:510,:)';      %ѵ������������
train_T=y0(1:510,:)';      %ѵ�����������
%��������
test_P=x0(511:end,:)';    %���Լ���������
test_T=y0(511:end,:)';    %���Լ��������  
N=size(test_P,2);
      
[train_p,ps_input]=mapminmax(train_P,0,1);%ѵ�����������ݹ�һ��
test_p=mapminmax('apply',test_P,ps_input);%���Լ��������ݹ�һ��
[train_t,ps_output]=mapminmax(train_T,0,1);%ѵ����������ݹ�һ��

%�Ŵ��㷨�Ż�
np=40; %��Ⱥ��ģ
pc=0.7;%�������
pm=0.01;%�������
ng=10;%��������
ggap=0.95;%��Ⱥ����
len=10;%���볤��
FieldD=[repmat(len,1,S);repmat([-0.5;0.5],1,S);repmat([1;0;1;1],1,S)];%����������
%   FieldD = [len; lb; ub; code; scale; lbin; ubin]
%   len�ǰ�����Chrom�е�ÿ���Ӵ��ĳ��ȣ�ע��sum(len)=size(Chrom, 2)
%   lb��ub�ֱ���ÿ���������½���Ͻ硣
%   codeָ���Ӵ�����������ģ�1Ϊ��׼�Ķ����Ʊ��룬 0Ϊ���ױ���
%   scaleָ��ÿ���Ӵ���ʹ�õĿ̶ȣ�0��ʾ�����̶ȣ�1��ʾ�����̶ȡ�
%   lbin��ubinָ����ʾ��Χ���Ƿ�����߽硣0��ʾ��������1��ʾ�����߽�
x=round(rand(np,len*S));%��ʼ����Ⱥ
X=bs2rv(x,FieldD);%�����ʼ��Ⱥ��ʮ����ת��
objv=objfun(X,train_p,train_t);
for k=1:ng
    fx=1./objv;
    Objv=fx;
    preObjV(k)=max(Objv);%��¼ÿ��������
    a=find(fx==max(fx));%�ҳ�������Ӧ���������������
    Youx(k,:)=x(a(1),:);%�ҳ����ŵĸ���
    Zuiyou=Youx(k,:);
    if k>=2&&preObjV(k)<preObjV(k-1)
       preObjV(k)=preObjV(k-1);
       Zuiyou=Youx(k-1,:);%�ҳ����ŵĸ���
    end
    figure(1)%�Ż�����ͼ
    if k>=2
       line([k-1,k],[1/preObjV(k-1),1/preObjV(k)]);
       xlabel('�Ŵ�����');
       ylabel('���');
       title('��������');
    end
    nx=XuanZe(x,ggap,np,fx);%ѡ�����
    nx=JiaoCha(nx,pc);%�������
    nx=BianYi(nx,pm);%�������
    nx=Reins(x,nx,Objv);%��ȫ��Ⱥ
    x=nx;
    X=bs2rv(x,FieldD);
    objv=objfun(X,train_p,train_t);
end
Youx=Zuiyou;%�ҳ����ŵĸ���
x=bs2rv(Youx,FieldD);
%����������
net=feedforwardnet(S2);
net=configure(net,train_p,train_t);
net.layers{2}.transferFcn='logsig';
%�������������
net.trainparam.epochs=1000;%��������
net.trainparam.lr=0.1;%ѧϰ��
net.trainparam.goal=0.01;%ѵ�������

%��ֵ��������
 w1num=S1*S2;                             %����㵽�������Ȩֵ����
 w2num=S2*S3;                             %�����㵽������Ȩֵ����
 w1=x(1:w1num);                           %����㵽������Ȩֵ
 b1=x(w1num+1:w1num+S2);                  %��������Ԫ��ֵ
 w2=x(w1num+S2+1:w1num+S2+w2num);         %�����㵽�����Ȩֵ
 b2=x(w1num+S2+w2num+1:w1num+S2+w2num+S3);%�������Ԫ��ֵ
 net.iw{1,1}=reshape(w1,S2,S1);
 net.lw{2,1}=reshape(w2,S3,S2);
 net.b{1}=reshape(b1,S2,1);
 net.b{2}=reshape(b2,S3,1);

%�����µ�Ȩֵ����ֵ����ѵ��
net=train(net,train_p,train_t);
%�������
s_ga=sim(net,test_p);%�Ŵ��Ż���ķ�����
%����һ��
sim_T=mapminmax('reverse',s_ga,ps_output);
%��������
err=norm(sim_T-test_T);
%R^2����ϵ��
%����ϵ��R^2
R2=(N*sum(sim_T.*test_T)-sum(sim_T)*sum(test_T))^2/((N*sum((sim_T).^2)-(sum(sim_T))^2)*(N*sum((test_T).^2)-(sum(test_T))^2));
%����Ա�
result=[test_T',sim_T',err']
figure(2)
plot(1:N,test_T,'b:*',1:N,sim_T,'r-o')
legend('��ʵֵ','Ԥ��ֵ')
xlabel('Ԥ������')
ylabel('Ԥ��ֵ')
string={'���Լ�Ԥ�����Ա�';['R^2=',num2str(err)]};
title(string)

function nx=XuanZe(x,GGAP,np,fx)%ѡ��
       GGAP1=floor(GGAP*np);  %ѡ�����µĸ�����      
       [~,temp]=sort(fx,'descend');%������Ӧ�Ƚ�������
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
            nx(i,Rand)=rand(1)-0.5;%�����������λ��
        end
    end
end
function nx=Reins(x,nx,Objv)%��ȫ��Ⱥ
    NIND=size(x,1);%��ʼ��Ⱥ��С
    NSel=size(nx,1);%ѡ���������Ⱥ�Ĵ�С
    [~,index]=sort(Objv,'descend');%���Ŵ�����ǰ�ľ�������
    nx=[x(index(1:NIND-NSel),:);nx];%���Ŵ�����ǰ���ŵ�һЩ���屣������
end
function err=Bp(x,train_p,train_t)
    S1=size(train_p,1);        %������Ԫ����
    S2=2*size(train_p,1)+1;    %������Ԫ����
    S3=size(train_t,1);        %�����Ԫ����
   
    %�½�BP����
    %net=newff(train_p,train_t,S1)
    net=feedforwardnet(S2);
    net=configure(net,train_p,train_t);
    net.layers{2}.transferFcn='logsig';
    
    %�������������
    net.trainparam.epochs=1000;%��������
    net.trainparam.goal=0.01;%ѵ�������
    net.trainparam.lr=0.1;%ѧϰ��
    net.trainparam.show=10;%��ʵƵ�ʣ���������Ϊûѵ��10����ʾһ��
    net.trainparam.showWindow=0;%��n�����ں���ʾһ���������ߵı仯
    
    %���ó�ʼȨֵ����ֵ
    w1num=S1*S2;                             %����㵽�������Ȩֵ����
    w2num=S2*S3;                             %�����㵽������Ȩֵ����
    w1=x(1:w1num);                           %����㵽������Ȩֵ
    b1=x(w1num+1:w1num+S2);                  %��������Ԫ��ֵ
    w2=x(w1num+S2+1:w1num+S2+w2num);         %�����㵽�����Ȩֵ
    b2=x(w1num+S2+w2num+1:w1num+S2+w2num+S3);%�������Ԫ��ֵ
    net.iw{1,1}=reshape(w1,S2,S1);
    net.lw{2,1}=reshape(w2,S3,S2);
    net.b{1}=reshape(b1,S2,1);
    net.b{2}=reshape(b2,S3,1);
    
    %ѵ��������
    net=train(net,train_p,train_t);

    %��������
    Y=sim(net,train_p);
    err=norm(Y-train_t);
end
function Obj=objfun(x,train_p,train_t)
    [M,N]=size(x);
    Obj=zeros(M,1);
    for i=1:M
        Obj(i)=Bp(x(i,:),train_p,train_t);
    end
end
