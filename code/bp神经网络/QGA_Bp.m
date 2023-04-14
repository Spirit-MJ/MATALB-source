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
train_P=x0(1:9,:)';      %ѵ������������
train_T=y0(1:9,:)';      %ѵ�����������
%��������
test_P=x0(10:end,:)';    %���Լ���������
test_T=y0(10:end,:)';    %���Լ��������   
N=size(test_P,2);  

[train_p,ps_input]=mapminmax(train_P,0,1);%ѵ�����������ݹ�һ��
test_p=mapminmax('apply',test_P,ps_input);%���Լ��������ݹ�һ��
[train_t,ps_output]=mapminmax(train_T,0,1);%ѵ����������ݹ�һ��
np=10;%��Ⱥ��С
ng=100;%��������
variatenum=S;%��������
len=10;%���볤��
Youx=struct('fitx',0,'X',[],'Erjinzhi',[],'chrom',[]);%��¼�����Ӧ��ֵ��ʮ����ֵ,�����Ʊ��룬���ӱ��ر���
for i=1:2*np
    for j=1:variatenum*len
        x(i,j)=1/sqrt(2);%��ʼ����Ⱥ
    end
end
Erjinzhi=Celiang(x);%����Ⱥʵʩһ�β������õ������Ʊ���
Field=[repmat(len,1,S);repmat([-1;1],1,S);repmat([1;0;1;1],1,S)];%����������
%   FieldD = [len; lb; ub; code; scale; lbin; ubin]
%   len�ǰ�����Chrom�е�ÿ���Ӵ��ĳ��ȣ�ע��sum(len)=size(Chrom, 2)
%   lb��ub�ֱ���ÿ���������½���Ͻ硣
%   codeָ���Ӵ�����������ģ�1Ϊ��׼�Ķ����Ʊ��룬 0Ϊ���ױ���
%   scaleָ��ÿ���Ӵ���ʹ�õĿ̶ȣ�0��ʾ�����̶ȣ�1��ʾ�����̶ȡ�
%   lbin��ubinָ����ʾ��Χ���Ƿ�����߽硣0��ʾ��������1��ʾ�����߽�

X=bs2rv(Erjinzhi,Field);%�õ���Ⱥ�����ʮ����
fitx=objfun(X,train_p,train_t);%������Ӧ��ֵ
[best.fitx bestindex]=max(fitx);%�ҳ����ֵ
best.Erjinzhi=Erjinzhi(bestindex,:);%�ҳ������Ӧ�ȶ�Ӧ�Ķ�����
best.X=X(bestindex,:);%�ҳ������Ӧ�ȶ�Ӧ��ʮ����
best.x=x((2*bestindex-1):(2*bestindex),:);%�ҳ������Ӧ�ȶ�Ӧ�����ӱ�����
Youfitx(1)=best.fitx;
for k=2:ng
    Erjinzhi=Celiang(x);%����Ⱥʵʩһ�β���
    X=bs2rv(Erjinzhi,Field);%�õ���Ⱥ�����ʮ����
    fitx=objfun(X,train_p,train_t);%������Ӧ��
    x=XuanZhuan(x,fitx,best.Erjinzhi,Erjinzhi,best.fitx);%������ת��
    [newbestfitx,newbestindex]=max(fitx);%�ҳ��˴����
    if newbestfitx>best.fitx
       best.fitx=newbestfitx;
       best.Erjinzhi=Erjinzhi(newbestindex,:);
       best.x=x((2*bestindex-1):(2*bestindex),:);
       best.X=X(newbestindex,:);
    end
    Youfitx(k)=best.fitx;
    figure(1)
    line([k-1,k],[1/Youfitx(k-1),1/Youfitx(k)]);
    grid on
    xlabel('�Ŵ�����')
    ylabel('���')
    title('��������')
end
R2=0
%while R2<0.9
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
 w1=best.X(1:w1num);                           %����㵽������Ȩֵ
 b1=best.X(w1num+1:w1num+S2);                  %��������Ԫ��ֵ
 w2=best.X(w1num+S2+1:w1num+S2+w2num);         %�����㵽�����Ȩֵ
 b2=best.X(w1num+S2+w2num+1:w1num+S2+w2num+S3);%�������Ԫ��ֵ
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
%����ϵ��R^2
%R2=(N*sum(sim_T.*test_T)-sum(sim_T)*sum(test_T))^2/((N*sum((sim_T).^2)-(sum(sim_T))^2)*(N*sum((test_T).^2)-(sum(test_T))^2));
figure(2)
plot(1:N,test_T,'b:*',1:N,sim_T,'r-o')
legend('��ʵֵ','Ԥ��ֵ')
xlabel('Ԥ������')
ylabel('Ԥ��ֵ')
string={'���Լ�Ԥ�����Ա�';['R^2=',num2str(R2)]};
title(string)
%end
%����Ա�
result=[test_T',sim_T']
function Erjinzhi=Celiang(chrom)%��������������Ⱥʵʩ�������õ������Ʊ���
   [m,n]=size(chrom);
   for i=1:m/2
       for j=1:n
           Rand=rand(1);%����һ�������
           if Rand>(chrom(2*i-1,j)^2)%�����������������ڦ���ƽ��
               Erjinzhi(i,j)=1;
           else
               Erjinzhi(i,j)=0;
           end
       end
   end
end
function chrom=XuanZhuan(chrom,fitx,Youx,Erjinzhi,bestfitx)%������ת�ź���
   sizepop=size(chrom,1)/2;%��Ⱥ��С
   lenchrom=size(Erjinzhi,2);%���ӱ��ر��볤��
   for i=1:sizepop
       for j=1:lenchrom
           arfa=chrom(2*i-1,j);  %��
           beita=chrom(2*i,j);   %��
           x=Erjinzhi(i,j);      %�����������Ⱥ��ֵ��x
           b=Youx(j);
           if((x==0)&(b==1))||((x==1)&(b==1))
               delta=0;          %��ת�ǵĴ�С
               s=0;              %��ת�ǵķ��ţ���ת����
           elseif(x==0)&(b==1)&(fitx(i)<bestfitx)
               delta=0.01*pi;
               if arfa*beita>0
                   s=1;
               elseif arfa*beita<0
                   s=-1;
               elseif arfa==0
                   s=0;
               elseif beita==0
                   s=sign(randn);
               end
           elseif (x==0)&(b==1)&(fitx(i)>=bestfitx)
               delta=0.01*pi;
               if arfa*beita>0
                   s=-1;
               elseif arfa*beita<0
                   s=1;
               elseif arfa==0
                   s=sign(randn);
               elseif beita==0
                   s=0;
               end
           elseif (x==1)&(b==0)&(fitx(i)<bestfitx)
               delta=0.01*pi;
               if arfa*beita>0
                   s=-1;
               elseif arfa*beita<0
                   s=1;
               elseif arfa==0
                   s=sign(randn);
               elseif beita==0
                   s=0;
               end
           elseif (x==1)&(b==0)&(fitx(i)>=bestfitx)
               delta=0.01*pi;
               if arfa*beita>0
                   s=1;
               elseif arfa*beita<0
                   s=-1;
               elseif arfa==0
                   s=0;
               elseif beita==0
                   s=sign(randn);
               end
           else
               delta=0.01*pi;
               s=1;
           end
           fai=s*delta;      %��ת��
           U=[cos(fai) -sin(fai)
              sin(fai)  cos(fai)];%������ת��
          y=U*[arfa;beita];   %���º������λ
          chrom(2*i-1,j)=y(1);
          chrom(2*i,j)=y(2);
       end
   end        
end
function Obj=objfun(x,train_p,train_t)
    [M,N]=size(x);
    Obj=zeros(M,1);
    for i=1:M
        Obj(i)=Bp(x(i,:),train_p,train_t);
    end
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
    err=1/max(max(abs(Y-train_t)));
    %err=1./norm(Y-train_t);
end
