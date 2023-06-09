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
S1=size(x0,2);        %输入神经元个数
S2=2*size(x0,2)+1;    %隐层神经元个数
S3=size(y0,2);        %输出神经元个数

w1num=S1*S2;          %输入层到隐含层的权值个数
w2num=S2*S3;          %隐含层到输出层的权值个数
S=w1num+S2+w2num+S3;  %优化变量个数

%训练数据
train_P=x0(1:9,:)';      %训练集输入数据
train_T=y0(1:9,:)';      %训练集输出数据
%测试数据
test_P=x0(10:end,:)';    %测试集输入数据
test_T=y0(10:end,:)';    %测试集输出数据   
N=size(test_P,2);  

[train_p,ps_input]=mapminmax(train_P,0,1);%训练集输入数据归一化
test_p=mapminmax('apply',test_P,ps_input);%测试集输入数据归一化
[train_t,ps_output]=mapminmax(train_T,0,1);%训练集输出数据归一化
np=10;%种群大小
ng=100;%进化代数
variatenum=S;%变量个数
len=10;%编码长度
Youx=struct('fitx',0,'X',[],'Erjinzhi',[],'chrom',[]);%记录最佳适应度值，十进制值,二进制编码，量子比特编码
for i=1:2*np
    for j=1:variatenum*len
        x(i,j)=1/sqrt(2);%初始化种群
    end
end
Erjinzhi=Celiang(x);%对种群实施一次测量，得到二进制编码
Field=[repmat(len,1,S);repmat([-1;1],1,S);repmat([1;0;1;1],1,S)];%区域描述器
%   FieldD = [len; lb; ub; code; scale; lbin; ubin]
%   len是包含在Chrom中的每个子串的长度，注意sum(len)=size(Chrom, 2)
%   lb和ub分别是每个变量的下界和上界。
%   code指明子串是怎样编码的，1为标准的二进制编码， 0为格雷编码
%   scale指明每个子串所使用的刻度，0表示算术刻度，1表示对数刻度。
%   lbin和ubin指明表示范围中是否包含边界。0表示不包含，1表示包含边界

X=bs2rv(Erjinzhi,Field);%得到种群个体的十进制
fitx=objfun(X,train_p,train_t);%计算适应度值
[best.fitx bestindex]=max(fitx);%找出最大值
best.Erjinzhi=Erjinzhi(bestindex,:);%找出最大适应度对应的二进制
best.X=X(bestindex,:);%找出最大适应度对应的十进制
best.x=x((2*bestindex-1):(2*bestindex),:);%找出最大适应度对应的量子比特码
Youfitx(1)=best.fitx;
for k=2:ng
    Erjinzhi=Celiang(x);%对种群实施一次测量
    X=bs2rv(Erjinzhi,Field);%得到种群个体的十进制
    fitx=objfun(X,train_p,train_t);%计算适应度
    x=XuanZhuan(x,fitx,best.Erjinzhi,Erjinzhi,best.fitx);%量子旋转门
    [newbestfitx,newbestindex]=max(fitx);%找出此代最佳
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
    xlabel('遗传代数')
    ylabel('误差')
    title('进化过程')
end
R2=0
%while R2<0.9
%创建神经网络
net=feedforwardnet(S2);
net=configure(net,train_p,train_t);
net.layers{2}.transferFcn='logsig';
%设置神经网络参数
net.trainparam.epochs=1000;%迭代次数
net.trainparam.lr=0.1;%学习率
net.trainparam.goal=0.01;%训练误差率
%赋值给神经网络
 w1num=S1*S2;                             %输入层到隐含层的权值个数
 w2num=S2*S3;                             %隐含层到输出层的权值个数
 w1=best.X(1:w1num);                           %输入层到隐含层权值
 b1=best.X(w1num+1:w1num+S2);                  %隐含层神经元阈值
 w2=best.X(w1num+S2+1:w1num+S2+w2num);         %隐含层到输出层权值
 b2=best.X(w1num+S2+w2num+1:w1num+S2+w2num+S3);%输出层神经元阈值
 net.iw{1,1}=reshape(w1,S2,S1);
 net.lw{2,1}=reshape(w2,S3,S2);
 net.b{1}=reshape(b1,S2,1);
 net.b{2}=reshape(b2,S3,1);

%利用新的权值和阈值进行训练
net=train(net,train_p,train_t);
%仿真测试
s_ga=sim(net,test_p);%遗传优化后的仿真结果
%反归一化
sim_T=mapminmax('reverse',s_ga,ps_output);
%性能评价
err=norm(sim_T-test_T);
%决定系数R^2
%R2=(N*sum(sim_T.*test_T)-sum(sim_T)*sum(test_T))^2/((N*sum((sim_T).^2)-(sum(sim_T))^2)*(N*sum((test_T).^2)-(sum(test_T))^2));
figure(2)
plot(1:N,test_T,'b:*',1:N,sim_T,'r-o')
legend('真实值','预测值')
xlabel('预测样本')
ylabel('预测值')
string={'测试集预测结果对比';['R^2=',num2str(R2)]};
title(string)
%end
%结果对比
result=[test_T',sim_T']
function Erjinzhi=Celiang(chrom)%测量函数，对种群实施测量，得到二进制编码
   [m,n]=size(chrom);
   for i=1:m/2
       for j=1:n
           Rand=rand(1);%产生一个随机数
           if Rand>(chrom(2*i-1,j)^2)%如果产生的随机数大于α的平方
               Erjinzhi(i,j)=1;
           else
               Erjinzhi(i,j)=0;
           end
       end
   end
end
function chrom=XuanZhuan(chrom,fitx,Youx,Erjinzhi,bestfitx)%量子旋转门函数
   sizepop=size(chrom,1)/2;%种群大小
   lenchrom=size(Erjinzhi,2);%量子比特编码长度
   for i=1:sizepop
       for j=1:lenchrom
           arfa=chrom(2*i-1,j);  %α
           beita=chrom(2*i,j);   %β
           x=Erjinzhi(i,j);      %将测量完的种群赋值给x
           b=Youx(j);
           if((x==0)&(b==1))||((x==1)&(b==1))
               delta=0;          %旋转角的大小
               s=0;              %旋转角的符号，旋转方向
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
           fai=s*delta;      %旋转角
           U=[cos(fai) -sin(fai)
              sin(fai)  cos(fai)];%量子旋转门
          y=U*[arfa;beita];   %更新后的量子位
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
    S1=size(train_p,1);        %输入神经元个数
    S2=2*size(train_p,1)+1;    %隐层神经元个数
    S3=size(train_t,1);        %输出神经元个数
   
    %新建BP网络
    %net=newff(train_p,train_t,S1)
    net=feedforwardnet(S2);
    net=configure(net,train_p,train_t);
    net.layers{2}.transferFcn='logsig';
    
    %设置神经网络参数
    net.trainparam.epochs=1000;%迭代次数
    net.trainparam.goal=0.01;%训练误差率
    net.trainparam.lr=0.1;%学习率
    net.trainparam.show=10;%现实频率，这里设置为没训练10次显示一次
    net.trainparam.showWindow=0;%是n个周期后显示一下收敛曲线的变化
    
    %设置初始权值和阈值
    w1num=S1*S2;                             %输入层到隐含层的权值个数
    w2num=S2*S3;                             %隐含层到输出层的权值个数
    w1=x(1:w1num);                           %输入层到隐含层权值
    b1=x(w1num+1:w1num+S2);                  %隐含层神经元阈值
    w2=x(w1num+S2+1:w1num+S2+w2num);         %隐含层到输出层权值
    b2=x(w1num+S2+w2num+1:w1num+S2+w2num+S3);%输出层神经元阈值
    net.iw{1,1}=reshape(w1,S2,S1);
    net.lw{2,1}=reshape(w2,S3,S2);
    net.b{1}=reshape(b1,S2,1);
    net.b{2}=reshape(b2,S3,1);
    
    %训练神经网络
    net=train(net,train_p,train_t);

    %测试网络
    Y=sim(net,train_p);
    err=1/max(max(abs(Y-train_t)));
    %err=1./norm(Y-train_t);
end
