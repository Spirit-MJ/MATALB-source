close all
clear
clc
np=40;%种群大小
ng=200;%进化代数
variatenum=2;%变量个数
len=20;%编码长度
for i=1:2*np
    for j=1:variatenum*len
        x(i,j)=1/sqrt(2);%初始化种群
    end
end
Erjinzhi=Celiang(x);%对种群实施一次测量，得到二进制编码
Field=[repmat(len,1,variatenum);[-3.0 4.1;12.1 5.8];repmat([1;0;1;1],1,variatenum)];%区域描述器
%   FieldD = [len; lb; ub; code; scale; lbin; ubin]
%   len是包含在Chrom中的每个子串的长度，注意sum(len)=size(Chrom, 2)
%   lb和ub分别是每个变量的下界和上界。
%   code指明子串是怎样编码的，1为标准的二进制编码， 0为格雷编码
%   scale指明每个子串所使用的刻度，0表示算术刻度，1表示对数刻度。
%   lbin和ubin指明表示范围中是否包含边界。0表示不包含，1表示包含边界
X=bs2rv(Erjinzhi,Field);%得到种群个体的十进制
fitx=fitness(Erjinzhi,Field);%计算适应度值
[best.fitx bestindex]=max(fitx);%找出最大值
best.Erjinzhi=Erjinzhi(bestindex,:);%找出最大适应度对应的二进制
best.X=X(bestindex,:);%找出最大适应度对应的十进制
best.x=x((2*bestindex-1):(2*bestindex),:);%找出最大适应度对应的量子比特码
Youfitx(1)=best.fitx;
for k=2:ng
    Erjinzhi=Celiang(x);%对种群实施一次测量
    X=bs2rv(Erjinzhi,Field);%得到种群个体的十进制
    fitx=fitness(Erjinzhi,Field);%计算适应度
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
    line([k-1,k],[Youfitx(k-1),Youfitx(k)]);
    grid on
    xlabel('进化代数')
    ylabel('最优值')
    title('进化过程')
end
best.fitx
best.X
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
   sita=0.01;
   for i=1:sizepop
       for j=1:lenchrom
           arfa=chrom(2*i-1,j);  %α
           beita=chrom(2*i,j);   %β
           x=Erjinzhi(i,j);      %将测量完的种群赋值给x
           b=Youx(j);
           if((x==0)&(b==0))||((x==1)&(b==1))%如果最优的基因与每个个体的每位基因相同，则不做变化
               delta=0;          %旋转角的大小
               s=0;              %旋转角的符号，旋转方向
           elseif(x==0)&(b==1)&(fitx(i)<bestfitx)
               delta=sita*pi;
               if arfa*beita>0
                   s=1;
               elseif arfa*beita<0
                   s=-1;
               elseif arfa==0
                   s=0;
               elseif beita==0
                   s=sign(randn);%随机一个方向旋转
               end
           elseif (x==0)&(b==1)&(fitx(i)>=bestfitx)
               delta=sita*pi;
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
               delta=sita*pi;
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
               delta=sita*pi;
               if arfa*beita>0
                   s=1;
               elseif arfa*beita<0
                   s=-1;
               elseif arfa==0
                   s=0;
               elseif beita==0
                   s=sign(randn);
               end
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
function fitx=fitness(Erjinzhi,Field)%适应度函数
    sizepop=size(Erjinzhi,1);  %种群大小
    for i=1:sizepop
        fitx(i)=objfunction(bs2rv(Erjinzhi(i,:),Field));
    end
end
function fx=objfunction(x)
    fx=sin(4*pi*x(1))*x(1)+sin(20*pi*x(2))*x(2);
end