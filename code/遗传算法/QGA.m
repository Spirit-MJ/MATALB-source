close all
clear
clc
np=40;%��Ⱥ��С
ng=200;%��������
variatenum=2;%��������
len=20;%���볤��
for i=1:2*np
    for j=1:variatenum*len
        x(i,j)=1/sqrt(2);%��ʼ����Ⱥ
    end
end
Erjinzhi=Celiang(x);%����Ⱥʵʩһ�β������õ������Ʊ���
Field=[repmat(len,1,variatenum);[-3.0 4.1;12.1 5.8];repmat([1;0;1;1],1,variatenum)];%����������
%   FieldD = [len; lb; ub; code; scale; lbin; ubin]
%   len�ǰ�����Chrom�е�ÿ���Ӵ��ĳ��ȣ�ע��sum(len)=size(Chrom, 2)
%   lb��ub�ֱ���ÿ���������½���Ͻ硣
%   codeָ���Ӵ�����������ģ�1Ϊ��׼�Ķ����Ʊ��룬 0Ϊ���ױ���
%   scaleָ��ÿ���Ӵ���ʹ�õĿ̶ȣ�0��ʾ�����̶ȣ�1��ʾ�����̶ȡ�
%   lbin��ubinָ����ʾ��Χ���Ƿ�����߽硣0��ʾ��������1��ʾ�����߽�
X=bs2rv(Erjinzhi,Field);%�õ���Ⱥ�����ʮ����
fitx=fitness(Erjinzhi,Field);%������Ӧ��ֵ
[best.fitx bestindex]=max(fitx);%�ҳ����ֵ
best.Erjinzhi=Erjinzhi(bestindex,:);%�ҳ������Ӧ�ȶ�Ӧ�Ķ�����
best.X=X(bestindex,:);%�ҳ������Ӧ�ȶ�Ӧ��ʮ����
best.x=x((2*bestindex-1):(2*bestindex),:);%�ҳ������Ӧ�ȶ�Ӧ�����ӱ�����
Youfitx(1)=best.fitx;
for k=2:ng
    Erjinzhi=Celiang(x);%����Ⱥʵʩһ�β���
    X=bs2rv(Erjinzhi,Field);%�õ���Ⱥ�����ʮ����
    fitx=fitness(Erjinzhi,Field);%������Ӧ��
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
    line([k-1,k],[Youfitx(k-1),Youfitx(k)]);
    grid on
    xlabel('��������')
    ylabel('����ֵ')
    title('��������')
end
best.fitx
best.X
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
   sita=0.01;
   for i=1:sizepop
       for j=1:lenchrom
           arfa=chrom(2*i-1,j);  %��
           beita=chrom(2*i,j);   %��
           x=Erjinzhi(i,j);      %�����������Ⱥ��ֵ��x
           b=Youx(j);
           if((x==0)&(b==0))||((x==1)&(b==1))%������ŵĻ�����ÿ�������ÿλ������ͬ�������仯
               delta=0;          %��ת�ǵĴ�С
               s=0;              %��ת�ǵķ��ţ���ת����
           elseif(x==0)&(b==1)&(fitx(i)<bestfitx)
               delta=sita*pi;
               if arfa*beita>0
                   s=1;
               elseif arfa*beita<0
                   s=-1;
               elseif arfa==0
                   s=0;
               elseif beita==0
                   s=sign(randn);%���һ��������ת
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
           fai=s*delta;      %��ת��
           U=[cos(fai) -sin(fai)
              sin(fai)  cos(fai)];%������ת��
          y=U*[arfa;beita];   %���º������λ
          chrom(2*i-1,j)=y(1);
          chrom(2*i,j)=y(2);
       end
   end        
end
function fitx=fitness(Erjinzhi,Field)%��Ӧ�Ⱥ���
    sizepop=size(Erjinzhi,1);  %��Ⱥ��С
    for i=1:sizepop
        fitx(i)=objfunction(bs2rv(Erjinzhi(i,:),Field));
    end
end
function fx=objfunction(x)
    fx=sin(4*pi*x(1))*x(1)+sin(20*pi*x(2))*x(2);
end