close all
clear all
clc
np=200;%������Ŀ
mp=200;%��Ⱥ��Ŀ
ng=20;%��������
ggap=0.9;%��Ⱥ����
variatenum=2;%��������
len=20;%ÿ�������ı��볤��
k0=0;%��ʼ��������
k=0;
variatex=-3:0.01:12.1;
variatey=4.1:0.01:5.8;
for i=1:length(variatex)
    for j=1:length(variatey)
        z(i,j)=21.5+variatex(i)*sin(4*pi*variatex(i))+variatey(j)*sin(20*pi*variatey(j));
    end
end
figure(1)
mesh(variatex',variatey',z')
FieldD=[repmat(len,1,variatenum);[-3,4.1;12.1,5.8];repmat([1;0;1;1],1,variatenum)];%����������
%   FieldD = [len; lb; ub; code; scale; lbin; ubin]
%   len�ǰ�����Chrom�е�ÿ���Ӵ��ĳ��ȣ�ע��sum(len)=size(Chrom, 2)
%   lb��ub�ֱ���ÿ���������½���Ͻ硣
%   codeָ���Ӵ�����������ģ�1Ϊ��׼�Ķ����Ʊ��룬 0Ϊ���ױ���
%   scaleָ��ÿ���Ӵ���ʹ�õĿ̶ȣ�0��ʾ�����̶ȣ�1��ʾ�����̶ȡ�
%   lbin��ubinָ����ʾ��Χ���Ƿ�����߽硣0��ʾ��������1��ʾ�����߽�
for i=1:mp
    x{i}=round(rand(np,len*variatenum));%������ʼ��Ⱥ
end
%x{1}(1,:)=zeros(1,len);%��ʼ��
pc=0.7+(0.9-0.7)*rand(mp,1);%��[0.7,0.9]��������������������
pm=0.05+(0.25-0.05)*rand(mp,1);%��[0.05,0.25]��������������������
maxfx=0;
Youfx=zeros(mp,1);%��¼������Ⱥ
Youx=zeros(mp,len*variatenum);%��¼������Ⱥ�ı���
for i=1:mp
    fx{i}=objectfunction(bs2rv(x{i},FieldD));%�����ʼ��ȺĿ�꺯��ֵ
end
while k0<ng
    k=k+1;
    for j=1:mp
        fitx{j}=fx{j};%�������Ⱥ����Ӧ��
        nx{j}=XuanZe(x{j},ggap,np,fitx{j});%ѡ�����
        nx{j}=JiaoCha(nx{j},pc(j));%�������
        nx{j}=BianYi(nx{j},pm(j));%�������
        objvsel=objectfunction(bs2rv(nx{j},FieldD));%�����Ӵ�Ŀ�꺯��ֵ
        [x{j},fx{j}]=Reins(x{j},nx{j},fitx{j},fx{j},objvsel);%��ȫ��Ⱥ
    end
    [x,fx]=Yimin(x,fx,mp);%�������
    [Youfx,Youx]=PeopleSelect(x,fx,Youfx,Youx,mp);%�˹�ѡ�񾫻���Ⱥ
    Zuiyoufx(k)=max(Youfx);%�ҳ�������Ⱥ�����ŵĸ���
    if Zuiyoufx(k)>maxfx %�жϵ�ǰ�Ż�ֵ�Ƿ��ǰһ����
        maxfx=Zuiyoufx(k);
        k0=0;
    else
        k0=k0+1;
    end  
    if k>=2
       figure(2)
       line([k-1,k],[Zuiyoufx(k-1),Zuiyoufx(k)]);
       xlabel('��������')
       ylabel('���Ž�仯')
       title('��������')
    end
end
function fx=objectfunction(X)
  hang=size(X,1);
  for i=1:hang
      fx(i,1)=21.5+X(i,1)*sin(4*pi*X(i,1))+X(i,2)*sin(20*pi*X(i,2));
  end
end
function nx=XuanZe(x,GGAP,np,fitx)%ѡ��
       GGAP1=floor(GGAP*np);  %ѡ�����µĸ�����      
       [~,temp]=sort(fitx,'descend');%������Ӧ�Ƚ�������
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
            nx(i,Rand)=~nx(i,Rand);%�����������λ��
        end
    end
end
function [x,objv]=Reins(x,nx,fitx,fx,objvsel)%��ȫ��Ⱥ
    NIND=size(x,1);%��ʼ��Ⱥ��С
    NSel=size(nx,1);%ѡ���������Ⱥ�Ĵ�С
    [~,index]=sort(fitx,'descend');%������Ӧ�Ƚ�������
    %[~,index]=sort(fx,'ascend');%��������
    x=[nx;x(index(1:NIND-NSel),:)];%���Ŵ�����ǰ���ŵ�һЩ���屣������
    objv=[fx(index(1:NIND-NSel));objvsel];%�Ӵ�����ֵ 
end
function [chrom,objv]=Yimin(chrom,objv,mp)%�������
for i=1:mp
    [~,Youi]=max(objv{i});%�ҳ���i����Ⱥ�����ŵĸ���
    mubiaoi=i+1;%Ŀ����Ⱥ(�������)
    if mubiaoi>mp
        mubiaoi=1;
    end
    [~,Chai]=min(objv{mubiaoi});%�ҳ�Ŀ����Ⱥ�����ӵĸ���
    chrom{mubiaoi}(Chai,:)=chrom{i}(Youi,:);%Ŀ����Ⱥ���Ӹ����滻ΪԴ��Ⱥ�����ӵĸ���
    objv{mubiaoi}(Chai)=objv{i}(Youi);
end
end
function [maxobjv,maxchrom]=PeopleSelect(x,fx,maxobjv,maxchrom,mp)%�˹�ѡ������
for i=1:mp
    [Youo,Youi]=max(fx{i});%�ҳ���i��Ⱥ�����ŵĸ���
    if Youo>maxobjv(i)
       maxobjv(i)=Youo;%��¼����Ⱥ�ľ�������
       maxchrom(i,:)=x{i}(Youi,:);
    end
end
end
