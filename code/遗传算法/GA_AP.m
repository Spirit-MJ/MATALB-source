close all
clear all
clc
C=[6 7 11 2
   4 5 9  8
   3 1 10 4
   5 9 8  2]
np=1000;%��Ⱥ��С
pc=0.9;%�������
pm=0.1;%�������
ng=100;%��������
ggap=0.8;%��Ⱥ����
len=length(C)*length(C);%���볤��
x=round(rand(np,len));%��ʼ����Ⱥ
x(1,:)=[1 zeros(1,4) 1 zeros(1,4) 1 zeros(1,4) 1];%������ʼ��
for k=1:ng
    fx=fitness(C,x);
    Objv=fx;
    [preObjV(k),index]=max(Objv);%��¼ÿ��������
    Youx=x(index,:);
    figure(1)%�Ż�����ͼ
    if k>=2
       line([k-1,k],[preObjV(k-1),preObjV(k)]);
       xlabel('��������');
       ylabel('��ֵ');
       title('��������');
    end
    nx=XuanZe(x,ggap,np,fx);%ѡ�����
    nx=JiaoCha(nx,pc);%�������
    nx=BianYi(nx,pm);%�������
    nx=Reins(x,nx,Objv);%��ȫ��Ⱥ
    x=nx;
end
You=reshape(Youx,length(C),length(C))'
preObjV(k)
function fx=fitness(C,x)
    for i=1:size(x,1)
        x0=reshape(x(i,:),length(C),length(C))';   
        if sum(x0,1)==ones(1,size(C,2))&sum(x0,2)==ones(size(C,1),1)
           fx(i)=1/(sum(sum(x0.*C)));
        else
           fx(i)=0; 
        end
    end
end
function nx=XuanZe(x,GGAP,np,fx)%ѡ��
       GGAP1=floor(GGAP*np);  %ѡ�����µĸ�����
%        cumfx=cumsum(fx)./sum(fx);%���̶�
%        for j=1:GGAP1
%            sita=rand(1);
%            for i=1:length(fx)
%                if cumfx(i)==0
%                   nx(j,:)=x(randperm(np,1),:);
%                else sita<=cumfx(i)
%                    nx(j,:)=x(i,:);
%                    break;
%                end
%            end
%        end
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
            nx(i,Rand)=~nx(i,Rand);%�����������λ��
        end
    end
end
function nx=Reins(x,nx,Objv)%��ȫ��Ⱥ
    NIND=size(x,1);%��ʼ��Ⱥ��С
    NSel=size(nx,1);%ѡ���������Ⱥ�Ĵ�С
    [~,index]=sort(Objv,'descend');%���Ŵ�����ǰ�ľ�������
    nx=[x(index(1:NIND-NSel),:);nx];%���Ŵ�����ǰ���ŵ�һЩ���屣������
end