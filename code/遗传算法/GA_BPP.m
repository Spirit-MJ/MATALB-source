close all
clear all
clc
Weight=[71 34 82 23 1 88 12 57 10 68 5 33 37 69 98 24 26 83 16 26 18 43 52 71 22 65 68 8 40 40 24 72 16 34 10 19 28 13 34 98 29 31 79 33 60 74 44 56 54 17];
Value=[26 59 30 19 66 85 94 8 3 44 5 1 41 82 76 1 12 81 73 32 74 54 62 41 19 10 65 53 56 53 70 66 58 22 72 33 96 88 68 45 44 61 78 78 6 66 11 59 83 48];
WeightLimit=300;%����
np=1000;%��Ⱥ��С
pc=0.9;%�������
pm=0.3;%�������
ng=300;%��������
ggap=0.8;%��Ⱥ����
len=length(Weight);%���볤��
x=round(rand(np,len));%��ʼ����Ⱥ
x(1,:)=zeros(1,len);%������ʼ��
for k=1:ng
    fx=fitness(Weight,Value,x,WeightLimit);
    Objv=fx;
    preObjV(k)=max(Objv);%��¼ÿ��������
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
preObjV(ng)
function fx=fitness(Weight,Value,x,WeightLimit)
    for i=1:size(x,1)
        fx(i)=sum(x(i,:).*Value);
        if sum(x(i,:).*Weight)>=WeightLimit
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