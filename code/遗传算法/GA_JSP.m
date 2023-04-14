close all
clear
clc
%load('C:\Users\idiots\Desktop\MATLAB�����㷨30����������Դ����\chapter11-book\scheduleData')
%  Ԫ������  Jm  �д����� �д����� ÿ��ÿ�д���ÿ��������ÿ�������¿��õĻ�������
%  JmNumber  ��������
%  Ԫ������  T �д����� �д�����ÿ��ÿ�д���ÿ��������ÿ��������ʹ�û�����ʱ��ļ��� 
Jm={[3 10] [1] [2] [4 7] [6 8] [5]
    [2] [3] [5 8] [6 7] [1] [4 10]
    [3 9] [4 7] [6 8] [1] [2 10] [5]
    [4] [1 9] [3 7] [2 8] [5] [6]
    [5] [2 7] [3 10] [6 9] [1] [4 8]
    [2] [4 7] [6 9] [1] [5 8] [3]};
JmNumber=10;
T={[3 5] [10] [9] [5 4] [3 3] [10]
    [6] [8] [1 4] [5 6] [3] [3 3]
    [1 4] [5 7] [5 6] [5] [9 11] [1]
    [7] [4 3] [4 6] [3 5] [1] [3]
    [6] [10 12] [7 9] [8 8] [5] [4 7]
    [2] [4 7] [6 9] [1] [5 8] [3]};
np=500;        %��Ⱥ��С
ng=200;        %��������
GGAP=0.95;       %����
pc=0.8;       %�������
pm=0.1;        %�������
gen=1;          %��������2
PNumber=size(Jm,1);          %PNumber �������� 
trace=zeros(2, ng);      %Ѱ�Ž���ĳ�ʼֵ
por=everyGongxu(Jm);         %���ÿ�������Ĺ�����
WNumber=sum(por);            %�����ܸ���
Number=zeros(1,PNumber);     %1�У� PNumber���������У�����ÿ�������Ĺ������
for i=1:PNumber
    Number(i)=por(i);       %ÿ�������� MNumber �������
end
x=zeros(np,2*sum(Number)); %��ʼ����Ⱥ
for i=1:np      %������ʼ��Ⱥ
    WPNumberTemp=Number;     %��ÿ��������Ҫ�Ĺ�������ֵ�� WPNumberTemp
    for j=1:WNumber
        val=randperm(PNumber,1);    %�������һ������
        while WPNumberTemp(val)==0   %��ֹһ�������Ĺ�����ɺ��������Ĺ���
            val=randperm(PNumber,1);
        end
        x(i,j)=val;          %��ʼ����Ⱥ
        WPNumberTemp(val)=WPNumberTemp(val)-1;  %��������Ĺ������һ 
        Temp=Jm{val,Number(val)-WPNumberTemp(val)}; %�ҵ��ù����ù����µĿ��û�������
        SizeTemp=length(Temp);                  %�õ����û������ϵ�����
        x(i,j+WNumber)=randperm(SizeTemp,1);  %���ѡ��ù����ù���Ŀ��û���
    end
end
%�Ŵ�����
while gen<ng
    [PVal,fx,P,S]=cal(x,JmNumber,T,Jm,por);%�����ÿ�������Objv(�����ʱ��)Pval(���Ž�����ʱ�����)P(��ѵ��ȹ���)S(���Ⱦɫ��)
    fitx=1./fx; %������Ӧ��ֵ 
    nx=XuanZe(x,GGAP,fitx);%ѡ�����
    nx=Jiaocha(nx,pc,PNumber);%�������      
    nx=Bianyi(nx,pm,Jm,T);%������� 
    nx=Reverse(nx);%������ת����
    [newPVal,newfx,newP,newS]=cal(nx,JmNumber,T,Jm,por);%��������Ⱥ����Ӧ��ֵ
    x=Reins(x,nx,fitx);%���²�������Ⱥ  
    if min(fx)<min(newfx)
       trace(1, gen)=min(fx);
       trace(2, gen)=mean(fx);
       Val1=PVal;
       Val2=P;
       MinVal=trace(1,gen);
       STemp=S;
    else
       trace(1, gen)=min(newfx);
       trace(2, gen)=mean(newfx);
       Val1=newPVal;
       Val2=newP;
       MinVal=trace(1,gen);
       STemp=newS; 
    end
    %����ı仯
    figure(1)
    if gen>=2
       draw1=line([gen-1,gen],[trace(1, gen-1),trace(1, gen)]);
       hold on;
       draw2=line([gen-1,gen],[trace(2, gen-1),trace(2, gen)]);
       grid on
       set(draw1,'color',[0 0 1]);
       set(draw2,'color',[1 0 0]);
       xlabel('��������')
       ylabel('���Ź����ʱ��')
       title('��������')
       legend('���Ž�ı仯','��Ⱥ��ֵ�ı仯')
    end
    gen=gen+1;%������������       
end
PVal=Val1 %��ѹ���ʱ��
P=Val2  %��ѹ��� 
S=STemp; %��ѵ��Ȼ��򺬻�������
MP=S(1,sum(Number)+1:2*sum(Number));%ȡ�������Ļ���
function por=everyGongxu(jm)
for i=1:size(jm,1)
    temp=1;
    for j=1:size(jm,2)  
        temp=length(jm{i,j});
        if temp==0
            por(i)=j-1;
            break;
        else
            por(i)=j;
        end
    end
end
end%���ÿ�������Ĺ�����
function [PVal,ObjV,P,S]=cal(Chrom,JmNumber,T,Jm,por)
% ����˵����       ���ݻ���Ⱥ,�������Ⱥ��ÿ������ĵ��ȹ���ʱ�䣬
%                 ������Сʱ��ĵ��ȹ���͵��ȹ���ʱ��
% ���������
%       Chrom     Ϊ������Ⱥ  
%       T         Ϊ������������ʹ�õ�ʱ�� 
%       Jm        Ϊ������������ʹ�õĻ��� 
%       por       Ϊ�������Ĺ�����
% �������:
%       PVal      Ϊ��ѵ��ȹ���ʱ�� 
%       P         Ϊ�������ĵ��ȹ��� 
%       ObjV      ΪȺ��ÿ������ĵ��ȹ���ʱ��
%       S         Ϊ�������ĵ��Ȼ���
NIND=size(Chrom,1);      %��Ⱥ�Ĵ�С
ObjV=zeros(NIND,1);      %��¼ÿ�����н�����ʱ��
PNumber=size(Jm,1);   %PNumber �������� MNumber  �������
for i=1:NIND  
    S=Chrom(i,:);             %ȡһ������
    P=calp(S,PNumber);        %���ݻ��򣬼�����ȹ���   
    PVal=caltime(S,P,JmNumber,T,Jm,por);   %���ݵ��ȹ��򣬼�������ȹ���ʱ��  
    TVal=max(max(PVal));         %ȡ���ʱ��  
    ObjV(i,1)=TVal;              %�������ʱ��
    %��ʼ��
    if i==1
        Val1=PVal;               
        Val2=P;
        MinVal=ObjV(i,1);
        STemp=S;
    end
    %��¼ ��С�ĵ��ȹ���ʱ�䡢��ѵ��ȹ���ʱ�� �������ĵ��ȹ���
    if MinVal>ObjV(i,1)
        Val1=PVal;
        Val2=P;
        MinVal=ObjV(i,1);
        STemp=S;
    end   
end 
%��ѵ��ȹ���ʱ�� �������ĵ��ȹ���
 PVal=Val1;
 P=Val2;
 S=STemp;
end%����ÿ�����ʱ���
function PVal=caltime(S,P,JmNumber,T,Jm,por)%��������ȹ���ʱ��
% ����˵����    ���ݵ��ȹ���,��������ȹ���ʱ��
% ���������
%        P     Ϊ���ȹ���  
%        JmNumber    Ϊ��������
%        T     Ϊ������������ļӹ�ʱ�� 
%        Jm    Ϊ������������ʹ�õĻ��� 
% �������:
%        PVal  Ϊ���ȹ���ʼ�ӹ�ʱ�估���ʱ��

PNumber=size(Jm,1);    %PNumber �������� 
M=S(1,sum(por)+1:2*sum(por));   %ȡ��������ȡ����ĺ�һ�뼴ʹ�õĻ���
WNumber=length(P);             %�����ܸ���
%��ʼ��
TM=zeros(1,JmNumber);          %1�У�����������
TP=zeros(1,PNumber);           %1�У�����������
PVal=zeros(2,WNumber);         %2�У������ܸ�����
for i=1:WNumber                %������ȹ���ʱ��
    val=P(1,i);               %ȡ������
    a=(mod(val,100));          %������Ӧ�Ĺ���
    b=((val-a)/100);           %����
    Temp=Jm{b,a};              %�ҳ��ù����ڸù����µĿ��û�������   
    m=Temp(M(i));            %�������õ��ù����ڸù�����ʹ�õĻ���   
    
    Temp=T{b,a};               %�ҳ��ù����ڸù����µĿ��û����ļӹ�ʱ�伯��
    t=Temp(M(i));            %�������õ��ù����ڸù�����ʹ�õĻ�����ʱ��
%ȡ�����ӹ�������Ŀ�ʼʱ���ǰ��һ����������ʱ��
    TMval=TM(1,m);
    TPval=TP(1,b); 
    if TMval>TPval    %�����ӹ��������Ŀ�ʼʱ����ڹ���ǰ��һ����������ʱ��
       val=TMval;     %ȡ�����ӹ�������Ŀ�ʼʱ��  
    else
       val=TPval;     %ȡǰ��һ����������ʱ��
    end   
    %����ʱ��
    PVal(1,i)=val;    %��ʼʱ��
    PVal(2,i)=val+t;  %���ʱ��
    %��¼���ι���Ļ���ʱ��͹���ʱ��
    TM(1,m)=PVal(2,i);             %���¼ӹ�������Ŀ�ʼʱ��
    TP(1,b)=PVal(2,i);             %������һ���������ʱ��
end
end
function P=calp(S,PNumber)%������ȹ���
% ����˵����          ���ݻ���S,������ȹ���P
% ���������
%        S           Ϊ����  
%        PNumber     Ϊ�������� 
% �������: 
%        P           Ϊ����ĵ��ȹ��� 
WNumber=length(S)/2;    %�ܹ�����
S=S(1,1:WNumber);     %ȡ�������ȡ�����һ��
%��ʼ��
temp=zeros(1,PNumber);  %һ�У��������� ��
P=zeros(1,WNumber);     %һ�У��ܹ����� �� 
for i=1: WNumber       %�������ɵ��ȹ���
  temp(S(i))=temp(S(i))+1;    %�����+1
  P(i)=S(i)*100+temp(S(i));   %��������201������2��������1�����һ������
end
end
function nx=Jiaocha(Chrom,pc,PNumber)

% Chrom=[1 3 2 3 1 2 1 3 2; 
%     1 1 2 3 3 1 2 3 2;
%     1 3 2 3 2 2 1 3 1;
%     1 3 3 3 1 2 1 2 2;
% ]; 
%   XOVR=0.7;
[NIND,WNumber]=size(Chrom);  %ѡ�����֮����Ⱥ�Ĵ�СΪ NIND�����򳤶�ΪWNumber
WNumber=WNumber/2;           %�ܹ������
nx=Chrom;              %��ʼ������Ⱥ
%PNumber=size(Jm,1);  %PNumber �������� 
Number=zeros(1,PNumber);     %һ�У�����������
for i=1:PNumber
  Number(i)=1;               %��¼ÿ�������Ĺ���
end
%���ѡ�񽻲����(ϴ�ƽ���)
SelNum=randperm(NIND);       %����Ⱥ��˳�����  
Num=floor(NIND/2);           %���������������������1��������
for i=1:2:Num
    if pc>rand             %������㽻�����
        Pos=randperm(WNumber,1);%�������һ������λ��
        while Pos==1
            Pos=randperm(WNumber,1);%��ֹ�����λ���ڵ�һ������
        end
        %ȡ������ĸ���
        S1=Chrom(SelNum(i),1:WNumber);    %ȡ��ǰ�ܹ���Ļ���
        S2=Chrom(SelNum(i+1),1:WNumber); 
        %��ʼ���µĸ���
        S11=S2;
        S22=S1;      
        %���濪ʼ      
        S11(1:Pos)=S1(1:Pos);      
        S22(1:Pos)=S2(1:Pos);        
        %�Ƚ�S11���S1,S22���S2�����ȱʧ�Ļ���
        S3=S11;    %�����ĵ�һ��Ⱦɫ��
        S4=S1;     %����֮ǰ�ĵ�һ��Ⱦɫ��
        S5=S22;    %�����ĵڶ���Ⱦɫ��
        S6=S2;     %����֮ǰ�ĵڶ���Ⱦɫ��
        for j=1:WNumber         
           Pos1=find(S4==S3(j),1);  %�ҳ������ĵ�һ��Ⱦɫ��ĵ�iλ������ԭ�������Ⱦɫ���λ��
           Pos2=find(S6==S5(j),1);  %�ҳ������ĵڶ���Ⱦɫ��ĵ�iλ������ԭ�������Ⱦɫ���λ��
           if Pos1>0                %��������Ⱦɫ���ϵĻ�����ԭ�����н����Ⱦɫ�����ҵõ��Ļ�
               S3(j)=0;             %�Ͱѽ�����Ⱦɫ��Ļ���ֵΪ0
               S4(Pos1)=0;          %��ԭ��
           end                         
           if Pos2>0
               S5(j)=0;
               S6(Pos2)=0;
           end
        end
        for j=1:WNumber          
          if S3(j)~=0                 %����Ļ���          
            Pos1=find(S11==S3(j),1);  %�ҳ���������ĵ�һ������λ��      
            Pos2=find(S4,1);          %����ȱʧ�Ļ���
            S11(Pos1)=S4(Pos2);       %��ȱʧ�Ļ����޲�����Ļ���
            S4(Pos2)=0;       
          end 
          if S5(j)~=0              
            Pos1=find(S22==S5(j),1); 
            Pos2=find(S6,1);           
            S22(Pos1)=S6(Pos2);
            S6(Pos2)=0;          
          end  
        end                         
        % ���潻��ǰ�Ļ��� ����
        S1=Chrom(SelNum(i),:);
        S2=Chrom(SelNum(i+1),:); 
       
        for k=1:WNumber            
            Pos1=Find(S11(k),S1);           
            S11(WNumber+k)=S1(WNumber+Pos1);
            S1(Pos1)=0;
            
            Pos1=Find(S22(k),S2);           
            S22(WNumber+k)=S2(WNumber+Pos1);
            S2(Pos1)=0;
        end    
        %�����µ���Ⱥ
        nx(SelNum(i),:)=S11;
        nx(SelNum(i+1),:)=S22;
    end
end
end
function  Pos=Find(FindVal,S)
% S=[1 3 2 3 1 2 1 3 2];
% FindVal=3;
[m n]=size(S);
Pos=-1;
for i=1:n 
    if FindVal==S(i)
      Pos=i;
      break;
    end
end
end
function ChromNew=Bianyi(Chrom,MUTR,Jm,T)
[NIND,WNumber]=size(Chrom);%ѡ�񽻲��������Ⱥ�Ĵ�СNIND ���򳤶�WNumber
WNumber=WNumber/2; %�ܹ�����
ChromNew=Chrom;    %��ʼ��
[PNumber MNumber]=size(Jm); %PNumber �������� MNumber  �������
Number=zeros(1,PNumber);    %һ�У�����������
for i=1:PNumber
  Number(i)=1;
end
for i=1:NIND              
    S=Chrom(i,:);    %ȡ��һ������          
       WPNumberTemp=Number;       
       for j=1:WNumber        
          JMTemp=Jm{S(j), WPNumberTemp(S(j))};    %��j���������Ĺ����͹���
          SizeTemp=length(JMTemp);    %��ʹ�õĻ�������    
          if MUTR>rand %�Ƿ����
                %ѡ����������ѡ��
                %S(j+WNumber)=randperm(SizeTemp,1); 
                %ѡ��������ӹ�ʱ���ٵ�ѡ���ʴ�
                if SizeTemp==1      
                       S(j+WNumber)=1;      %ֻ��ѡ��һ������
                else
                    S(j+WNumber)=roulette(T{S(j),WPNumberTemp(S(j))});
                end
          end
            WPNumberTemp(S(j))=WPNumberTemp(S(j))+1;  %�����1
        end         
    ChromNew(i,:)=S;
end
end%�������
function bit=roulette(S_T)%���̶Ĳ���
T=1./S_T;
cumfx=cumsum(T)./sum(T);%���̶�
sita=rand(1);
for i=1:length(S_T)
    if sita<=cumfx(i)
        bit=i;
        break;
    end
end
end
function nx=XuanZe(x,GGAP,fx)%ѡ��
       np=size(x,1);
       GGAP1=floor(GGAP*np);  %ѡ�����µĸ�����
       [~,temp]=sort(fx,'descend');%������Ӧ�Ƚ�������
       for j=1:GGAP1
           nx(j,:)=x(temp(j),:);%��ʤ��̭������Ӧ�ȴ�ĸ��屣������
       end
end
function SelCh=Reverse(SelCh)%������ת
    [row,col]=size(SelCh);   %�����������Ⱥ�Ĵ�С
     col=col/2;              %�ܹ�����
     SelCh1=SelCh;
    for i=1:row
        r=randsrc(1,2,1:col);%�������������
        mininverse=min(r);   %��С����
        maxinverse=max(r);   %�ϴ����
        SelCh1(i,mininverse:maxinverse)=SelCh1(i,maxinverse:-1:mininverse);%��ת����,�������Ի�λ��
        SelCh1(i,col+mininverse:col+maxinverse)=SelCh1(i,col+maxinverse:-1:col+mininverse);%��ת����,��������Ի�λ��
    end
end
function nx=Reins(x,nx,Objv)%��ȫ��Ⱥ
    NIND=size(x,1);%��ʼ��Ⱥ��С
    NSel=size(nx,1);%ѡ���������Ⱥ�Ĵ�С
    [~,index]=sort(Objv,'descend');%���Ŵ�����ǰ����Ӧ�Ƚ��н�������
    nx=[x(index(1:NIND-NSel),:);nx];%���Ŵ�����ǰ���ŵ�һЩ���屣������
end
