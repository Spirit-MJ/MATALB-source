close all
clear
clc
%load('C:\Users\idiots\Desktop\MATLAB智能算法30个案例分析源代码\chapter11-book\scheduleData')
%  元胞数组  Jm  行代表工件 列代表工序 每行每列代表每个工件在每个工序下可用的机器集合
%  JmNumber  机器个数
%  元胞数组  T 行代表工件 列代表工序，每行每列代表每个工件在每个工序下使用机器的时间的集合 
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
np=500;        %种群大小
ng=200;        %进化代数
GGAP=0.95;       %代沟
pc=0.8;       %交叉概率
pm=0.1;        %变异概率
gen=1;          %代计数器2
PNumber=size(Jm,1);          %PNumber 工件个数 
trace=zeros(2, ng);      %寻优结果的初始值
por=everyGongxu(Jm);         %求出每个工件的工序数
WNumber=sum(por);            %工序总个数
Number=zeros(1,PNumber);     %1行， PNumber工件个数列，储存每个工件的工序个数
for i=1:PNumber
    Number(i)=por(i);       %每个工件的 MNumber 工序个数
end
x=zeros(np,2*sum(Number)); %初始化种群
for i=1:np      %产生初始种群
    WPNumberTemp=Number;     %将每个工件需要的工序数赋值给 WPNumberTemp
    for j=1:WNumber
        val=randperm(PNumber,1);    %随机产生一个工件
        while WPNumberTemp(val)==0   %防止一个工件的工序完成后产生多余的工序
            val=randperm(PNumber,1);
        end
        x(i,j)=val;          %初始化种群
        WPNumberTemp(val)=WPNumberTemp(val)-1;  %这个工件的工序就少一 
        Temp=Jm{val,Number(val)-WPNumberTemp(val)}; %找到该工件该工序下的可用机器集合
        SizeTemp=length(Temp);                  %得到可用机器集合的数量
        x(i,j+WNumber)=randperm(SizeTemp,1);  %随机选择该工件该工序的可用机器
    end
end
%遗传操作
while gen<ng
    [PVal,fx,P,S]=cal(x,JmNumber,T,Jm,por);%记算出每个个体的Objv(总完成时间)Pval(最优解的完成时间情况)P(最佳调度工序)S(最佳染色体)
    fitx=1./fx; %分配适应度值 
    nx=XuanZe(x,GGAP,fitx);%选择操作
    nx=Jiaocha(nx,pc,PNumber);%交叉操作      
    nx=Bianyi(nx,pm,Jm,T);%变异操作 
    nx=Reverse(nx);%进化逆转操作
    [newPVal,newfx,newP,newS]=cal(nx,JmNumber,T,Jm,por);%计算新种群的适应度值
    x=Reins(x,nx,fitx);%重新插入新种群  
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
    %描绘解的变化
    figure(1)
    if gen>=2
       draw1=line([gen-1,gen],[trace(1, gen-1),trace(1, gen)]);
       hold on;
       draw2=line([gen-1,gen],[trace(2, gen-1),trace(2, gen)]);
       grid on
       set(draw1,'color',[0 0 1]);
       set(draw2,'color',[1 0 0]);
       xlabel('进化代数')
       ylabel('最优工序的时间')
       title('进化过程')
       legend('最优解的变化','种群均值的变化')
    end
    gen=gen+1;%代计数器增加       
end
PVal=Val1 %最佳工序时间
P=Val2  %最佳工序 
S=STemp; %最佳调度基因含机器基因
MP=S(1,sum(Number)+1:2*sum(Number));%取出机器的基因
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
end%求出每个工件的工序数
function [PVal,ObjV,P,S]=cal(Chrom,JmNumber,T,Jm,por)
% 功能说明：       根据基因群,计算出个群中每个个体的调度工序时间，
%                 保存最小时间的调度工序和调度工序时间
% 输入参数：
%       Chrom     为基因种群  
%       T         为各工件各工序使用的时间 
%       Jm        为各工件各工序使用的机器 
%       por       为各工件的工序数
% 输出参数:
%       PVal      为最佳调度工序时间 
%       P         为最佳输出的调度工序 
%       ObjV      为群中每个个体的调度工序时间
%       S         为最佳输出的调度基因
NIND=size(Chrom,1);      %种群的大小
ObjV=zeros(NIND,1);      %记录每个可行解的完成时间
PNumber=size(Jm,1);   %PNumber 工件个数 MNumber  工序个数
for i=1:NIND  
    S=Chrom(i,:);             %取一个个体
    P=calp(S,PNumber);        %根据基因，计算调度工序   
    PVal=caltime(S,P,JmNumber,T,Jm,por);   %根据调度工序，计算出调度工序时间  
    TVal=max(max(PVal));         %取完成时间  
    ObjV(i,1)=TVal;              %保存完成时间
    %初始化
    if i==1
        Val1=PVal;               
        Val2=P;
        MinVal=ObjV(i,1);
        STemp=S;
    end
    %记录 最小的调度工序时间、最佳调度工序时间 最佳输出的调度工序
    if MinVal>ObjV(i,1)
        Val1=PVal;
        Val2=P;
        MinVal=ObjV(i,1);
        STemp=S;
    end   
end 
%最佳调度工序时间 最佳输出的调度工序
 PVal=Val1;
 P=Val2;
 S=STemp;
end%计算每个解的时间等
function PVal=caltime(S,P,JmNumber,T,Jm,por)%计算出调度工序时间
% 功能说明：    根据调度工序,计算出调度工序时间
% 输入参数：
%        P     为调度工序  
%        JmNumber    为机器个数
%        T     为各工件各工序的加工时间 
%        Jm    为各工件各工序使用的机器 
% 输出参数:
%        PVal  为调度工序开始加工时间及完成时间

PNumber=size(Jm,1);    %PNumber 工件个数 
M=S(1,sum(por)+1:2*sum(por));   %取机器基因，取基因的后一半即使用的机器
WNumber=length(P);             %工序总个数
%初始化
TM=zeros(1,JmNumber);          %1行，机器个数列
TP=zeros(1,PNumber);           %1行，工件个数列
PVal=zeros(2,WNumber);         %2行，工序总个数列
for i=1:WNumber                %计算调度工序时间
    val=P(1,i);               %取机器号
    a=(mod(val,100));          %工件对应的工序
    b=((val-a)/100);           %工件
    Temp=Jm{b,a};              %找出该工件在该工序下的可用机器集合   
    m=Temp(M(i));            %解码基因得到该工件在该工序下使用的机器   
    
    Temp=T{b,a};               %找出该工件在该工序下的可用机器的加工时间集合
    t=Temp(M(i));            %解码基因得到该工件在该工序下使用的机器的时间
%取机器加工本工序的开始时间和前面一道工序的完成时间
    TMval=TM(1,m);
    TPval=TP(1,b); 
    if TMval>TPval    %机器加工本工件的开始时间大于工件前面一道工序的完成时间
       val=TMval;     %取机器加工本工序的开始时间  
    else
       val=TPval;     %取前面一道工序的完成时间
    end   
    %计算时间
    PVal(1,i)=val;    %开始时间
    PVal(2,i)=val+t;  %完成时间
    %记录本次工序的机器时间和工序时间
    TM(1,m)=PVal(2,i);             %更新加工本工序的开始时间
    TP(1,b)=PVal(2,i);             %更新上一道工序完成时间
end
end
function P=calp(S,PNumber)%解码调度工序
% 功能说明：          根据基因S,计算调度工序P
% 输入参数：
%        S           为基因  
%        PNumber     为工件个数 
% 输出参数: 
%        P           为输出的调度工序 
WNumber=length(S)/2;    %总工序数
S=S(1,1:WNumber);     %取工序基因，取基因的一半
%初始化
temp=zeros(1,PNumber);  %一行，工件个数 列
P=zeros(1,WNumber);     %一行，总工序数 列 
for i=1: WNumber       %解码生成调度工序
  temp(S(i))=temp(S(i))+1;    %工序加+1
  P(i)=S(i)*100+temp(S(i));   %生成形如201的数，2代表工件，1代表第一道工序
end
end
function nx=Jiaocha(Chrom,pc,PNumber)

% Chrom=[1 3 2 3 1 2 1 3 2; 
%     1 1 2 3 3 1 2 3 2;
%     1 3 2 3 2 2 1 3 1;
%     1 3 3 3 1 2 1 2 2;
% ]; 
%   XOVR=0.7;
[NIND,WNumber]=size(Chrom);  %选择操作之后种群的大小为 NIND，基因长度为WNumber
WNumber=WNumber/2;           %总工序个数
nx=Chrom;              %初始化新种群
%PNumber=size(Jm,1);  %PNumber 工件个数 
Number=zeros(1,PNumber);     %一行，工件个数列
for i=1:PNumber
  Number(i)=1;               %记录每个工件的工序
end
%随机选择交叉个体(洗牌交叉)
SelNum=randperm(NIND);       %把种群的顺序打乱  
Num=floor(NIND/2);           %交叉个体配对数，如果多余1个不交叉
for i=1:2:Num
    if pc>rand             %如果满足交叉概率
        Pos=randperm(WNumber,1);%随机产生一个交叉位置
        while Pos==1
            Pos=randperm(WNumber,1);%防止交叉的位置在第一个基因处
        end
        %取两交叉的个体
        S1=Chrom(SelNum(i),1:WNumber);    %取出前总工序的基因
        S2=Chrom(SelNum(i+1),1:WNumber); 
        %初始化新的个体
        S11=S2;
        S22=S1;      
        %交叉开始      
        S11(1:Pos)=S1(1:Pos);      
        S22(1:Pos)=S2(1:Pos);        
        %比较S11相对S1,S22相对S2多余和缺失的基因
        S3=S11;    %交叉后的第一个染色体
        S4=S1;     %交叉之前的第一个染色体
        S5=S22;    %交叉后的第二个染色体
        S6=S2;     %交叉之前的第二个染色体
        for j=1:WNumber         
           Pos1=find(S4==S3(j),1);  %找出交叉后的第一个染色体的第i位基因在原来交叉的染色体的位置
           Pos2=find(S6==S5(j),1);  %找出交叉后的第二个染色体的第i位基因在原来交叉的染色体的位置
           if Pos1>0                %如果交叉后染色体上的基因在原来进行交叉的染色体上找得到的话
               S3(j)=0;             %就把交叉后的染色体的基因赋值为0
               S4(Pos1)=0;          %把原来
           end                         
           if Pos2>0
               S5(j)=0;
               S6(Pos2)=0;
           end
        end
        for j=1:WNumber          
          if S3(j)~=0                 %多余的基因          
            Pos1=find(S11==S3(j),1);  %找出不等于零的第一个基因位置      
            Pos2=find(S4,1);          %查找缺失的基因
            S11(Pos1)=S4(Pos2);       %用缺失的基因修补多余的基因
            S4(Pos2)=0;       
          end 
          if S5(j)~=0              
            Pos1=find(S22==S5(j),1); 
            Pos2=find(S6,1);           
            S22(Pos1)=S6(Pos2);
            S6(Pos2)=0;          
          end  
        end                         
        % 保存交叉前的机器 基因
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
        %生成新的种群
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
[NIND,WNumber]=size(Chrom);%选择交叉操作后种群的大小NIND 基因长度WNumber
WNumber=WNumber/2; %总工序数
ChromNew=Chrom;    %初始化
[PNumber MNumber]=size(Jm); %PNumber 工件个数 MNumber  工序个数
Number=zeros(1,PNumber);    %一行，工件个数列
for i=1:PNumber
  Number(i)=1;
end
for i=1:NIND              
    S=Chrom(i,:);    %取出一个个体          
       WPNumberTemp=Number;       
       for j=1:WNumber        
          JMTemp=Jm{S(j), WPNumberTemp(S(j))};    %第j个基因代表的工件和工序
          SizeTemp=length(JMTemp);    %能使用的机器数量    
          if MUTR>rand %是否变异
                %选择机器（随机选择）
                %S(j+WNumber)=randperm(SizeTemp,1); 
                %选择机器（加工时间少的选择几率大）
                if SizeTemp==1      
                       S(j+WNumber)=1;      %只能选择一个机器
                else
                    S(j+WNumber)=roulette(T{S(j),WPNumberTemp(S(j))});
                end
          end
            WPNumberTemp(S(j))=WPNumberTemp(S(j))+1;  %工序加1
        end         
    ChromNew(i,:)=S;
end
end%变异操作
function bit=roulette(S_T)%轮盘赌操作
T=1./S_T;
cumfx=cumsum(T)./sum(T);%轮盘赌
sita=rand(1);
for i=1:length(S_T)
    if sita<=cumfx(i)
        bit=i;
        break;
    end
end
end
function nx=XuanZe(x,GGAP,fx)%选择
       np=size(x,1);
       GGAP1=floor(GGAP*np);  %选择留下的个体数
       [~,temp]=sort(fx,'descend');%按照适应度降序排列
       for j=1:GGAP1
           nx(j,:)=x(temp(j),:);%优胜略汰，将适应度大的个体保留下来
       end
end
function SelCh=Reverse(SelCh)%进化逆转
    [row,col]=size(SelCh);   %变异操作后种群的大小
     col=col/2;              %总工序数
     SelCh1=SelCh;
    for i=1:row
        r=randsrc(1,2,1:col);%产生两个随机数
        mininverse=min(r);   %较小的数
        maxinverse=max(r);   %较大的数
        SelCh1(i,mininverse:maxinverse)=SelCh1(i,maxinverse:-1:mininverse);%逆转操作,工序基因对换位置
        SelCh1(i,col+mininverse:col+maxinverse)=SelCh1(i,col+maxinverse:-1:col+mininverse);%逆转操作,机器基因对换位置
    end
end
function nx=Reins(x,nx,Objv)%补全种群
    NIND=size(x,1);%初始种群大小
    NSel=size(nx,1);%选择操作后种群的大小
    [~,index]=sort(Objv,'descend');%将遗传操作前的适应度进行降序排序
    nx=[x(index(1:NIND-NSel),:);nx];%将遗传操作前最优的一些个体保留下来
end
