close all
clear
clc
figure
axes;
set(gcf,'DoubleBuffer','on');%双缓冲设置在动画的制作中比较常用，这样设置的目的是为了防止在不断循环画动画的时候产生闪烁的现象。
S=ones(300);
S(randperm(300,1),randperm(300,1))=2;%随机地点发生一个火灾
Sk=ones(302);%防止越界
Sk(2:301,2:301)=S;%添加随机产生1 2的矩阵
% 红色表示正在燃烧(S中等于2的位置)
% 绿色表示绿树(S中等于1的位置)
C=zeros(302,302,3);
R=zeros(300);
G=zeros(300);
R(S==2)=1;%红色(S中等于2的位置)表示正在燃烧
G(S==1)=1;%绿色(S中等于1的位置)表示绿树
C(2:301,2:301,1)=R;%红色表示正在燃烧
C(2:301,2:301,2)=G;%绿色表示绿树
Ci=imshow(C);%将初始状态绘制出来
ti=0;%时间
tp=title(['T = ',num2str(ti)]);
while 1
    ti=ti+1;
    St=Sk;
    Sf=Sk;
    Sf(2:301,2:301)=Sf(1:300,2:301)+Sf(2:301,1:300)+Sf(2:301,3:302)+Sf(3:302,2:301);
    St(Sf>4)=2;
    Sk=St;
    S=Sk(2:end-1,2:end-1);
    if (sum(sum(S)))==2*300*300
        break;
    end
    R=zeros(302);
    G=zeros(302);
    R(Sk==2)=1;
    G(Sk==1)=1;
    C(:,:,1)=R;
    C(:,:,2)=G;
    set(Ci,'CData',C);
    set(tp,'string',['T = ',num2str(ti)])
    pause(0.01);
end
 
