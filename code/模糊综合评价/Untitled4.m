close all
clear all
clc
W=[0.5 0.2 0.3];%确定各评价因素的权重
R1=[160/600 380/600 20/600  40/600
    180/600 250/600 130/600 40/600
    130/600 270/600 130/600 70/600]%确定评价矩阵
R2=[170/650 410/650 10/650  60/650
    200/650 310/650 120/650 20/650
    110/650 320/650 120/650 100/650]
S1=W*R1%主因素突出型算子
S2=W*R2
Dengji=[4 3 2 1];%将评价等级赋值
u1=S1*Dengji'
u2=S2*Dengji'




