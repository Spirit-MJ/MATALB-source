function cengcifenxi
close
clear
clc
A=[1 1 1 4 1;1 1 2 4 1; 1 1/2 1 5 3;1/4 1/4 1/5 1 1/3;1 1 1/3 3 1];
B1=[1 1/4 1/2;4 1 3;2 1/3 1];
B2=[1 1/4 1/4;4 1 1/2;4 2 1];
B3=[1 3 4;1/3 1 1;1/4 1 1];
B4=[1 1/3 5;3 1 7;1/5 1/7 1];
B5=[1 1 5;1 1 3;1/5 1/3 1];
%A=input('请输入成对比较矩阵：回车键退出\n');
%B1=input('请输入运输时间矩阵：回车键退出\n');
%B2=input('请输入货物特性矩阵：回车键退出\n');
%B3=input('请输入运输成本矩阵：回车键退出\n');
%B4=input('请输入运输效率矩阵：回车键退出\n');
%B5=input('请输入运输安全矩阵：回车键退出\n');
[x1,y1,CR1]=tezheng(A);
[x2,y2,CR2]=tezheng(B1);
[x3,y3,CR3]=tezheng(B2);
[x4,y4,CR4]=tezheng(B3);
[x5,y5,CR5]=tezheng(B4);
[x6,y6,CR6]=tezheng(B5);
disp(['对应矩阵A、B1、B2、B3、B4、B5的特征值为：',num2str(y1),num2str(y2),num2str(y3),num2str(y4),num2str(y5),num2str(y6)]);
disp(['对应矩阵A的特征向量为：',num2str(x1)]);
disp(['对应矩阵B1的特征向量为：',num2str(x2)]);
disp(['对应矩阵B2的特征向量为：',num2str(x3)]);
disp(['对应矩阵B3的特征向量为：',num2str(x4)]);
disp(['对应矩阵B4的特征向量为：',num2str(x5)]);
disp(['对应矩阵B5的特征向量为：',num2str(x6)]);
W3=[x2',x3',x4',x5',x6'];
W=W3*x1';
disp(W);
end
function  [x,y,CR]=tezheng(a)
[x0,y0]=eig(a);
d=diag(y0);
y=max(d);
i=find(d==y);
w=x0(:,i);
x=w';
CI=(y-length(a))/(length(a)-1);
if CI>0
i=length(a);
switch i
    case 3
        RI=0.58;
    case 4
        RI=0.9;
    case 5
        RI=1.12;
    case 6
        RI=1.24;
    case 7
        RI=1.32;
    case 8
        RI=1.41;
    case 9
        RI=1.45;
    otherwise
        RI=input('请输入RI:');
end
CR=CI/RI;
if CR<0.1
    disp('通过一致性检验');
else
    disp('未通过一致性检验');
end
else
    disp('通过一致性检验');
end
end