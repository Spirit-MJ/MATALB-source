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
%A=input('������ɶԱȽϾ��󣺻س����˳�\n');
%B1=input('����������ʱ����󣺻س����˳�\n');
%B2=input('������������Ծ��󣺻س����˳�\n');
%B3=input('����������ɱ����󣺻س����˳�\n');
%B4=input('����������Ч�ʾ��󣺻س����˳�\n');
%B5=input('���������䰲ȫ���󣺻س����˳�\n');
[x1,y1,CR1]=tezheng(A);
[x2,y2,CR2]=tezheng(B1);
[x3,y3,CR3]=tezheng(B2);
[x4,y4,CR4]=tezheng(B3);
[x5,y5,CR5]=tezheng(B4);
[x6,y6,CR6]=tezheng(B5);
disp(['��Ӧ����A��B1��B2��B3��B4��B5������ֵΪ��',num2str(y1),num2str(y2),num2str(y3),num2str(y4),num2str(y5),num2str(y6)]);
disp(['��Ӧ����A����������Ϊ��',num2str(x1)]);
disp(['��Ӧ����B1����������Ϊ��',num2str(x2)]);
disp(['��Ӧ����B2����������Ϊ��',num2str(x3)]);
disp(['��Ӧ����B3����������Ϊ��',num2str(x4)]);
disp(['��Ӧ����B4����������Ϊ��',num2str(x5)]);
disp(['��Ӧ����B5����������Ϊ��',num2str(x6)]);
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
        RI=input('������RI:');
end
CR=CI/RI;
if CR<0.1
    disp('ͨ��һ���Լ���');
else
    disp('δͨ��һ���Լ���');
end
else
    disp('ͨ��һ���Լ���');
end
end