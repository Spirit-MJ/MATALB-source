close all
clear
clc
syms x1 x2 ; 
f=(x1-2)^2+2*(x2-1)^2;
figure(1)
hold on
fmesh((x1-2)^2+2*(x2-1)^2,[-100,100,-100,100])
x=[100;100];
plot3(x(1),x(2),(x(1)-2)^2+2*(x(2)-1)^2,'k*')
e=10^(-20);
d=-[diff(f,x1);diff(f,x2)];  %�ֱ���x1��x2��ƫ���������½��ķ���
epoch=50;  % ѧϰ����
lr = 0.1;   % ѧϰ��
loss = zeros(epoch);
for i = 1:epoch
    
    % �����ݶ�ֵ
    d_temp=subs(d,x1,x(1));     
    d_temp=subs(d_temp,x2,x(2)); 
    
    %�����ݶ�
    x=x+lr*d_temp; 
    
    % ����loss
    f_temp = subs(f,x1,x(1)); 
    f_temp = subs(f_temp,x2,x(2)); 
    loss(i) = f_temp;

    if i>=2
        figure(1)
        hold on
        plot3(x(1),x(2),(x(1)-2)^2+2*(x(2)-1)^2,'k*');
        figure(2)
        line([i-1,i],[loss(i-1),loss(i)]);
    end
end
ender=double(x)  %�յ�


