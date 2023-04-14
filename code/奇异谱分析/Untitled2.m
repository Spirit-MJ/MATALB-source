close all
clear
clc
t = linspace(0,2*pi*10,1000);
x = sin(t);
N = length(t);
x_ba = mean(x);


rx(x,N,x_ba)

function rx(x,N,x_ba)
for T = 1:100
    temp0 = 0;
    temp1 = 0;
    for i = 1:N-T
        temp0 = temp0 + (x(i)-x_ba)*(x(i+T)-x_ba);
        temp1 = temp1 + x(i)*x(i+T);
    end
    fenmu = 1/N * sum((x - x_ba).^2);
    Rxx(T) = 1/N*temp0/fenmu;  %% 选择第一次为零的值
    Rxx0(T) = 1/N*temp1;  %% 选择下降到初始值的1-e-1的的值
end

tag0 = (1-1/exp(1))*Rxx0(1);

x1 = x(1:N-15);
x2 = x(16:N);
x3 = x(1:N-25);
x4 = x(26:N);

figure(1)
plot(1:length(Rxx),Rxx)
figure(2)
subplot(121)
plot(x2,x1,'r-o')
subplot(122)
plot(x4,x3,'r-o')
end

