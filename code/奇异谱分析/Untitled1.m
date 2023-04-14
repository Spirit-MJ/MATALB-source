close all
clear
clc

x = linspace(0,2*pi*10,1000);
noise_true = normrnd(0,0.3,[1,length(x)]);
y_true = sin(x);
y = sin(x) + noise_true;


N = length(y);
windowLen = 100;  % 窗口长度
K = N - windowLen + 1;
X = zeros(K, windowLen);

%  建立轨迹矩阵
for i=1:K
    X(i,1:windowLen) = y(i:windowLen+i-1);
end

%  SVD分解
[U,S,V] = svd(X);
lamda = diag(S);
per_lamda = lamda/sum(lamda); % 累计贡献率

%  分组
V = V';
r = 2; % 前r个奇异值
New_X0 = U(:, 1:r)*S(1:r, 1:r)*V(1:r, :);  % 信号
New_X1 = U(:, r+1:windowLen)*S(r+1:windowLen, r+1:windowLen)*V(r+1:windowLen, :); % 噪声

% 重构
y_new0 = Recory(N, New_X0);  % 信号
y_new1 = Recory(N, New_X1);  % 噪声

figure(1)
subplot(121)
box on
plot(x, y_true, 'b', 'LineWidth', 1)
axis([0, x(end),-2, 2])
legend('真实数据')
subplot(122)
box on
plot(x, y, 'r', 'LineWidth',1)
axis([0, x(end),-2, 2])
legend('带噪声的数据')

figure(2)
subplot(121)
box on
plot(x, y_true, 'b', 'LineWidth',1)
axis([0, x(end),-2, 2])
legend('真实数据')
subplot(122)
box on
plot(x, y_new0, 'r', 'LineWidth',1)
axis([0, x(end),-2, 2])
legend('SVD分解')

figure(3)
subplot(121)
box on
plot(x, y_new1, 'b', 'LineWidth',1)
axis([0, x(end),-2, 2])
legend('真实噪声')
subplot(122)
box on
plot(x, noise_true, 'r', 'LineWidth',1)
legend('SVD分解')
axis([0, x(end),-2, 2])

figure(4)
pareto(per_lamda)  %主成分贡献率
xlabel('奇异值')
ylabel('贡献率(%)')
title('贡献率')


function N_X = Recory(n, New_X)
    y_temp = zeros(1,n);
    y_num = zeros(1,n);
    for i = 1:size(New_X,1)
        for j  = 1:size(New_X,2)
            y_temp(j + i - 1) = y_temp(j + i - 1) + New_X(i,j);
            y_num(j + i - 1) = y_num(j + i - 1) + 1;
        end
    end
    N_X = y_temp./y_num;
end

