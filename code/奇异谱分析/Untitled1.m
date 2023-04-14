close all
clear
clc

x = linspace(0,2*pi*10,1000);
noise_true = normrnd(0,0.3,[1,length(x)]);
y_true = sin(x);
y = sin(x) + noise_true;


N = length(y);
windowLen = 100;  % ���ڳ���
K = N - windowLen + 1;
X = zeros(K, windowLen);

%  �����켣����
for i=1:K
    X(i,1:windowLen) = y(i:windowLen+i-1);
end

%  SVD�ֽ�
[U,S,V] = svd(X);
lamda = diag(S);
per_lamda = lamda/sum(lamda); % �ۼƹ�����

%  ����
V = V';
r = 2; % ǰr������ֵ
New_X0 = U(:, 1:r)*S(1:r, 1:r)*V(1:r, :);  % �ź�
New_X1 = U(:, r+1:windowLen)*S(r+1:windowLen, r+1:windowLen)*V(r+1:windowLen, :); % ����

% �ع�
y_new0 = Recory(N, New_X0);  % �ź�
y_new1 = Recory(N, New_X1);  % ����

figure(1)
subplot(121)
box on
plot(x, y_true, 'b', 'LineWidth', 1)
axis([0, x(end),-2, 2])
legend('��ʵ����')
subplot(122)
box on
plot(x, y, 'r', 'LineWidth',1)
axis([0, x(end),-2, 2])
legend('������������')

figure(2)
subplot(121)
box on
plot(x, y_true, 'b', 'LineWidth',1)
axis([0, x(end),-2, 2])
legend('��ʵ����')
subplot(122)
box on
plot(x, y_new0, 'r', 'LineWidth',1)
axis([0, x(end),-2, 2])
legend('SVD�ֽ�')

figure(3)
subplot(121)
box on
plot(x, y_new1, 'b', 'LineWidth',1)
axis([0, x(end),-2, 2])
legend('��ʵ����')
subplot(122)
box on
plot(x, noise_true, 'r', 'LineWidth',1)
legend('SVD�ֽ�')
axis([0, x(end),-2, 2])

figure(4)
pareto(per_lamda)  %���ɷֹ�����
xlabel('����ֵ')
ylabel('������(%)')
title('������')


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

