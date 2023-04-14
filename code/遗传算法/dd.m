close all
clear
clc

or_data_alt = xlsread('D:\΢��\WeChat Files\wxid_wtv61giv94jk22\FileStorage\File\2023-02\ģ������\��ѡ������.xlsx');
or_data_need = xlsread('D:\΢��\WeChat Files\wxid_wtv61giv94jk22\FileStorage\File\2023-02\ģ������\�����.xlsx');
R = 6371; % ����뾶6371km

xyz_alt = wesntoxyz(or_data_alt(:,2), or_data_alt(:, 3), R);
xyz_need = wesntoxyz(or_data_need(:,2), or_data_need(:, 3), R);
d = distanse(xyz_need, xyz_alt);
D = 10; % ����뾶
p = 40; % ���ñ�ѡ���������
cover_threshold = 0.8; % ����㸲����
w = zeros(length(or_data_need), 1); % �����Ȩ��
w(or_data_need(:, 4) == 1) = 0.3;
w(or_data_need(:, 4) == 2) = 0.6;
w(or_data_need(:, 4) == 3) = 0.9;

% �㷨������
np = 50; % ��Ⱥ��С
pc = 0.95; % �������
pm = 0.1; % �������
ng = 200; % ��������
GGAP = 0.9; % ��Ⱥ����
best_fx = zeros(ng, 3);
best_x = zeros(ng, length(or_data_alt));

x = zeros(np, length(or_data_alt));
for i = 1:np
    x(i, :) = round(rand(1, length(or_data_alt)));
    flag = judge_const(x(i, :), p, d, D, cover_threshold);
    while flag == 1
        x(i, :) = round(rand(1, length(or_data_alt)));
        flag = judge_const(x(i, :), p, d, D, cover_threshold);
    end
end

for k = 1:ng
    [fx, fx1, fx2] = fitness(x, d, D, w);
    nx = select(x, GGAP, fx);
    best_fx(k, 1) = max(fx);
    index1 = find(fx == max(fx));
    best_fx(k, 2) = fx1(index1(1));
    best_fx(k, 3) = fx2(index1(1));
    best_x(k, :) = nx(1, :); 
    nx = cross(nx, pc, p, d, D, cover_threshold);
    nx = variation(nx, pm, p, d, D, cover_threshold);
    x = reins(x, nx, fx);
    if k>=2
        figure(1)
        line([k-1,k],[best_fx(k-1, 1),best_fx(k, 1)]);
        xlabel('�Ŵ�����');
        ylabel('��Ӧ��ֵ');
        title('��������');
        grid on
    end
end

index_bext_x = find(best_x(end, :) == 1);
disp('ѡ�еı�ѡ��Ϊ��')
disp(index_bext_x)
best_pop = d(:, best_x(end, :)==1); 
res = sum(best_pop<=D, 2);
index_bext_y = find(res>0);
disp('���ǵ������Ϊ��')
disp(index_bext_y')

function xyz = wesntoxyz(lng, lat, R) % ��γ��ת����ά����
    xyz = zeros(length(lng), 3);
    xyz(:, 1) = R * cos(lat /180 * pi) .* cos(lng /180 * pi);
    xyz(:, 2) = R * cos(lat /180 * pi) .* sin(lng /180 * pi);
    xyz(:, 3) = R * sin(lat /180 * pi);
end

function D = distanse(need, alt) % ����������
    row = length(need);
    col = length(alt);
    D = zeros(row, col);
    for i = 1:row
        for j = 1:col
            D(i,j) = ((need(i,1) - alt(j,1))^2 + ...
                      (need(i,2) - alt(j,2))^2 + ...
                      (need(i,3) - alt(j,3))^2)^0.5;
        end
    end
end

function flag = judge_const(pop, p, d, D, threshold)
    flag = 0;
    if sum(pop) > p % Լ��1
        flag = 1;
    end
    
    pop1 = d(:, pop==1);  % Լ��2
    res = sum(pop1<=D, 2);
    res(res>0) = 1;
    if sum(res) <= (size(d, 1) * threshold)
        flag = 1;
    end
end

function [fx, Max, Min] = fitness(Chrom, d, D, w) % ������Ӧ��
    fx = zeros(size(Chrom, 1), 1);
    Max = zeros(size(Chrom, 1), 1);
    Min = zeros(size(Chrom, 1), 1);
    for i = 1:size(Chrom, 1)
        pop = d(:, Chrom(i, :)==1); 
        res = sum(pop<=D, 2);
        res(res>0) = 1;
        Max(i) = sum(res .* w);
        Min(i) = sum(Chrom(i, :));
        fx(i, 1) =  sum(res .* w) / (sum(Chrom(i, :)) + 1e-5);
    end
end

function nx = select(Chrom, GGAP, fx) % ѡ������
       GGAP1 = floor(GGAP * size(Chrom, 1));
       nx = zeros(GGAP1, size(Chrom, 2));
       [~, temp] = sort(fx, 'descend'); 
       for j = 1:GGAP1
           nx(j, :) = Chrom(temp(j), :); 
       end
end

function Chrom = cross(Chrom, pc, p, d, D, threshold) %��������
    i = 1;
    while i <= size(Chrom, 1)
        if pc>=rand(1)
            rand_loc = randperm(size(Chrom, 1), 2);
            rand_num = randperm(size(Chrom, 2), 2);
            temp = Chrom(rand_loc(1), min(rand_num):max(rand_num));
            Chrom(rand_loc(1), min(rand_num):max(rand_num)) = ...
            Chrom(rand_loc(2), min(rand_num):max(rand_num));
            Chrom(rand_loc(2), min(rand_num):max(rand_num)) = temp;
           if (judge_const(Chrom(rand_loc(1), :), p, d, D, threshold) == 1) || ...
              (judge_const(Chrom(rand_loc(2), :), p, d, D, threshold) == 1)
              temp = Chrom(rand_loc(1), min(rand_num):max(rand_num));
              Chrom(rand_loc(1), min(rand_num):max(rand_num)) = ...
              Chrom(rand_loc(2), min(rand_num):max(rand_num));
              Chrom(rand_loc(2), min(rand_num):max(rand_num)) = temp; 
              i = i - 1;
           end
        end
        i = i + 1;
    end
end

function Chrom = variation(Chrom, pm, p, d, D, threshold) % �������� 
    i = 1;
    while i <= size(Chrom, 1)
        if pm >= rand(1)
            rand_loc = randperm(size(Chrom, 1), 1);
            rand_num = randperm(size(Chrom, 2), 1);
            Chrom(rand_loc, rand_num) = ~Chrom(rand_loc, rand_num);
            while (judge_const(Chrom(rand_loc, :), p, d, D, threshold) == 1) || ...
                  (judge_const(Chrom(rand_loc, :), p, d, D, threshold) == 1)
              Chrom(rand_loc, rand_num) = ~Chrom(rand_loc, rand_num);
              i = i - 1;
            end
        end
        i = i + 1;
    end
end

function Chrom = reins(Chrom, nx, fx) % ��ȫ����
    init_len = size(Chrom, 1);
    new_len = size(nx, 1);
    [~, index] = sort(fx, 'descend');
    Chrom = [Chrom(index(1:init_len-new_len), :); nx]; 
end

