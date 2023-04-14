close all
clear
clc
data_origin = xlsread('E:\����\��γ�ȸ߳�.xls');
origin_start = [120.767,52.183,419];  % ��ʼ��
origin_end = [121.413,52.360,675];  % �յ�

[X,Y,Z]=griddata(data_origin(:,2),data_origin(:,3),data_origin(:,1), ...
    linspace(min(data_origin(:,2)),max(data_origin(:,2)),100)', ...
    linspace(min(data_origin(:,3)),max(data_origin(:,3)),100));  % �������
figure(1)
mesh(X,Y,Z)
xlabel('����');
ylabel('γ��');
zlabel('����');
title('�������')

Z(ismissing(Z)) = max(max(Z));

%  ��ͼ��ģ
LevelGrid = 500;  % ���߶Ȼ��ֳ�500��
LevelGrid_h = (max(max(Z)) - min(min(Z))) / LevelGrid;  % �����ε����֮��������룩

PortGrid = 100;  % ��x��y���ֳ�100��
PortGrid_x = (max(data_origin(:,2)) - min(data_origin(:,2))) / PortGrid;  % �����ȵ����֮��������룩
PortGrid_y = (max(data_origin(:,3)) - min(data_origin(:,3))) / PortGrid;  % ��ά�ȵ����֮��������룩

Grid_Z = floor((Z - min(data_origin(:,1))) / LevelGrid_h);  % ��ģ����Ժ��Σ�
Grid_Z(Grid_Z == 501) = 500;
start_x = round((origin_start(1) - min(data_origin(:,2))) / PortGrid_x);  % ��ģ����ʼ��λ��
start_y = round((origin_start(2) - min(data_origin(:,3))) / PortGrid_y);  % ��ģ����ʼ��λ��
start_h = round((origin_start(3) - min(min(Z))) / LevelGrid_h);  % ��ģ����ʼ��λ��
end_x = round((origin_end(1) - min(data_origin(:,2))) / PortGrid_x);  % ��ģ���յ�λ��
end_y = round((origin_end(2) - min(data_origin(:,3))) / PortGrid_y);  % ��ģ���յ�λ��
end_h = round((origin_end(3) - min(min(Z))) / LevelGrid_h);  % ��ģ���յ�λ��

%�㷨����
ant_num = 50;  % ��������
rou = 0.2;
T = ones(PortGrid,PortGrid,LevelGrid);  % ��ʼ����Ϣ��
mingen = 2;  % ����������ʼֵ
maxgen = 300;  % �����������ֵ
zuiyou_chang = zeros(1,maxgen);  % ��¼ÿ�ε��������ž���
zuiyou_route = zeros(maxgen, 2 * (end_x-start_x + 1));  % ��¼ÿ�ε��������ž���

% ��ʼ����·��
[path,T] = find_path(ant_num,T,PortGrid,start_x,start_y,start_h,end_x,end_y,end_h, ...
    Grid_Z,Z,data_origin,LevelGrid_h,PortGrid_x,PortGrid_y);  % ��ʼ����·��
fitness = comput_fit(path,data_origin,PortGrid_x,PortGrid_y,LevelGrid_h,Z);  %��Ӧ�ȼ���
[min_chang,min_index] = min(fitness);           %�����Ӧ��
zuiyou_route(1,:) = path(min_index,:);     %���·��
zuiyou_chang(1) = min_chang;               %��Ӧ��ֵ��¼

% ��Ϣ�ظ���
cfit = 100 / zuiyou_chang(1);
for i = 2:59
    T(i,zuiyou_route(1,i*2-1),zuiyou_route(1,i*2))= ...
        (1-rou)*T(i,zuiyou_route(1,i*2-1),zuiyou_route(1,i*2))+ rou * cfit;
end


while mingen <= maxgen
    % ·������
    [path,T] = find_path(ant_num,T,PortGrid,start_x,start_y,start_h,end_x,end_y,end_h, ...
    Grid_Z,Z,data_origin,LevelGrid_h,PortGrid_x,PortGrid_y); 
    
    % ��Ӧ��ֵ�������
    fitness = comput_fit(path,data_origin,PortGrid_x,PortGrid_y,LevelGrid_h,Z);                               
    [min_chang,min_index] = min(fitness);
    
    if min_chang < zuiyou_chang(mingen-1)  % �жϴ˴ξ����Ƿ���ϴξ����
         zuiyou_chang(mingen) = min_chang;  % �ǵĻ���¼�˴���̾���
         zuiyou_route(mingen,:) = path(min_index,:);  % �ǵĻ���¼�˴�����·��
    else
         zuiyou_chang(mingen) = zuiyou_chang(mingen-1);  % ����Ļ���¼��һ����̾���
         zuiyou_route(mingen,:) = zuiyou_route(mingen-1,:);  % ����Ļ���¼��һ������·��
    end
    
    % ������Ϣ��
    cfit = 100 / zuiyou_chang(mingen);
    for i = 2:59
        T(i,zuiyou_route(mingen,i*2-1),zuiyou_route(mingen,i*2)) = (1-rou) * ...
            T(i,zuiyou_route(mingen,i*2-1),zuiyou_route(mingen,i*2)) + rou*cfit;
    end

    figure(2)  % ��������ͼ��
    hold on
    if mingen>=2
        line([mingen-1,mingen],[zuiyou_chang(mingen-1),zuiyou_chang(mingen)]);
    end
    xlabel('��������')
    ylabel('����/km')
    title('��������')
    hold off
    mingen = mingen+ 1;
end
x_index = start_x:end_x;
true_x = x_index * PortGrid_x + min(min(data_origin(:,2)));
true_y = zuiyou_route(end,1:2:end) * PortGrid_y + min(min(data_origin(:,3)));
disp('����·��Ϊ��')
for i = 1:length(true_x)
    true_h(i) = Z(zuiyou_route(end,2*i-1),x_index(i));
    disp([true_x(i),true_y(i),true_h(i)])
end
%  ����ת��
true_xyz_x = (6371000 + true_h) .* cos(true_y ./ 180 * pi) .* cos(true_x ./ 180 * pi);
true_xyz_y = (6371000 + true_h) .* cos(true_y ./ 180 * pi) .* sin(true_x ./ 180 * pi);
true_xyz_z = (6371000 + true_h) .* sin(true_y ./ 180 * pi);
% ����ŷ�Ͼ���
distance = 0;
for i = 1:length(true_xyz_x)-1
    distance = distance + sqrt((true_xyz_x(i)-true_xyz_x(i+1))^2+(true_xyz_y(i)-true_xyz_y(i+1))^2+(true_xyz_z(i)-true_xyz_z(i+1))^2);
end

disp(['�ܾ���Ϊ��',num2str(distance/1000),'ǧ��'])
disp(['��ʱ��Ϊ��',num2str(distance/5/3600),'Сʱ'])

figure(1)
hold on
plot3(true_x,true_y,true_h,'r-*','linewidth',2)
plot3(true_x(1),true_y(1),true_h(1),'--o','LineWidth',2,...
                       'MarkerEdgeColor','k',...
                       'MarkerFaceColor','g',...
                       'MarkerSize',10)
plot3(true_x(end),true_y(end),true_h(end),'--o','LineWidth',2,...
                       'MarkerEdgeColor','k',...
                       'MarkerFaceColor','g',...
                       'MarkerSize',10)
text(true_x(1),true_y(1),true_h(1),'���');
text(true_x(end),true_y(end),true_h(end),'�յ�');


function [path, T] = find_path(ant_num,T,PortGrid,start_x,start_y,start_h,end_x,end_y,end_h,Grid_Z,Z,data_origin,LevelGrid_h,PortGrid_x,PortGrid_y)
% ��������
y_max = 8;   % ����γ�ȷ�����ӷ�Χ
path = zeros(ant_num, 2 * (end_x-start_x));
for ii = 1:ant_num

    path(ii,1:2) = [start_y, start_h];  % ��¼��ʼ����
    temp = 1;
    now_point=[start_y, start_h];      % ��ǰ����
    
    % �������Ӧ��ֵ
    for loc_x = start_x+1:end_x-1  % ������Ϊ���ȷ���
        temp = temp + 1;
        %�����������ݵ��Ӧ����Ӧ��ֵ
        kk = 1;
        for i = -y_max:y_max
            if (now_point(1)+i<=PortGrid) && (now_point(1)+i>0)  % ��֤��������������
                to_j = Grid_Z(now_point(1)+i,loc_x);
                next_point(kk,:) = [now_point(1)+i,to_j];  % ���ž��ȷ������·��
                qfz(kk) = comput_qfz(next_point(kk,1),next_point(kk,2),now_point(1),now_point(2),end_y,end_h,loc_x,end_x,Z,data_origin,LevelGrid_h,PortGrid_x,PortGrid_y);  % ��������ֵ
                qz(kk) = qfz(kk) * T(loc_x,next_point(kk,1),next_point(kk,2));  % ����������е�ĸ���
                kk = kk+1;
            else
                qz(kk) = 0;  % �����������������
                kk = kk + 1;
            end
        end
        
        %���̶�ѡ����ʵ���һ���ص�
        p = qz/sum(qz);
        pc = cumsum(p); 
        target_index = find(pc>=rand(1));  
        index = target_index(1);  % �ҳ�������һ�������ȥ�ĵ�
        
        old_point = next_point(index,:);

        %������Ϣ��
        T(loc_x+1,old_point(1),old_point(2)) = 0.9 * T(loc_x+1,old_point(1),old_point(2));
        
        %·������
        path(ii,temp*2-1:temp*2) = [old_point(1),old_point(2)];
        now_point = old_point; 
    end
    path(ii,(temp+1)*2-1:(temp+1)*2) = [end_y,end_h];  % �յ�
end
end

function qfz = comput_qfz(next_y,next_h,now_y,now_h,end_y,end_h,loc_x,end_x,Z,data_origin,LevelGrid_h,PortGrid_x,PortGrid_y)

% ��γ���뺣����ʵֵ
true_now_x = (loc_x-1) * PortGrid_x + min(data_origin(:,2));
true_now_y = now_y * PortGrid_y + min(data_origin(:,3));
true_now_h = now_h * LevelGrid_h + min(min(Z));

true_next_x = loc_x * PortGrid_x + min(data_origin(:,2));
true_next_y = next_y * PortGrid_y + min(data_origin(:,3));
true_next_h = next_h * LevelGrid_h + min(min(Z));


true_end_x = end_x * PortGrid_x + min(data_origin(:,2));
true_end_y = end_y * PortGrid_y + min(data_origin(:,3));
true_end_h = end_h * LevelGrid_h + min(min(Z));


% ת��ֱ������
true_now_xyz_x = (6371000 + true_now_h) * cos(true_now_y / 180 * pi) * cos(true_now_x / 180 * pi);
true_now_xyz_y = (6371000 + true_now_h) * cos(true_now_y / 180 * pi) * sin(true_now_x / 180 * pi);
true_now_xyz_z = (6371000 + true_now_h) * sin(true_now_y / 180 * pi);

true_next_xyz_x = (6371000 + true_next_h) * cos(true_next_y / 180 * pi) * cos(true_next_x / 180 * pi);
true_next_xyz_y = (6371000 + true_next_h) * cos(true_next_y / 180 * pi) * sin(true_next_x / 180 * pi);
true_next_xyz_z = (6371000 + true_next_h) * sin(true_next_y / 180 * pi);

true_end_xyz_x = (6371000 + true_end_h) * cos(true_end_y / 180 * pi) * cos(true_end_x / 180 * pi);
true_end_xyz_y = (6371000 + true_end_h) * cos(true_end_y / 180 * pi) * sin(true_end_x / 180 * pi);
true_end_xyz_z = (6371000 + true_end_h) * sin(true_end_y / 180 * pi);

% �ж��¸����Ƿ�ɴ�

slope = abs(true_now_xyz_z-true_next_xyz_z) / sqrt((true_now_xyz_x-true_next_xyz_x)^2+(true_now_xyz_y-true_next_xyz_y)^2);
if (atan(slope) / pi * 180) <= 30  % �¶�Լ��
    S = 1;
else
    S = 0;
end

%D����
D = 10000 / (sqrt((true_next_xyz_x-true_now_xyz_x)^2 + (true_next_xyz_y-true_now_xyz_y)^2+(true_next_xyz_z-true_now_xyz_z)^2) + ...  
    sqrt((true_end_xyz_x-true_next_xyz_x)^2+(true_end_xyz_y-true_next_xyz_y)^2+(true_end_xyz_z-true_next_xyz_z)^2));  % ��һ���е㵽�յ�ľ��� 
%��������ֵ
qfz = D * S;
end

function fitness = comput_fit(path,data_origin,PortGrid_x,PortGrid_y,LevelGrid_h,Z)
[n,m] = size(path); 
fitness = zeros(n,1);
for i = 1:n
    fitness(i) = 0;
    for j = 2:m/2
        
        % ��γ���뺣����ʵֵ
        true_path_x0 =  PortGrid_x + min(data_origin(:,2));
        true_path_x1 =  2 * PortGrid_x + min(data_origin(:,2));
        
        true_path_y0 = path(i,j*2-1) * PortGrid_y + min(data_origin(:,3));
        true_path_y1 = path(i,(j-1)*2-1) * PortGrid_y + min(data_origin(:,3));
    
        true_path_h0 = path(i,j*2) * LevelGrid_h + min(min(Z));
        true_path_h1 = path(i,(j-1)*2) * LevelGrid_h + min(min(Z));
        
        % ת��ֱ������
        true_xyz_x0 = (6371000 + true_path_h0) * cos(true_path_y0 / 180 * pi) * cos(true_path_x0 / 180 * pi);
        true_xyz_y0 = (6371000 + true_path_h0) * cos(true_path_y0 / 180 * pi) * sin(true_path_x0 / 180 * pi);
        true_xyz_z0 = (6371000 + true_path_h0) * sin(true_path_y0 / 180 * pi);
        
        true_xyz_x1 = (6371000 + true_path_h1) * cos(true_path_y1 / 180 * pi) * cos(true_path_x1 / 180 * pi);
        true_xyz_y1 = (6371000 + true_path_h1) * cos(true_path_y1 / 180 * pi) * sin(true_path_x1 / 180 * pi);
        true_xyz_z1 = (6371000 + true_path_h1) * sin(true_path_y1 / 180 * pi);
        
        %��Ӧ��ֵΪ���ȼӸ߶�
        fitness(i) = fitness(i)+ sqrt((true_xyz_x1-true_xyz_x0)^2 + (true_xyz_y1-true_xyz_y0)^2 + (true_xyz_z1-true_xyz_z0)^2)/1000;  % + abs(path(i,j*2));
    end
end
end

