close all
clear
clc

figure(1)
plot([1,200],[1,200],'.'); % �滮�ռ�
hold on
xlabel('km','fontsize',12)
ylabel('km','fontsize',12)
title('��ά�滮�ռ�','fontsize',12)
% ���������յ�
Start = [20,180];
End = [160,90];
plot([Start(1),End(1)],[Start(2),End(2)],'.');

% ͼ�α�ע
text(Start(1)+2,Start(2),'���');  % ���
text(End(1)+2,End(2),'�յ�');  % �յ�

position = [40,120;40,170;80,170;80,120;
            100,120;100,150;140,150;140,120;
            150,100;150,120;170,120;170,100;
            80,100;80,110;120,110;120,100];  % �ϰ���
        
% ����ϰ���ͼ��
fill(position(1:4,1),position(1:4,2),[0,0,0]);
fill(position(5:8,1),position(5:8,2),[0,0,0]);
fill(position(9:12,1),position(9:12,2),[0,0,0]);
fill(position(13:16,1),position(13:16,2),[0,0,0]);

%�㷨����
ant_num = 100;  % ��������
rou = 0.2;
view_h = 10;  % ����������С
view_w = abs(End(2) - Start(2)) + 1;  % ����������С 
num_gen = floor((End(1)-Start(1)-1)/view_h);  % ����ʵĵ�x������
T_fiil0 = ones(200);  % ��ʼ����Ϣ��
for ii = 1:length(position)/4
    for i = position(ii*4-3,1):position(ii*4-1,1)
        for j = position(ii*4-3,2):position(ii*4-1,2)
            T_fiil0(i,j) = 0;
        end
    end
end  
T_fill = zeros(num_gen, view_w);
for i0 = 1:num_gen
    T_fill(i0,:) = T_fiil0(Start(1)+i0*view_h,90:180);
end
T = ones(num_gen, view_w);  % ��ʼ����Ϣ��
mingen = 2;  % ����������ʼֵ
maxgen = 500;  % �����������ֵ
zuiyou_chang = zeros(1,maxgen);  % ��¼ÿ�ε��������ž���
zuiyou_route = zeros(maxgen,num_gen+2);  % ��¼ÿ�ε���������·��

% ��ʼ����·��
[path,T] = find_path(ant_num, num_gen, Start, End, view_h, T, T_fill);  % ��ʼ����·��
fitness = comput_fit(path,view_h);  %��Ӧ�ȼ���
[min_chang,min_index] = min(fitness);           %�����Ӧ��
zuiyou_route(1,:) = path(min_index,:);     %���·��
zuiyou_chang(1) = min_chang;               %��Ӧ��ֵ��¼

% ��Ϣ�ظ���
cfit = 1 / zuiyou_chang(1);
mm = 1;
for i = 2:length(zuiyou_route(1,:))-1
    T(mm,zuiyou_route(1,i)-90+1) = (1-rou)*T(mm,zuiyou_route(1,i)-90+1)+ rou * cfit;
    mm = mm + 1;
end

while mingen <= maxgen
    % ·������
    [path,T] = find_path(ant_num, num_gen, Start, End, view_h, T, T_fill);  % ��ʼ����·��
    
    % ��Ӧ��ֵ�������
    fitness = comput_fit(path,view_h);  %��Ӧ�ȼ���
    [min_chang,min_index] = min(fitness);           %�����Ӧ��
    
    if min_chang < zuiyou_chang(mingen-1)  % �жϴ˴ξ����Ƿ���ϴξ����
         zuiyou_chang(mingen) = min_chang;  % �ǵĻ���¼�˴���̾���
         zuiyou_route(mingen,:) = path(min_index,:);  % �ǵĻ���¼�˴�����·��
    else
         zuiyou_chang(mingen) = zuiyou_chang(mingen-1);  % ����Ļ���¼��һ����̾���
         zuiyou_route(mingen,:) = zuiyou_route(mingen-1,:);  % ����Ļ���¼��һ������·��
    end
    
    % ������Ϣ��
    cfit = 1 / zuiyou_chang(mingen);
    mm = 1;
    for i = 2:length(zuiyou_route(mingen,:))-1
        T(mm,zuiyou_route(mingen,i)-90+1) = (1-rou)*T(mm,zuiyou_route(mingen,i)-90+1)+ rou * cfit;
        mm = mm + 1;
    end

    figure(2)  % ��������ͼ��
    hold on
    line([mingen-1,mingen],[zuiyou_chang(mingen-1),zuiyou_chang(mingen)]);
    xlabel('��������')
    ylabel('����')
    title('��������')
    hold off
    mingen = mingen+ 1;
end
loc_x = Start(1);
for i = 1:length(zuiyou_route(end,:))-2
    tenp = Start(1) + view_h * i;
    loc_x = [loc_x,tenp];
end
loc_x = [loc_x,End(1)];
figure(1)
plot(loc_x,zuiyou_route(end,:),'r-.'); % �滮�ռ�
disp(['���ų�Ϊ��',num2str(zuiyou_chang(end))])

function [path, T] = find_path(ant_num, num_gen, Start, End, view_h, T, T_fill)
path = zeros(ant_num, num_gen+2);
for ii = 1:ant_num

    path(ii,1) = Start(2);  % ��¼��ʼ����
    temp = 1;
    now_point = Start(2);      % ��ǰ����
    
    % �������Ӧ��ֵ
    
    for loc_x = 1:num_gen  % ���������
        temp = temp + 1;
        %�����������ݵ��Ӧ����Ӧ��ֵ
        kk = 1;
        for i = 90:180  % ��������
            next_point(kk) = i;  % �����������·��
            qfz(kk) = comput_qfz(Start(1)+(loc_x-1)*view_h, now_point, next_point(kk),End(1),End(2),view_h);  % ��������ֵ
            qz(kk) = qfz(kk) * T(loc_x,next_point(kk)-90+1) * T_fill(loc_x,next_point(kk)-90+1);  % ����������е�ĸ���
            kk = kk+1;
        end
        %���̶�ѡ����ʵ���һ���ص�
        p = qz/sum(qz);
        pc = cumsum(p); 
        target_index = find(pc>=rand(1));  
        index = target_index(1);  % �ҳ�������һ�������ȥ�ĵ�
        
        old_point = next_point(index);

        %������Ϣ��
        T(loc_x,old_point-90+1) = 0.9 * T(loc_x,old_point-90+1);  % ��Ϣ�ػӷ�
        
        %·������
        path(ii,temp) = old_point;
        now_point = old_point; 
    end
    path(ii,temp+1) = End(2);  % �յ�
end
end

function qfz = comput_qfz(now_x,now_y,next_y,end_x,end_y,view_h)
next_x = now_x + view_h;

%D����
D = 1 / (sqrt(view_h^2 + (now_y-next_y)^2) + ...   % ��ǰ������һ�����
    sqrt((next_x-end_x)^2+(next_y-end_y)^2));  % ��һ���е㵽�յ�ľ��� 
%��������ֵ
qfz = D;
end

function fitness = comput_fit(path,view_h)
[n,m] = size(path); 
fitness = zeros(n,1);
for i = 1:n
    fitness(i) = 0;
    for j = 2:m
        fitness(i) = fitness(i)+ sqrt(view_h^2 + (path(i,j)-path(i,j-1))^2);  
    end
end
end
% mmm = 0;
% for i = 1:length(zuiyou_route(end,:))-1
%     mmm = mmm + sqrt(view_h^2+(zuiyou_route(end,i) + zuiyou_route(end,i+1))^2);
% end
