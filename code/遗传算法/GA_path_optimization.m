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

np = 1000;  % ��Ⱥ��С
pc = 0.9;  % �������
pm = 0.2;  % �������
ng = 200;  % ��������
GGAP = 0.95;  % ��Ⱥ����
view = 9;  % �������С
num_gen = floor((End(1)-Start(1)-1)/view);  % ÿ��Ⱦɫ���ϻ������
zuiyou_fit = zeros(ng,1);
zuiyou_path = zeros(ng,num_gen); 

x = zeros(np,num_gen);
for i=1:np
    x(i,:) = ceil(90+90*rand(1,num_gen));  % ��ʼ����Ⱥʵ������
end

for k = 1:ng
    fx = fitness(x,Start,End,view);  % ������Ӧ��
    [zuiyou,index] = max(fx);
    zuiyou_fit(k) = 1/zuiyou;
    zuiyou_path(k,:) = x(index,:);
    nx = XuanZe(x,GGAP,np,fx);  % ѡ�����
    nx = JiaoCha(nx,pc);  % �������
    nx = BianYi(nx,pm);  % ����
    x = Reins(x,nx,fx);  % ��ȫ��Ⱥ
    figure(2)
    if k>=2
        line([k-1,k],[zuiyou_fit(k-1),zuiyou_fit(k)]);
        xlabel('�Ŵ�����');
        ylabel('����');
        title('��������');
        grid on
    end
end
figure(1)
hold on
loc_y = [Start(2),zuiyou_path(end,:),End(2)];
loc_x = Start(1);
for i = 1:length(zuiyou_path(end,:))
    tenp = Start(1) + view * i;
    loc_x = [loc_x,tenp];
end
loc_x = [loc_x,End(1)];
plot(loc_x,loc_y,'r-.'); % �滮�ռ� 
disp(['���ų�Ϊ',num2str(zuiyou_fit(end)),'km'])

function fit = fitness(x0,Start,End,view)  % ������Ӧ��
fit = zeros(size(x0,1),1);
for j = 1:size(x0,1)
    fit(j) = sqrt(view^2 + (Start(2) - x0(j,1))^2);
    for i = 1:size(x0,2)-1
        fit(j) = fit(j) + sqrt(view^2 + (x0(j,i) - x0(j,i+1))^2);
    end
    fit(j) = 1/(fit(j) + sqrt(view^2 + (End(2) - x0(j,end))^2));
end

end

function nx = XuanZe(x0,GGAP,np,fx)  % ѡ������
       GGAP1 = floor(GGAP*np);  %ѡ�����µĸ�����      
       [~,temp] = sort(fx,'descend');  % ������Ӧ�Ƚ�������
       for j = 1:GGAP1
           nx(j,:) = x0(temp(j),:);  % ��ʤ��̭������Ӧ�ȴ�ĸ��屣������
       end
end

function x0 = JiaoCha(x0,Pc)
 [NSel,L] = size(x0);  % ����ѡ���������Ⱥ��С
    for i = 1:NSel/2
        if Pc>=rand(1)
            rand_2 = randperm(L,2);  % ���ѡ�񽻲��λ��
            temp = x0(2*i-1,min(rand_2):max(rand_2));
            x0(2*i-1,min(rand_2):max(rand_2)) = x0(2*i,min(rand_2):max(rand_2));  % ����
            x0(2*i,min(rand_2):max(rand_2)) = temp;  % ����
        end
    end
end

function x0 = BianYi(x0,Pm)  % ����
    [NSel,L] = size(x0);  % �����������Ⱥ�Ĵ�С
    for i = 1:NSel
        if Pm>=rand
            R0 = randperm(L,1);  % ���ѡ�����Ļ���λ��
            R1 = ceil(90+90*rand(1));  % ���ѡ�����Ļ���
            x0(i,R0) = R1;  % ����
        end
    end
end
 
function Chrom = Reins(Chrom,SelCh,ObjV)  % ��ȫ��Ⱥ
    NIND = size(Chrom,1);  % ��ʼ��Ⱥ��С
    NSel = size(SelCh,1);  % ѡ���������Ⱥ�Ĵ�С
    [~,index] = sort(ObjV,'descend');  % ���Ŵ�����ǰ�ľ�������
    size(Chrom(index(1:NIND-NSel),:));
    Chrom = [SelCh;Chrom(index(1:NIND-NSel),:)];  % ���Ŵ�����ǰ���ŵ�һЩ���屣������
end

