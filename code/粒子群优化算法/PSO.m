close all
clear
clc
tic
x=-1.5:0.01:1.5;
y=-2:0.01:2;
for i=1:length(y)
    for j=1:length(x)
        z(i,j)=sin(sqrt((y(i))^2+x((j))^2))/sqrt((y(i))^2+(x(j))^2)+exp((cos(2*pi*x(j))+cos(2*pi*y(i)))/2)-2.71289;
    end
end
figure(1)
mesh(x,y,z);
grid on
%��ʼ������
c1=1.49445;
c2=1.49445;
maxgen=100;   % ��������  
sizepop=100;   %��Ⱥ��ģ
Vmax=0.5;  %�����Ա���ȡֵ  
Vmin=-0.5;
popmax=2;%����ֵ��Χ
popmin=-2;
%������ʼ���Ӻ��ٶ�
for i=1:sizepop
    % �������һ����Ⱥ
    pop(i,:)=2*rands(1,2);  %��ʼ��Ⱥ���������-2����2����
    V(i,:)=0.5*rands(1,2);  %��ʼ���ٶȣ��������-0.5��0.5֮�����
    %������Ӧ��
    fitx(i)=fitness(pop(i,:));  %���㺯��ֵ 
end
%���弫ֵ��Ⱥ�弫ֵ
[bestfitness,bestindex]=max(fitx);%�ҳ���Ӧ������
Pbest=pop;   %��������
Gbest=pop(bestindex,:);   %ȫ����ѣ��ҵ�����x
fitnessPbest=fitx;   %���������Ӧ��ֵ���˴�������Ӧ��ֵ(Yֵ)
fitnessGbest=bestfitness;   %ȫ�������Ӧ��ֵ���˴����ŵ���Ӧ��ֵ
% VI. ����Ѱ��
for i=1:maxgen
    for j=1:sizepop
        % �ٶȸ���
        V(j,:)=V(j,:)+c1*rand(1)*(Pbest(j,:)-pop(j,:))+c2*rand(1)*(Gbest-pop(j,:));
        V(j,find(V(j,:)>Vmax))=Vmax;
        V(j,find(V(j,:)<Vmin))=Vmin;
        
        % ��Ⱥ����
        pop(j,:)=pop(j,:)+V(j,:);
        pop(j,find(pop(j,:)>popmax))=popmax;
        pop(j,find(pop(j,:)<popmin))=popmin;
        % ��Ӧ��ֵ����
        fitx(j)=fitness(pop(j,:)); 
    end
    for j=1:sizepop    
        % �������Ÿ���
        if fitx(j)>fitnessPbest(j) %������º����Ⱥ����ĳ��������Ӧ�ȴ���ǰһ�ζ�Ӧ����Ӧ��ֵ
            Pbest(j,:)=pop(j,:); %�ͽ��������(x��ֵ)��ֵ��gbest(�˴���Ӧ��xֵ)�������� 
            fitnessPbest(j)=fitx(j);%���´˴���Ӧ����Ӧ��ֵ(Yֵ)
        end
        % Ⱥ�����Ÿ���
        if fitx(j)>fitnessGbest %������º����Ⱥ����ĳ��������Ӧ�ȴ���ǰһ�����ŵ���Ӧ��ֵ
            Gbest=pop(j,:);%�����ŵ��������(x��ֵ)��¼������ֵ��zbest
            fitnessGbest=fitx(j);%�����ŵ���Ӧ��ֵ(y��ֵ)��¼������ֵ��fitnesszbest
        end
        fx(i)=fitnessGbest;
    end  
    figure(2)
    if i>1
        line([i-1,i],[fx(i-1),fx(i)]);
    end
    title('���Ÿ�����Ӧ��','fontsize',12);
    xlabel('��������','fontsize',12);
    ylabel('��Ӧ��','fontsize',12);
end
% VII. ����������ͼ
disp(['���������ֵΪ��',num2str(fx(end))])
figure(3)
hold on
mesh(x,y,z);
plot3(Gbest(1),Gbest(2),fx(end),'r*')
toc
function y=fitness(x)
y=sin(sqrt((x(1))^2+(x(2))^2))/sqrt((x(1))^2+(x(2))^2)+exp((cos(2*pi*x(1))+cos(2*pi*x(2)))/2)-2.71289;
end


