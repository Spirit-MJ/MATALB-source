%��ջ�������
close all
clear
clc
%��������
load iris_data.mat
%�������ѵ�����Ͳ��Լ�
train_P = [];
train_T = [];
test_P = [];
test_T = [];
for i = 1:3
    temp_input=features((i-1)*50+1:i*50,:);
    temp_output=classes((i-1)*50+1:i*50,:);
    n=randperm(50);
    %ѵ��������120������
    train_P=[train_P temp_input(n(1:40),:)'];
    train_T=[train_T temp_output(n(1:40),:)'];
    %���Լ�����30������
    test_P=[test_P temp_input(n(41:50),:)'];
    test_T=[test_T temp_output(n(41:50),:)'];
end
%ģ�ͽ��� 
train_p=train_P;
test_p=test_P;
%GRNN�������������
%��������
net_grnn=newgrnn(train_p,train_T);
%�������
t_sim_grnn=sim(net_grnn,test_p);
T_sim_grnn=round(t_sim_grnn);
result_grnn=T_sim_grnn';
%PNN�������������
Tc_train=ind2vec(train_T);
%��������
net_pnn=newpnn(train_p,Tc_train);
%�������
Tc_test=ind2vec(test_T);
t_sim_pnn=sim(net_pnn,test_p);
T_sim_pnn=vec2ind(t_sim_pnn);
result_pnn=T_sim_pnn';
%��������
% ��ȷ��accuracy
accuracy_grnn=length(find(result_grnn==test_T'))/length(test_T);
accuracy_pnn=length(find(result_pnn==test_T'))/length(test_T);
%����Ա�
result=[test_T' result_grnn result_pnn]
%��ͼ
figure(1)
plot(1:30,test_T,'bo',1:30,result_grnn,'r-*',1:30,result_pnn,'k:^')
grid on
xlabel('���Լ��������')
ylabel('���Լ��������')
string={'���Լ�Ԥ�����Ա�(GRNN vs PNN)';['��ȷ��:' num2str(accuracy_grnn*100) '%(GRNN) vs ' num2str(accuracy_pnn*100) '%(PNN)']};
title(string)
legend('��ʵֵ','GRNNԤ��ֵ','PNNԤ��ֵ')
