close all
clear all
clc
%����������ȡ
img=imread('1.jpg');

%�Ҷ�ת������Ȩƽ���㷨���лҶȴ���
[R,G,B]=size(img);
img_greymean=zeros(R,G);
img_greymean=uint8(img_greymean);
for i=1:R
    for j=1:G
        img_greymean(i,j)=sum(img(i,j,:))/3;
    end
end
figure(1)
% img_grey=rgb2gray(img);
% imshow(img_grey);
imshow(img_greymean);
title('�Ҷ�ͼ��')
%[m,n]=size(img_bw);
% for i=1:m
%     for j=1:n
%        x(i,j)=img(i,j,1);
figure(2)
img_2=im2bw(img)
imshow(img_2);
title('��ֵ����ͼ��')
% if c==3
%    fid=fopen('С��.txt','w');
%     for i=1:m
%         for j=1:n
%             fprintf(fid,'%d         %d         %d\n',img(i,j,1),img(i,j,2),img(i,j,3));
%         end
%         fprintf(fid,'\n');
%     end
%     fclose(fid);
% end