close all
clear all
clc
%车牌样本读取
img=imread('1.jpg');

%灰度转化，加权平均算法进行灰度处理：
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
title('灰度图像')
%[m,n]=size(img_bw);
% for i=1:m
%     for j=1:n
%        x(i,j)=img(i,j,1);
figure(2)
img_2=im2bw(img)
imshow(img_2);
title('二值化后图像')
% if c==3
%    fid=fopen('小狗.txt','w');
%     for i=1:m
%         for j=1:n
%             fprintf(fid,'%d         %d         %d\n',img(i,j,1),img(i,j,2),img(i,j,3));
%         end
%         fprintf(fid,'\n');
%     end
%     fclose(fid);
% end