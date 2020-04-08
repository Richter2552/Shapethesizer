sh=imread('1_strict_box_4x4_0_1.bmp');
gr=rgb2gray(sh);
dbl=im2double(gr);
vhs3=dbl;% matrix for the shape. Measures like corr2, immse etc can be applied on it. 
dbl=imresize(dbl,20);%scale the image to 8 times for better resolution
[m,n]=size(dbl);
dbl=~dbl;

%for creating shapes with grating
dbl=imrotate(dbl,210,'bicubic','crop');
grat=zeros(m,n);
cols = 0:n-1;
grat(:,find(mod(cols,640)<320))=1;
%bg=imrotate(grat,45,'bicubic','crop');
%imshow(bg);
bg=grat;
fig=bg.*dbl;
fig(:)=~fig;
%fig=imrotate(fig,-45,'bicubic','crop');
fig(fig==1)=0.5;
figure,imshow(fig)

% for creating colored shapes
bg=ones(m,n);
fig=bg.*dbl;
fig(:)=~fig;
figure,imshow(fig)
map = [0 1 1
       1 1 1];
colormap(map);
axis off


%for using figure as surface of the shapes
img=imread('craig.jpg');
img=rgb2gray(img);
img=im2double(img);
bg=imresize(img,[m,n]);
fig=bg.*dbl;
fig(fig==0)=1;
figure,imshow(fig)

