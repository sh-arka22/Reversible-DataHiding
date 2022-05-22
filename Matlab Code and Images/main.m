clc;
clear all;
close all;
prompt = 'What should be the left threshold for 1st level RDH?(Enter a -ve value between -1 and -30) ';
Tl = input(prompt)
% Tl=-1;
prompt = 'What should be the right threshold for 1st level RDH?(Enter a value between 1 and 30) ';
Tr = input(prompt)
% Tr=1;
Tl_r=Tl;
Tr_r=Tr;
%% input image
input=imread('lena.jpg');
figure,imshow(input);title('input image');
figure,imhist(input);title('input image histogram');
%% IWT
[xdec1,ss1, ds1, sd1, dd1]=HaarFwd(input);
%% taking only the high frequency bands
ds1_11=ds1(1:64,1:128);
ds1_12=ds1(65:128,1:128);
I1dec1=[ds1_11,ds1_12;sd1,dd1];
% I_1dec1=I1dec1(:);
% [N_1,edges_1]=histcounts(I_1dec1);
% figure,histogram(I_1dec1,edges_1);
% axis([-200 200 0 50000]);
%% embedding
%% histogram shifting for data hiding
[m1,n1]=size(I1dec1);
% Tl=-1;
% Tr=1;
cap1=length(I1dec1(I1dec1>=Tl & I1dec1<Tr));  %number of bits that can be embedded
dt1=randi([0 1],1,cap1);       %data to be hidden
k1=1;
for i=1:m1
    for j=1:n1
        if (I1dec1(i,j)>=Tl)&&(I1dec1(i,j)<Tr)
            y1(i,j)=2*I1dec1(i,j)+dt1(k1);
            k1=k1+1;
        elseif I1dec1(i,j)>=Tr
            y1(i,j)=I1dec1(i,j)+Tr;
        elseif I1dec1(i,j)<Tl
            y1(i,j)=I1dec1(i,j)+Tl;
        end
    end
end
%% rearranging
y1_1=y1(1:64,1:128);
y2_1=y1(1:64,129:256);
y3_1=y1(65:192,1:128);
y4_1=y1(65:192,129:256);
y5_1=[y1_1;y2_1];
I2dec1=zeros(256,256);
I2dec1=[ss1,y5_1;y3_1,y4_1];
figure ,
subplot(2,1,1)
imshow(uint8(ss1));title('Decomposed Image');
subplot(2,1,2)
imhist(uint8(ss1));title('Decomposed Image Histogram');

%% Inverse IWT
[stego_1]=HaarInv(ss1, y5_1, y3_1, y4_1);
[op1,index1,p11,q11]=preventingOverflow(stego_1);
stego1=op1;
figure,imshow(uint8(stego1));title('stego image');
figure,imhist(uint8(stego1));title('stego image histogram');
%% Encryption level 1
[stego1_row,stego1_col]=size(stego1);

elliptic_matrix=[[1 6];[1 25];[3 8];[3 23];[4 3];[4 28];[5 3];[5 28];[6 15];[6 16];[9 11];[9 20];[12 10];[12 21];[14 8];[14 23];[15 13];[15 18];[17 2];[17 29];[18 5];[18 26];[20 5];[20 26];[21 4];[21 27];[22 3];[22 28];[23 14];[23 17];[24 5];[24 26];[26 11];[26 20];[27 11];[27 20];[28 2];[28 29];[30 1];[30 30]];


 Generator_point=6;
private_key_A=1;
n_A=private_key_A;
private_key_B=2;

n_B=private_key_B;
a=1;p=31;
b=3;

x1=elliptic_matrix(Generator_point,1);
y1=elliptic_matrix(Generator_point,2);

public_key_A=multell([x1 y1],n_A,a,b,p);
public_key_B=multell([x1 y1],n_B,a,b,p);
intial_key_K1=multell(public_key_B,n_A,a,b,p);
intial_key_K2=multell(public_key_A,n_B,a,b,p);

K1=multell([x1 y1],intial_key_K1(1,1),a,b,p);
K2=multell([x1 y1],intial_key_K1(1,2),a,b,p);

K11=[K1;K2];
K12=eye(2)-K11;
K12=mod(K12,256);
K21=eye(2)+K11;
K22=-K11;
K22=mod(K22,256);

K=[K11 K12;K21 K22];
%------------------------------------------------
%Encryption
%K=[4 28 253 228;15 18 241 239;5 28 252 228;15 19 241 238];
encrypted1=zeros(stego1_row,stego1_col);
for i=1:stego1_row
     l=1;
    for j=1:4:stego1_col
        P_1=double(stego1(i,j:j+3));
        P_1=P_1';
        C_1=K * P_1;
        C_1=mod(C_1,256);
        encrypted1(l:l+3,i)=C_1;
        l=l+4;
    end
end
encrypted_11z=encrypted1';
%% Encryption level 2
encrypted_1z=encrypted_11z(:);
p=3.628;
k(1)=0.632;
for n=1:stego1_row*stego1_col-1
    k(n+1) = cos(p*acos(k(n)));
end
k = abs(round(k*255));
k=k';
key = bitxor(uint8(k),uint8(encrypted_1z));
encrypted = reshape(key,[stego1_row stego1_col]);
% encrypted=uint8(encrypted);
figure,imshow(uint8(encrypted));title('encrypted image');
figure,imhist(uint8(encrypted));title('encrypted image histogram');


%% Decryption level 2
encrypted_z=encrypted(:);
decrypted = bitxor(uint8(k),uint8(encrypted_z));
decrypted = reshape(decrypted,[stego1_row stego1_col]);
isequal(decrypted,encrypted_11z);
%% Decryption level 1
decrypted_2=zeros(stego1_row,stego1_col);
for i=1:stego1_row
     l=1;
    for j=1:4:stego1_col
        P_1=double(decrypted(i,j:j+3));
        P_1=P_1';
        C_1=K * P_1;
        C_1=mod(C_1,256);
        decrypted_2(l:l+3,i)=C_1;
        l=l+4;
    end
end
decrypted_2=decrypted_2';
% decrypted_2=uint8(decrypted_2);
isequal(decrypted_2,stego1);
%% extracting
%% 1st level extraction
x5=overflowCompensation(decrypted_2,index1,p11,q11);
%% IWT
[xdec1_r,ss1_r, ds1_r, sd1_r, dd1_r]=HaarFwd(x5);
%% taking only the high frequency bands
ds1_11_r=ds1_r(1:64,1:128);
ds1_12_r=ds1_r(65:128,1:128);
r_I1dec1=[ds1_11_r,ds1_12_r;sd1_r,dd1_r];
[r_m1,r_n1]=size(r_I1dec1);
%% histogram shifting
Tl_r=-1;
Tr_r=1;
k=1;
z=r_I1dec1;
for i=1:r_m1
    for j=1:r_n1
        if (r_I1dec1(i,j)>=2*Tl_r)&&(r_I1dec1(i,j)<2*Tr_r)
            r_b(k)=mod(r_I1dec1(i,j),2);
            z(i,j)=floor(r_I1dec1(i,j)/2);
            k=k+1;
        elseif r_I1dec1(i,j)>=2*Tr_r
            z(i,j)=r_I1dec1(i,j)-Tr_r;
        elseif r_I1dec1(i,j)<2*Tl_r
            z(i,j)=r_I1dec1(i,j)-Tl_r;
        end
    end
end
isequal(dt1,r_b);
%% rearranging
y1_1=z(1:64,1:128);
y2_1=z(1:64,129:256);
y3_1=z(65:192,1:128);
y4_1=z(65:192,129:256);
y5_1=[y1_1;y2_1];
I2dec1_r=zeros(256,256);
I2dec1_r=[ss1_r,y5_1;y3_1,y4_1];
%% Inverse IWT
[retrieved]=HaarInv(ss1_r, y5_1, y3_1, y4_1);
isequal(uint8(retrieved),input);
% retrieved=uint8(retrieved);
figure,imshow(uint8(retrieved));title('retrieved image');
figure,imhist(uint8(retrieved));title('retrieved image histogram');


%% Result Analysis
%PSNR and SSIM of input and data hidden image
ssim_main_rdh=ssim(uint8(input),uint8(stego1))    %ssim
psnr_main_rdh=psnr(uint8(input),uint8(stego1))    %psnr
% embedding capacity of data hidden image
EC1=cap1/(256*256);
%%------------
%Entropy
E=entropy(encrypted);
disp('Entropy:')
disp(E)
%%---------------------------------------------
%NPCR
sum=0;
for i=1:stego1_row
    for j=1:stego1_col
        if(stego1(i,j)==encrypted(i,j))
           sum=sum+0;
        else
            sum=sum+1;
        end
        
    end
end
NPCR = (sum/(stego1_col*stego1_row));
disp('NPCR')
disp(NPCR)
%------------------------------
%UACI  Unified Average Changing Intensity
sum=0;
for i=1:stego1_row
    for j=1:stego1_col
        k=abs(double(stego1(i,j))-double(encrypted(i,j)));
        
        sum=sum+k/255;
    end
end
UACI=((sum/(256*256)))*100;
disp('UACI')
disp(UACI)
%------------------------------
%------------------------------
%PSNR Peak Signal Noise Ratio
sum=0;
for i=1:stego1_row
    for j=1:stego1_col
        k=(double(stego1(i,j))-double(encrypted(i,j)));
        sum=sum+k*k;
    end
end
MSE=((sum/(256*256)));

PSNR=20*log10(255)-10*log10(MSE);

disp('PSNR')
disp(PSNR);
%AdjancyCorrPixelRand

AdjancyCorrPixelRand(stego1,encrypted);

%------------------------------
%------------------------------

%correlation
  %Correlation coeff Original image horizontal direction
 I=stego1;
xi = I(:,1:end-1);  % original image
yi = I(:,2:end);  % original image
randIndex = randperm(numel(xi));                                   
randIndex = randIndex(1:2000);   
xRandi = xi(randIndex);            
yRandi = yi(randIndex); 
r_xyi = corr2(xRandi(:),yRandi(:));
disp('Correlation coeff Original image horizontal direction:')
disp(r_xyi);
% R=corrcoef(xRandi(:),yRandi(:));
% R_soh=R(2)^2;
figure(2)
subplot(3,2,1)
scatter(xRandi(:),yRandi(:),'.');
title('Scatter plot in horizontal direction for original image')
%Correlation coeff encrypted image horizontal direction
x = encrypted(:,1:end-1);  
y = encrypted(:,2:end);  
randIndex = randperm(numel(x));                               
randIndex = randIndex(1:2000);   
xRand = x(randIndex);  
yRand = y(randIndex);
% R=corrcoef(xRand(:),yRand(:));
% R_seh=R(2)^2;
% disp(R_seh);
r_xy = corr2(xRand(:),yRand(:));
disp('Correlation coeff encrypted image horizontal direction:')
disp(r_xy);
subplot(3,2,2)
scatter(xRand(:),yRand(:),'.');
title('Scatter plot in horizontal direction for encrypted image')
%Correlation coeff Original image vert direction
xiv = I(1:end-1,1:end-1,1);  
yiv = I(2:end,2:end,1); 
randIndex = randperm(numel(xiv));                                   
randIndex = randIndex(1:2000);   
xRandiv = xiv(randIndex);            
yRandiv = yiv(randIndex); 
% R=corrcoef(xRandiv(:),yRandiv(:));
% R_sov=R(2)^2;
% disp(R_sov);
r_xyiv =corr2(xRandiv(:),yRandiv(:));
disp('Correlation coeff Original image vertical direction:')
disp(r_xyiv);
subplot(3,2,3)
scatter(xRandiv(:),yRandiv(:),'.');
title('Scatter plot in vertical direction for original image')
%Correlation coeff encrypted image vert direction
xv = encrypted(1:end-1,1:end-1,1);  
yv = encrypted(2:end,2:end,1); 
randIndex = randperm(numel(xv));                                   
randIndex = randIndex(1:2000);   
xRandv = xv(randIndex);            
yRandv = yv(randIndex); 
% R=corrcoef(xRandv(:),yRandv(:));
% R_sev=R(2)^2;
% disp(R_sev);
r_xyv = corr2(xRandv(:),yRandv(:));
disp('Correlation coeff encrypted image vertical direction:')
disp(r_xyv);
subplot(3,2,4)
scatter(xRandv(:),yRandv(:),'.');
title('Scatter plot in vertical direction for encrypted image')
%Correlation coeff original image diagonal direction
xid = I(2:end,1:end-1,1);  
yid = I(1:end-1,2:end,1);
randIndex = randperm(numel(xid));                                   
randIndex = randIndex(1:2000);   
xRandid = xid(randIndex);            
yRandid = yid(randIndex); 
% R=corrcoef(xRandid(:),yRandid(:));
% R_sod=R(2)^2;
% disp(R_sod);
r_xyid =corr2(xRandid(:),yRandid(:));
disp('Correlation coeff Original image diagonal direction:')
disp(r_xyid);
subplot(3,2,5)
scatter(xRandid(:),yRandid(:),'.');
 title('Scatter plot in diagonal direction for original image')
%Correlation coeff encrypted image diagonal direction
xd = encrypted(2:end,1:end-1,1);  
yd = encrypted(1:end-1,2:end,1); 
randIndex = randperm(numel(xd));                                   
randIndex = randIndex(1:2000);   
xRandd = xd(randIndex);            
yRandd = yd(randIndex); 
%  R=corrcoef(xRandd(:),yRandd(:));
% R_sed=R(2)^2;
% disp(R_sed);
r_xyd = corr2(xRandd(:),yRandd(:));
disp('Correlation coeff encrypted image diagonal direction:')
subplot(3,2,6)
disp(r_xyd);
scatter(xRandd(:),yRandd(:),'.');
title('Scatter plot in diagonal direction for encrypted image')

%%-------------------------------
%Bits per pixel(BPP)

bpp=cap1/(128*128);

disp('BPP')
disp(bpp)
%%---------------------------------------------