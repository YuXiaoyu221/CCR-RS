function [PSNR,SSIM,M_SSIM,MS_SSIM,IW_SSIM,FD,CE,CC]= ResultEvaluation(OutData,DataMask0, RefData)
%%%函数基于原始无云参考影像，对去云结果进行评价，输出指标结果
%%%输入参数OutData,DataMask0,RefData，均为uint8类型，输入也为uint8型
%%%OutData是去云后结果,DataMask0云掩膜文件, RefData是原始无云参考影像
ibands=size(OutData,3);
OutData=double(OutData);
RefData=double(RefData);
DataMask0(DataMask0>0)=1;
DataMask=repmat(DataMask0,1,1,ibands);
ssimWindows=3;%%SSIM滑动窗口大小，3,5,7,9....奇数
% %-----------波段均值------------//据说CE交叉熵要将RGB转灰度影像计算，即rgb2gray()函数
% if(ibands==3)
OutData=double(rgb2gray(uint8(OutData)));
RefData=double(rgb2gray(uint8(RefData)));
ibands=1;
% end
DataMask=repmat(DataMask0,1,1,ibands);
% OutData1=double(OutData(:,:,1));
% RefData1=double(RefData(:,:,1));
% for ib=2:ibands
%     OutData1=OutData1+double(OutData(:,:,ib));
%     RefData1=RefData1+double(RefData(:,:,ib));
% end
% OutData=OutData1/ibands;
% RefData=RefData1/ibands;
% ibands=1;
% DataMask=repmat(DataMask0,1,1,ibands);
% clear OutData1 RefData1;
% %-----------波段均值------------
OutData1=OutData(DataMask>0);
RefData1=RefData(DataMask>0);
OutData1=reshape(OutData1,length(OutData1)/ibands,ibands);
RefData1=reshape(RefData1,length(RefData1)/ibands,ibands);%仅剩云像素参与计算
CloudsNum=length(find(DataMask0>0));
%% 峰值信噪比PSNR
MSE=double(OutData1)-double(RefData1);
MSE=MSE.*MSE;
MSE=sum((MSE))/CloudsNum;
PSNR=20*log10(255./sqrt(MSE));
PSNR=sum(PSNR)/ibands;
clear MSE;
%% 结构相似性指数SSIM
C1=(0.01*255)*(0.01*255);
C2=(0.03*255)*(0.03*255);
x=OutData1;
y=RefData1;
Ux=mean(x);
Uy=mean(y);
STDx2=Ux;STDy2=Ux;STDxy=Ux;SSIM=Ux;
for ib=1:ibands
    aa=abs(cov(x(:,ib),y(:,ib)));
    STDx2(1,ib)=aa(1,1);
    STDy2(1,ib)=aa(2,2);
    STDxy(1,ib)=aa(1,2);
    SSIM(1,ib)=((2*Ux(1,ib)*Uy(1,ib)+C1)*(2*STDxy(1,ib)+C2))/((Ux(1,ib)*Ux(1,ib)+Uy(1,ib)*Uy(1,ib)+C1)...
        *(STDx2(1,ib)+STDy2(1,ib)+C2));
end
SSIM=sum(SSIM)/ibands;
clear x y Ux Uy STDx2 STDy2 STDxy aa;
%% 平均结构相似性指数 M-SSIM
MatSSIM=zeros(size(OutData));
sWindows1=(ssimWindows-1)/2;
counts=0;
C1=(0.01*255)*(0.01*255);
C2=(0.03*255)*(0.03*255);
for i=1+sWindows1:size(OutData,1)-sWindows1
    for j=1+sWindows1:size(OutData,2)-sWindows1
        w_mask=DataMask0(i-sWindows1:i+sWindows1,j-sWindows1:j+sWindows1);
        if(sum(sum(w_mask))<ssimWindows*sWindows1)
            continue;%%如果当前窗口内原始云像素过少，则不计算（因为数量过少不具有一般性）
        end
        w_mask=repmat(w_mask,1,1,ibands);
        x=OutData(i-sWindows1:i+sWindows1,j-sWindows1:j+sWindows1,:);
        y=RefData(i-sWindows1:i+sWindows1,j-sWindows1:j+sWindows1,:);
        x=x(w_mask>0);x=reshape(x,length(x)/ibands,ibands);
        y=y(w_mask>0);y=reshape(y,length(y)/ibands,ibands);
        Ux=mean(x);
        Uy=mean(y);
        STDx2=Ux;STDy2=Ux;STDxy=Ux;
        for ib=1:ibands
            aa=abs(cov(x(:,ib),y(:,ib)));
            STDx2(1,ib)=aa(1,1);
            STDy2(1,ib)=aa(2,2);
            STDxy(1,ib)=aa(1,2);
            MatSSIM(i,j,ib)=((2*Ux(1,ib)*Uy(1,ib)+C1)*(2*STDxy(1,ib)+C2))/((Ux(1,ib)*Ux(1,ib)+Uy(1,ib)*Uy(1,ib)+C1)...
                *(STDx2(1,ib)+STDy2(1,ib)+C2));
        end
        counts=counts+1;
    end
end
M_SSIM=sum(sum(sum(abs(MatSSIM))))/ibands/counts;
clear x y Ux Uy STDx2 STDy2 STDxy aa MatSSIM w_mask;
%% 多尺度结构相似性指数 MS-SSIM。转化为灰度图计算的-rgb2gray
Belta=[0.0448,0.2856,0.3001,0.2363,0.1333];
C1=(0.01*255)*(0.01*255);
C2=(0.03*255)*(0.03*255);
sWindows1=(ssimWindows-1)/2;
MS_SSIM0=zeros(1,5);
if ibands==3 
    OutDatai=double(rgb2gray(uint8(OutData)));
    RefDatai=double(rgb2gray(uint8(RefData)));
else
    OutDatai=double(OutData);
    RefDatai=double(RefData);
end
DataMaski=double(DataMask0);
for iter=1:5
    MatSSIM=zeros(size(OutDatai));
    counts=0;
    for i=1+sWindows1:size(OutDatai,1)-sWindows1
        for j=1+sWindows1:size(OutDatai,2)-sWindows1
            w_mask=DataMaski(i-sWindows1:i+sWindows1,j-sWindows1:j+sWindows1);
            if(sum(sum(w_mask))<ssimWindows*sWindows1)
                continue;%%如果当前窗口内原始云像素过少，则不计算（因为数量过少不具有一般性）
            end
            x=OutDatai(i-sWindows1:i+sWindows1,j-sWindows1:j+sWindows1);
            y=RefDatai(i-sWindows1:i+sWindows1,j-sWindows1:j+sWindows1);
            x=x(w_mask>0);
            y=y(w_mask>0);
            Ux=mean(x);
            Uy=mean(y);
            aa=abs(cov(x,y));
            STDx2=aa(1,1);
            STDy2=aa(2,2);
            STDxy=aa(1,2);
            if(iter==5)
                MatSSIM(i,j)=(2*STDxy+C2)/(STDx2+STDy2+C2);
            else
                MatSSIM(i,j)=((2*Ux*Uy+C1)*(2*STDxy+C2))/((Ux*Ux+Uy*Uy+C1)*(STDx2+STDy2+C2));
            end
            counts=counts+1;
        end
    end
    MS_SSIM0(1,iter)=sum(sum(abs(MatSSIM)))/counts;
    OutDatai=double(LowPassFiltering_DownSampling(OutDatai,80,0.5));
    RefDatai=double(LowPassFiltering_DownSampling(RefDatai,80,0.5));
    DataMaski=double(imresize(DataMaski,0.5));
end
MS_SSIM=1;
for iter=1:5
    MS_SSIM=MS_SSIM*(MS_SSIM0(1,iter)^Belta(1,iter));
end
clear x y Ux Uy STDx2 STDy2 STDxy aa MatSSIM w_mask;
clear OutDatai RefDatai DataMaski MS_SSIM0;
%% 信息加权结构相似度指数 IW-SSIM
%Wang, Z.; Li, Q. Information Content Weighting for Perceptual Image Quality Assessment.IEEE Trans. Image Process. 2011, 20, 1185–1198. [CrossRef]

IW_SSIM=0;
%% 图像清晰度 FD
FD_OutData_row=OutData;%%行向梯度值
for i=1:size(OutData,1)-1
    FD_OutData_row(i,:)=(OutData(i,:)-OutData(i+1,:)).*(OutData(i,:)-OutData(i+1,:));
end
FD_OutData_col=OutData;%%列向梯度值
for i=1:size(OutData,2)-1
    FD_OutData_col(:,i)=(OutData(:,i)-OutData(:,i+1)).*(OutData(:,i)-OutData(:,i+1));
end
FD_OutData=sqrt((FD_OutData_col+FD_OutData_row)/2);
FD_OutData(DataMask==0)=0;
FD=sum(sum(FD_OutData))/CloudsNum;
FD=sum(FD)/ibands;
clear FD_OutData_row FD_OutData_col FD_OutData;
%% 交叉熵 CE
Freq_OutData=zeros(256,ibands);
Freq_RefData=zeros(256,ibands);
for ib=1:ibands
    for i=1:CloudsNum
        Freq_OutData(OutData1(i,ib)+1,ib)=Freq_OutData(OutData1(i,ib)+1,ib)+1;%第1个存储的是像素值为0的像素个数
        Freq_RefData(RefData1(i,ib)+1,ib)=Freq_RefData(RefData1(i,ib)+1,ib)+1;%第256个存储的是像素值为255的像素个数
    end
end
Freq_OutData=Freq_OutData/CloudsNum;%%次数转化为频率
Freq_RefData=Freq_RefData/CloudsNum;
CE=zeros(1,ibands);
for ib=1:ibands
    for i=1:256
        if(Freq_OutData(i,ib)~=0 && Freq_RefData(i,ib)~=0)
%             CE(1,ib)=CE(1,ib)+Freq_OutData(i,ib)*log2(Freq_OutData(i,ib)/Freq_RefData(i,ib));
           CE(1,ib)=CE(1,ib)+Freq_RefData(i,ib)*log2(Freq_RefData(i,ib)/Freq_OutData(i,ib));
        end
    end
end
CE=sum(CE)/ibands;
%% 相关系数CC
CC=zeros(1,ibands);
for ib=1:ibands
    cc0=corrcoef(OutData1(:,ib),RefData1(:,ib));
    CC(1,ib)=cc0(1,2);
end
CC=sum(CC)/ibands;
end