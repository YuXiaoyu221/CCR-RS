function OutputData2= Weighted_Linear_Regression_gray_2(inputData0_0,DataMask0_0, inputData1_0)
% figure,imshow(inputData1)
NR=80;%参考相似像素数目
swindow=30;%起始窗口数目，大小为2*swindow+1
max_swindow=180;%最大窗口数目
alpha=0.1;%小值
[nh,nw,nb]=size(inputData0_0);
inputData0_0=double(inputData0_0);inputData1_0=double(inputData1_0);
inputData0=ones(nh+max_swindow*2,nw+max_swindow*2,nb)*500;
inputData0(max_swindow+1:max_swindow+nh,max_swindow+1:max_swindow+nw,:)=inputData0_0;
inputData1=ones(nh+max_swindow*2,nw+max_swindow*2,nb)*500;
inputData1(max_swindow+1:max_swindow+nh,max_swindow+1:max_swindow+nw,:)=inputData1_0;
DataMask0=zeros(nh+max_swindow*2,nw+max_swindow*2);
DataMask0(max_swindow+1:max_swindow+nh,max_swindow+1:max_swindow+nw,:)=DataMask0_0;

inputData0_g=double(rgb2gray(uint8(inputData0)));inputData1_g=double(rgb2gray(uint8(inputData1)));
inputData0=double(inputData0);inputData1=double(inputData1);
[ch,cw]=find(DataMask0>0);

DataMask11=DataMask0;
DataMask11(inputData0(:,:,3)>185)=255;

OutputData2=inputData0;
for ci=1:length(ch)        
    if (ch(ci)==272 && cw(ci)==355)
            aaa=0;
     end
    %%确定相似像素
    bs=swindow;
    refi=inputData1_g(ch(ci),cw(ci),:);%参考图像对应位置像素值 
    patch_ref=inputData1_g(ch(ci)-bs:ch(ci)+bs,cw(ci)-bs:cw(ci)+bs,:);%参考图像,初始搜索窗口，大小为2*swindow+1
    patch_ori=inputData0_g(ch(ci)-bs:ch(ci)+bs,cw(ci)-bs:cw(ci)+bs,:);%原始图像,初始搜索窗口
    patch_mask=DataMask11(ch(ci)-bs:ch(ci)+bs,cw(ci)-bs:cw(ci)+bs,:);%原始图像,初始搜索窗口
    patch_ref_t=patch_ref(patch_mask==0);
    T=std(patch_ref_t)+0.01;%%确定阈值，用来判断相似像素
    dif=abs(patch_ref-refi);%%判断相似像素：差值小于阈值
    dif(patch_mask>0)=2*T;
    [sh,sw]=find(dif<=T);%%相似像素
    while ((bs<max_swindow) &&(length(sh)<NR || length(patch_ref_t)<NR))%%相似像素数目不满足要求，扩大搜索范围
        bs=bs+2;
        refi=inputData1_g(ch(ci),cw(ci),:);%参考图像对应位置像素值 
        patch_ref=inputData1_g(ch(ci)-bs:ch(ci)+bs,cw(ci)-bs:cw(ci)+bs,:);%参考图像,初始搜索窗口，大小为2*swindow+1
        patch_ori=inputData0_g(ch(ci)-bs:ch(ci)+bs,cw(ci)-bs:cw(ci)+bs,:);%原始图像,初始搜索窗口
        patch_mask=DataMask11(ch(ci)-bs:ch(ci)+bs,cw(ci)-bs:cw(ci)+bs,:);%原始图像,初始搜索窗口
        patch_ref_t=patch_ref(patch_mask==0);
        T=std(patch_ref_t)+0.01;%%确定阈值，用来判断相似像素
        dif=abs(patch_ref-refi);%%判断相似像素：差值小于阈值
        dif(patch_mask>0)=2*T;
        sh=find(dif<=T);%%相似像素
    end
    %%确定权重
    [ah,aw]=find(dif>-1);
    dif(dif>T)=-2*alpha;
    dif=reshape(dif,(2*bs+1)*(2*bs+1),1);
    Di=(dif+alpha).*((ah-bs-1).*(ah-bs-1)+(aw-bs-1).*(aw-bs-1));
    Di=Di(Di>0);
    iDi=1./Di;
    wi=iDi/(sum(iDi));
    %%取有效像素
for ib=1:nb
    patch_ori=inputData0(ch(ci)-bs:ch(ci)+bs,cw(ci)-bs:cw(ci)+bs,ib);
    patch_ref=inputData1(ch(ci)-bs:ch(ci)+bs,cw(ci)-bs:cw(ci)+bs,ib);
    patch_ori=reshape(patch_ori,(2*bs+1)*(2*bs+1),1);
    patch_ref=reshape(patch_ref,(2*bs+1)*(2*bs+1),1);
    patch_ori_t=patch_ori(dif>-alpha);
    patch_ref_t=patch_ref(dif>-alpha);
    %%最小二乘回归
    if length(sh)>=NR
        mOri=mean(patch_ori_t);
        Ori_mOri=patch_ori_t-mOri;
        mRef=mean(patch_ref_t);
        Ref_mRef=patch_ref_t-mRef;

        if (sum(abs(Ref_mRef))==0)
            ca=0;
        else
        ca=sum(wi.*Ori_mOri.*Ref_mRef)/(sum(wi.*Ref_mRef.*Ref_mRef));%%系数a
        end
%         ca=sum(wi.*Ori_mOri.*Ref_mRef)/sum(wi.*Ref_mRef.*Ref_mRef);%%系数a
        cb=mOri-ca*mRef;
    elseif length(sh)>=2
        mOri=mean(patch_ori_t);
        mRef=mean(patch_ref_t);
        ca=mOri/mRef;%%系数a
        cb=0;
    else
        ca=1;cb=0;
    end 
%     OutputData2(ch(ci),cw(ci),:)=ca*inputData1(ch(ci),cw(ci),:)+cb;
    OutputData2(ch(ci),cw(ci),ib)=ca*inputData1(ch(ci),cw(ci),ib)+cb;
end
%     if(sum(abs(OutputData2(ch(ci),cw(ci),:)-inputData1(ch(ci),cw(ci),:)))>60*3 )
%         OutputData2(ch(ci),cw(ci),:)=inputData1(ch(ci),cw(ci),:);
%     end
end
% OutputData2=uint8(OutputData2);
OutputData2=uint8(OutputData2(max_swindow+1:max_swindow+nh,max_swindow+1:max_swindow+nw,:));
end