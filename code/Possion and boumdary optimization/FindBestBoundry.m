function OutputData2= FindBestBoundry(inputData00,DataMask00, inputData10,BoundryType)
%%%函数对掩膜文件DataMask0进行边界调整，找到两景影像散度差异更小的新的边界点
%%%输入参数inputData0、DataMask0和inputData1，均为uint8类型，输入也为uint8型
%%%输出参数是调整边界后的mask，unit8型
%%%BoundryType输入1或2或3。1表示基于像素值差异来找最佳边界点，2表示基于原影像梯度值找最佳边界点，3表示基于散度值差异来找最佳边界点
BufferSize=3;
%% 提前扩充，防止云在影像边界的情况
[nh,nw,b]=size(inputData00);
inputData0=zeros(nh+20,nw+20,b)-1;inputData0(11:end-10,11:end-10,:)=inputData00;
inputData1=zeros(nh+20,nw+20,b)-1;inputData1(11:end-10,11:end-10,:)=inputData10;
DataMask0=zeros(nh+20,nw+20,1);DataMask0(11:end-10,11:end-10,:)=DataMask00;
%% 计算两个影像散度差异
dif_ImgDivergence=double(inputData0);
if(BoundryType==3)%计算两个影像散度差异
    inputData0=imfilter(double(inputData0), [0,1,0;1,-4,1;0,1,0]);
    inputData1=imfilter(double(inputData1), [0,1,0;1,-4,1;0,1,0]);
    dif_ImgDivergence=abs(inputData0-inputData1);
    dif_ImgDivergence=sum(dif_ImgDivergence,3);
    DataMask0=imfilter(DataMask0,[0,1,0;1,1,1;0,1,0]);%扩充操作：mask往外扩充1。边界点散度不可靠
elseif (BoundryType==2)%原影像散度
    inputData0=abs(imfilter(double(inputData0), [0,1,0;1,-4,1;0,1,0]));
    dif_ImgDivergence=sum(inputData0,3);
    DataMask0=imfilter(DataMask0,[0,1,0;1,1,1;0,1,0]);%扩充操作：mask往外扩充1。边界点梯度不可靠
elseif(BoundryType==1)%计算两个影像像素值差异
    dif_ImgDivergence=abs(double(inputData0)-double(inputData1));
    dif_ImgDivergence=sum(dif_ImgDivergence,3);
else
    disp('请确定BoundryType为1、2或3！')
end
% abc=dif_ImgDivergence;
% figure,imshow(DataMask0)
%% 获取原始边界点
OriBoundryMask=imfilter(DataMask0,[0,1,0;1,1,1;0,1,0]);%扩充操作：mask往外扩充1
OriBoundryMask(OriBoundryMask>0)=1;
OriBoundryMask(DataMask0>0)=0;
%%  创建掩膜文件缓冲区
DataMask=imfilter(DataMask0,ones(BufferSize*2+1,BufferSize*2+1));%扩充操作：mask往外扩充1
DataMask(DataMask>0)=1;
DataMask(DataMask0>0)=0;
dif_ImgDivergence(DataMask==0)=-1;
%% 确定mask文件缓冲区内的边界点
[Array_H,Array_W]=find(DataMask>0);%缓冲区内所有像素点
DataMask(DataMask0>0)=1;
index_Boundry=zeros(length(Array_W),2);
num_Boundry=0;
for i=1:length(Array_W)
%     if (Array_H(i)-1<=0 || Array_W(i)-1<=0)
%         aaaaaa=0;
%     end
   if ( DataMask(Array_H(i)-1,Array_W(i))==0 || DataMask(Array_H(i)+1,Array_W(i))==0 ...
           || DataMask(Array_H(i),Array_W(i)-1)==0 || DataMask(Array_H(i),Array_W(i)+1)==0 )
       num_Boundry=num_Boundry+1;
       index_Boundry(num_Boundry,1)=Array_H(i);
       index_Boundry(num_Boundry,2)=Array_W(i);
   end  
end
index_Boundry=index_Boundry(1:num_Boundry,:);%缓冲区边界点
%% 在mask文件的缓冲区内找散度差异更小的边界点，重新确定边界
for iter=1:BufferSize-1
    %对原缓冲区的边界点更新，若更新后边界点的散度值小于之前的，则进行更新；否则保留
    temp_Divergence=ones(5,5);
    index_NewBoundry=zeros(length(Array_W),2);
    num_NewBoundry=0;
    for i=1:num_Boundry
        %以当前边界点为中心的三邻域
        temp_Divergence(2:4,2:4)=dif_ImgDivergence(index_Boundry(i,1)-1:index_Boundry(i,1)+1,index_Boundry(i,2)-1:index_Boundry(i,2)+1);
        s0=temp_Divergence(3,3);
        %在三邻域内对中心边界点进行更新。由于更新后原边界点仍为边界点，即除了中心边界点外，只可能新增其他的边界点不可能减少其他的边界点
        %故将中心点更新为-1后，判断其上下左右有没有增加的边界点，求其散度和即可。
        temp_Divergence(3,3)=-1;
        s1=0;n1=0;
        jj=[1,0;-1,0;0,1;0,-1];
        if(OriBoundryMask(index_Boundry(i,1),index_Boundry(i,2))==1)%如果该点本来就是原边界值
            s0=-1;
        else
            for j0=1:4
                if(temp_Divergence(3+jj(j0,1),3+jj(j0,2))~=-1 && ((temp_Divergence(4+jj(j0,1),3+jj(j0,2))==-1)+(temp_Divergence(2+jj(j0,1),3+jj(j0,2))==-1)+...
                        (temp_Divergence(3+jj(j0,1),4+jj(j0,2))==-1)+(temp_Divergence(3+jj(j0,1),2+jj(j0,2))==-1)==1))
                    n1=n1+1;
                    temp_Boundry(n1,:)=[index_Boundry(i,1)+jj(j0,1),index_Boundry(i,2)+jj(j0,2)];
                    s1=s1+temp_Divergence(3+jj(j0,1),3+jj(j0,2));
                end
            end
        end
        
        if(s1<=s0) %若调整后散度值<原散度值，则对边界进行更新
            dif_ImgDivergence(index_Boundry(i,1),index_Boundry(i,2))=-1;
            DataMask(index_Boundry(i,1),index_Boundry(i,2))=0;
            if(n1==0) %%判断若中心点为非必要点，直接进行删除
                continue;
            end
            index_NewBoundry(num_NewBoundry+1:num_NewBoundry+n1,:)=temp_Boundry;
            num_NewBoundry=num_NewBoundry+n1;
        else %若调整后散度值>原散度值，则不进行更改
            num_NewBoundry=num_NewBoundry+1;
            index_NewBoundry(num_NewBoundry,:)=index_Boundry(i,:);
        end
        clear temp_Boundry;
    end
    index_NewBoundry=index_NewBoundry(1:num_NewBoundry,:);
    index_Boundry=index_NewBoundry;num_Boundry=num_NewBoundry;
%     iter=iter+1;
end
DataMask(DataMask>0)=255;
OutputData2=uint8(DataMask(11:end-10,11:end-10,:));
% figure,imshow(DataMask)
%%新边界只能在原边界上增加，不能减小。验证一下
DataMask0(DataMask0~=0)=255;
DataMask(DataMask~=0)=255;
aa=DataMask-DataMask0;
% bb=zeros(size(aa));bb(aa<0)=255;figure,imshow(bb)
if(length(find(aa<0))>0)
    disp('新边界只能在原边界上增加，不能减小！！！');
end
end