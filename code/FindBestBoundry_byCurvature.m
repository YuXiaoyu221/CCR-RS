function OutputData2= FindBestBoundry_byCurvature(inputData00,DataMask00, inputData10, BufferSize)
%%% the boundary for cloud coverage area eplacement is optimized before curvature reconstruction, 
%%%                                      to make the boundary pass through the pixels with minimum curvature difference. 
%%% the input data of inputData0¡¢DataMask0 and inputData1 are all uint8 type
%%% the output is updated mask(unit8),where the value 255 indicates the missing area, and the value 0 is the clear area.
% BufferSize=6;
Fsize=max(BufferSize*2+4,20);
%% image expansion
[nh,nw,b]=size(inputData00);
inputData0=zeros(nh+Fsize,nw+Fsize,b)-1;inputData0(Fsize/2+1:end-Fsize/2,Fsize/2+1:end-Fsize/2,:)=inputData00;
inputData1=zeros(nh+Fsize,nw+Fsize,b)-1;inputData1(Fsize/2+1:end-Fsize/2,Fsize/2+1:end-Fsize/2,:)=inputData10;
DataMask0=zeros(nh+Fsize,nw+Fsize,1);DataMask0(Fsize/2+1:end-Fsize/2,Fsize/2+1:end-Fsize/2,:)=DataMask00;
%% calculate the curvature differences
inputData0=imfilter(double(inputData0), [-1/16,5/16,-1/16;5/16,-1,5/16;-1/16,5/16,-1/16]);
inputData1=imfilter(double(inputData1), [-1/16,5/16,-1/16;5/16,-1,5/16;-1/16,5/16,-1/16]);
dif_ImgDivergence=abs(inputData0-inputData1);
dif_ImgDivergence=sum(dif_ImgDivergence,3);
DataMask0=imfilter(DataMask0,[1,1,1;1,1,1;1,1,1]);
%% Get the initial boundary pixels
OriBoundryMask=imfilter(DataMask0,[1,1,1;1,1,1;1,1,1]);
OriBoundryMask(OriBoundryMask>0)=1;
OriBoundryMask(DataMask0>0)=0;
%% Creat the buffer area
DataMask=imfilter(DataMask0,ones(BufferSize*2+1,BufferSize*2+1));
DataMask(DataMask>0)=1;
DataMask(DataMask0>0)=0;
dif_ImgDivergence(DataMask==0)=-1;
%% determine all the pixels in the buffer area
[Array_H,Array_W]=find(DataMask>0);
DataMask(DataMask0>0)=1;
index_Boundry=zeros(length(Array_W),2);
num_Boundry=0;
for i=1:length(Array_W)
    if (sum(sum(DataMask(Array_H(i)-1:Array_H(i)+1,Array_W(i)-1:Array_W(i)+1)))<9)
%    if ( DataMask(Array_H(i)-1,Array_W(i))==0 || DataMask(Array_H(i)+1,Array_W(i))==0 ...
%            || DataMask(Array_H(i),Array_W(i)-1)==0 || DataMask(Array_H(i),Array_W(i)+1)==0 )
       num_Boundry=num_Boundry+1;
       index_Boundry(num_Boundry,1)=Array_H(i);
       index_Boundry(num_Boundry,2)=Array_W(i);
   end  
end
index_Boundry=index_Boundry(1:num_Boundry,:);
%%  determine the optimal boundary with the minimum curvature cost
for iter=BufferSize-4:BufferSize
    % the initial boundary pixels are updated when the divergence of pixels get smaller 
    temp_Divergence=ones(5,5);
    index_NewBoundry=zeros(length(Array_W),2);
    num_NewBoundry=0;
    for i=1:num_Boundry
        % the 3*3 neighber
        temp_Divergence(2:4,2:4)=dif_ImgDivergence(index_Boundry(i,1)-1:index_Boundry(i,1)+1,index_Boundry(i,2)-1:index_Boundry(i,2)+1);
        s0=temp_Divergence(3,3);
        temp_Divergence(3,3)=-1;
        s1=0;n1=0;
        jj=[1,1;1,0;1,-1;0,1;0,-1;-1,1;-1,0;-1,-1];
        if(OriBoundryMask(index_Boundry(i,1),index_Boundry(i,2))==1)%if the pixel is belong to initial boundary
            s0=-1;
        else
            for j0=1:8
                if(temp_Divergence(3+jj(j0,1),3+jj(j0,2))~=-1 && ((temp_Divergence(3+jj(j0,1),4+jj(j0,2))==-1)+(temp_Divergence(3+jj(j0,1),2+jj(j0,2))==-1)+ ......
                        (temp_Divergence(4+jj(j0,1),4+jj(j0,2))==-1)+(temp_Divergence(4+jj(j0,1),3+jj(j0,2))==-1)+(temp_Divergence(4+jj(j0,1),2+jj(j0,2))==-1)+ ...
                        (temp_Divergence(2+jj(j0,1),4+jj(j0,2))==-1)+(temp_Divergence(2+jj(j0,1),3+jj(j0,2))==-1)+(temp_Divergence(2+jj(j0,1),2+jj(j0,2))==-1)==1))
                    n1=n1+1;
                    temp_Boundry(n1,:)=[index_Boundry(i,1)+jj(j0,1),index_Boundry(i,2)+jj(j0,2)];
                    s1=s1+temp_Divergence(3+jj(j0,1),3+jj(j0,2));
                end
            end
        end
        
        if(s1<=s0) %when the Divergence < original£¬update it
            dif_ImgDivergence(index_Boundry(i,1),index_Boundry(i,2))=-1;
            DataMask(index_Boundry(i,1),index_Boundry(i,2))=0;
            if(n1==0)
                continue;
            end
            index_NewBoundry(num_NewBoundry+1:num_NewBoundry+n1,:)=temp_Boundry;
            num_NewBoundry=num_NewBoundry+n1;
        else %when the Divergence > original£¬keep it
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
OutputData2=uint8(DataMask(Fsize/2+1:end-Fsize/2,Fsize/2+1:end-Fsize/2,:));
%% the new boundary pixels must be more that original
DataMask0(DataMask0~=0)=255;
DataMask(DataMask~=0)=255;
aa=DataMask-DataMask0;
if(length(find(aa<0))>0)
    disp('the new boundary pixels must be more that original £¡£¡£¡');
end
end