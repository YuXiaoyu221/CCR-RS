function OutputData2= Curvature_based_Image_Reconstruction(inputData00,DataMask00, inputData10)
%%% the input data of inputData0¡¢DataMask0 and inputData1 are all uint8 type
[nh,nw,b]=size(inputData00);
inputData0=zeros(nh+2,nw+2,b)-1;inputData0(2:end-1,2:end-1,:)=inputData00;
inputData1=zeros(nh+2,nw+2,b)-1;inputData1(2:end-1,2:end-1,:)=inputData10;
DataMask0=DataMask00;
DataMask0=imfilter(DataMask00, [1,1,1;1,1,1;1,1,1]);%image expansion
DataMask0(DataMask0>0)=1;
DataMask=zeros(nh+2,nw+2,1)-1;DataMask(2:end-1,2:end-1,:)=DataMask0;
DataMask0=zeros(nh+2,nw+2,1)-1;DataMask0(2:end-1,2:end-1,:)=DataMask00;
[Array_H,Array_W]=find(DataMask>0);
%%creat matrix DiagA
%%% the value of boundary pixels in DiagA is 1
DiagA_i=zeros(1,length(Array_H));
DiagA_j=zeros(1,length(Array_H));
DiagA_v=ones(1,length(Array_H));
% the contection between DiagA with [Array_H,Array_W]
ArrayWH_index=zeros(max(Array_H),max(Array_W));
for i=1:length(Array_W)
    ArrayWH_index(Array_H(i),Array_W(i))=i;
end
%%% creat matrix MatB£¬calculate Divergence
OutputData=imfilter(double(inputData1), [-1/16,5/16,-1/16;5/16,-1,5/16;-1/16,5/16,-1/16]);%-OutputData--------------------------------------
num_nozeros=1;% the count of non-zero pixesl
for i=1:length(Array_W)
    if ( DataMask0(Array_H(i),Array_W(i))>0)%%unknown pixels
      if( DataMask(Array_H(i)-1,Array_W(i))>0 && DataMask(Array_H(i)+1,Array_W(i))>0 && DataMask(Array_H(i),Array_W(i)-1)>0 && DataMask(Array_H(i),Array_W(i)+1)>0 && ...
       DataMask(Array_H(i)-1,Array_W(i)-1)>0 && DataMask(Array_H(i)+1,Array_W(i)-1)>0 && DataMask(Array_H(i)-1,Array_W(i)+1)>0 && DataMask(Array_H(i)+1,Array_W(i)+1)>0 ) 
        DiagA_i(num_nozeros)=i; DiagA_j(num_nozeros)=i; DiagA_v(num_nozeros)=-1; num_nozeros=num_nozeros+1;
        DiagA_i(num_nozeros)=i; DiagA_j(num_nozeros)=ArrayWH_index(Array_H(i)-1,Array_W(i)-1); DiagA_v(num_nozeros)=-1/16; num_nozeros=num_nozeros+1;
        DiagA_i(num_nozeros)=i; DiagA_j(num_nozeros)=i-1; DiagA_v(num_nozeros)=5/16; num_nozeros=num_nozeros+1;
        DiagA_i(num_nozeros)=i; DiagA_j(num_nozeros)=ArrayWH_index(Array_H(i)-1,Array_W(i)+1); DiagA_v(num_nozeros)=-1/16; num_nozeros=num_nozeros+1;
        DiagA_i(num_nozeros)=i; DiagA_j(num_nozeros)=ArrayWH_index(Array_H(i),Array_W(i)-1); DiagA_v(num_nozeros)=5/16; num_nozeros=num_nozeros+1;    
        DiagA_i(num_nozeros)=i; DiagA_j(num_nozeros)=ArrayWH_index(Array_H(i),Array_W(i)+1); DiagA_v(num_nozeros)=5/16; num_nozeros=num_nozeros+1;
        DiagA_i(num_nozeros)=i; DiagA_j(num_nozeros)=ArrayWH_index(Array_H(i)+1,Array_W(i)-1); DiagA_v(num_nozeros)=-1/16; num_nozeros=num_nozeros+1;
        DiagA_i(num_nozeros)=i; DiagA_j(num_nozeros)=i+1; DiagA_v(num_nozeros)=5/16; num_nozeros=num_nozeros+1;
        DiagA_i(num_nozeros)=i; DiagA_j(num_nozeros)=ArrayWH_index(Array_H(i)+1,Array_W(i)+1); DiagA_v(num_nozeros)=-1/16; num_nozeros=num_nozeros+1;
      else %%when the missing pixels in the boundary
          s=0;ls=0;
          if(DataMask(Array_H(i)-1,Array_W(i)-1)>0)
              DiagA_i(num_nozeros)=i; DiagA_j(num_nozeros)=ArrayWH_index(Array_H(i)-1,Array_W(i)-1); DiagA_v(num_nozeros)=-1/16; num_nozeros=num_nozeros+1;
              s=s+DiagA_v(num_nozeros-1);
              ls=ls+inputData1(Array_H(i)-1,Array_W(i)-1,:)*DiagA_v(num_nozeros-1);
          end
          if(DataMask(Array_H(i)-1,Array_W(i))>0)
              DiagA_i(num_nozeros)=i; DiagA_j(num_nozeros)=i-1; DiagA_v(num_nozeros)=5/16; num_nozeros=num_nozeros+1;
              s=s+DiagA_v(num_nozeros-1);
              ls=ls+inputData1(Array_H(i)-1,Array_W(i),:)*DiagA_v(num_nozeros-1);
          end
          if(DataMask(Array_H(i)-1,Array_W(i)+1)>0)
              DiagA_i(num_nozeros)=i; DiagA_j(num_nozeros)=ArrayWH_index(Array_H(i)-1,Array_W(i)+1); DiagA_v(num_nozeros)=-1/16; num_nozeros=num_nozeros+1;
              s=s+DiagA_v(num_nozeros-1);
              ls=ls+inputData1(Array_H(i)-1,Array_W(i)+1,:)*DiagA_v(num_nozeros-1);
          end
          if(DataMask(Array_H(i),Array_W(i)-1)>0)
              DiagA_i(num_nozeros)=i; DiagA_j(num_nozeros)=ArrayWH_index(Array_H(i),Array_W(i)-1); DiagA_v(num_nozeros)=5/16; num_nozeros=num_nozeros+1; 
              s=s+DiagA_v(num_nozeros-1);
              ls=ls+inputData1(Array_H(i),Array_W(i)-1,:)*DiagA_v(num_nozeros-1);
          end
          if(DataMask(Array_H(i),Array_W(i)+1)>0)   
              DiagA_i(num_nozeros)=i; DiagA_j(num_nozeros)=ArrayWH_index(Array_H(i),Array_W(i)+1); DiagA_v(num_nozeros)=5/16; num_nozeros=num_nozeros+1;
              s=s+DiagA_v(num_nozeros-1);
              ls=ls+inputData1(Array_H(i),Array_W(i)+1,:)*DiagA_v(num_nozeros-1);
          end
          if(DataMask(Array_H(i)+1,Array_W(i)-1)>0)
              DiagA_i(num_nozeros)=i; DiagA_j(num_nozeros)=ArrayWH_index(Array_H(i)+1,Array_W(i)-1); DiagA_v(num_nozeros)=-1/16; num_nozeros=num_nozeros+1;
              s=s+DiagA_v(num_nozeros-1);
              ls=ls+inputData1(Array_H(i)+1,Array_W(i)-1,:)*DiagA_v(num_nozeros-1);
          end
          if(DataMask(Array_H(i)+1,Array_W(i))>0)
              DiagA_i(num_nozeros)=i; DiagA_j(num_nozeros)=i+1; DiagA_v(num_nozeros)=5/16; num_nozeros=num_nozeros+1;
              s=s+DiagA_v(num_nozeros-1);
              ls=ls+inputData1(Array_H(i)+1,Array_W(i),:)*DiagA_v(num_nozeros-1);
          end
          if(DataMask(Array_H(i)+1,Array_W(i)+1)>0)
              DiagA_i(num_nozeros)=i; DiagA_j(num_nozeros)=ArrayWH_index(Array_H(i)+1,Array_W(i)+1); DiagA_v(num_nozeros)=-1/16; num_nozeros=num_nozeros+1;  
              s=s+DiagA_v(num_nozeros-1);
              ls=ls+inputData1(Array_H(i)+1,Array_W(i)+1,:)*DiagA_v(num_nozeros-1);
          end
          DiagA_i(num_nozeros)=i; DiagA_j(num_nozeros)=i; DiagA_v(num_nozeros)=-s; num_nozeros=num_nozeros+1;
          OutputData(Array_H(i),Array_W(i),:)=ls-s*inputData1(Array_H(i),Array_W(i),:);
      end
    else %%known pixels
        DiagA_i(num_nozeros)=i; DiagA_j(num_nozeros)=i; DiagA_v(num_nozeros)=1; num_nozeros=num_nozeros+1;
        OutputData(Array_H(i),Array_W(i),:)=inputData0(Array_H(i),Array_W(i),:);
    end
end
DiagA=sparse(DiagA_i,DiagA_j,DiagA_v);
clear ArrayWH_index DiagA_i DiagA_j DiagA_v;
DataMask=repmat(DataMask,1,1,b);
OutputData1=OutputData(DataMask>0);
clear DataMask0 OutputData inputData1
OutputData1=reshape(OutputData1,length(OutputData1)/b,b);
OutputData1=double(OutputData1);
for ib=1:b
    setup = struct('type','ilutp','droptol',1e-6);
    [L,U] = ilu(sparse(DiagA),setup);
   OutputData1(:,ib)=bicgstab(DiagA,OutputData1(:,ib),1e-4,100,L,U);
end
OutputData2=inputData0;
OutputData2(DataMask>0)=abs(OutputData1);
OutputData2=OutputData2(2:end-1,2:end-1,:);
OutputData2=uint8(OutputData2);
end 