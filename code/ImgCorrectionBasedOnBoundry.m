function OutputData3=ImgCorrectionBasedOnBoundry(OutData0,A0,A_Mask0, A_Mask1,weights)
%%% refine initial cloud-free result based on the deviation between true value and reconstructed value of these checkpoints
%%% OutData0 is the initial cloud removal result of original cloudy image A0
%%% A_Mask0 is the original cloud mask; A_Mask1 is the updated cloud mask
OutData0=double(OutData0);A0=double(A0);
OutputData2=OutData0;
A_Mask1(A_Mask1>0)=1;
A_Mask11=A_Mask1;
A_Mask11(A_Mask0>0)=0;%
% A_Mask1=repmat(A_Mask1,1,1,size(OutData0,3));
for ib=1:size(OutData0,3)
    OutData0_ib=OutData0(:,:,ib);
    OutData0_known=OutData0_ib(A_Mask11>0);% the checkpoints after reconstructed
    A0_ib=A0(:,:,ib);
    A0_known=A0_ib(A_Mask11>0);% the original values of these checkpoints
    p=polyfit(double(OutData0_known),double(A0_known),3);
    yi=polyval(p,double(OutData0_ib(A_Mask1>0)));
    OutputData2_ib=A0_ib;
    OutputData2_ib(A_Mask1>0)=yi;
%     yi=polyval(p,double(OutData0_ib(A_Mask0>0)));
%     OutputData2_ib=A0_ib;
%     OutputData2_ib(A_Mask0>0)=yi;
    OutputData2(:,:,ib)=OutputData2_ib;
end
OutputData2=double(weights(1)*double(OutputData2)+weights(2)*double(OutData0));

A_Mask03=repmat(A_Mask0,1,1,size(OutData0,3));
OutputData3=OutData0;
OutputData3(A_Mask03>0)=OutputData2(A_Mask03>0);
if(isa(OutData0,'uint8'))
    OutputData3=uint8(OutputData3);
end
end