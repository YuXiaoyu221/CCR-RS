function OutputData2= ImgCorrectionBasedOnBoundry_selected(OutData0,A0,A_Mask0, A_Mask1,weights)
%%% refine initial cloud-free result based on the deviation between true value and reconstructed value of these checkpoints
%%% OutData0 is the initial cloud removal result of original cloudy image A0
%%% A_Mask0 is the original cloud mask; A_Mask1 is the updated cloud mask
OutputData2=OutData0;
A_Mask1(A_Mask1>0)=1;
A_Mask11=A_Mask1;
A_Mask11(A_Mask0>0)=0;%A_Mask1
[ch,cw]=find(A_Mask0>0);%
[bh,bw]=find(A_Mask11>0);%
OutData0=rgb2gray(OutData0);
for i=1:length(ch)
    patchic=double(OutData0(ch(i)-1:ch(i)+1,cw(i)-1:cw(i)+1));%
    mcossim=-1;mj=-1;
    for j=1:length(bh)
        patchic_c=double(OutData0(bh(j)-1:bh(j)+1,bw(j)-1:bw(j)+1));%
         cossim=sum(sum(patchic.*patchic_c))/(sqrt(sum(sum(patchic.*patchic)))*sqrt(sum(sum(patchic_c.*patchic_c))));
   %     cossim=1.0/(sum(sum(abs(patchic_c-patchic)))+0.000001);
        if(cossim>mcossim)
            mcossim=cossim;
            mj=j;
        end
    end
    if (abs(OutData0(ch(i),cw(i))-OutData0(bh(mj),bw(mj)))<12 && mcossim>0.7)% 0.05
    %OutputData2(ch(i),cw(i),:)=A0(bh(mj),bw(mj),:);
     OutputData2(ch(i),cw(i),:)=A0(bh(mj),bw(mj),:)*weights(1)+OutputData2(ch(i),cw(i),:)*weights(2);
    end
end
% % OutputData2=double(0.8*double(OutputData2)+0.2*double(OutData0));
% if(isa(OutData0,'uint8'))
%     OutputData2=uint8(OutputData2);
% end
end