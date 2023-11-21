clc, clear
str_tif='.tif';
I_str='I0'; 
J_str='J0';
Imask_str='I0_mask';
Jmask_str='J0_mask';

%% =============================== data 1: 1_TrueColorComposition (Reference image J0 is Cloud-free)=============================== 
%%% the Whole Image Cloud Removal by Curvature
str_InputPath='..\data\1_TrueColorComposition';
I=imread([str_InputPath,'\',I_str,'_Clouds',str_tif]);
J=imread([str_InputPath,'\',J_str,str_tif]);
I_mask=imread([str_InputPath,'\',Imask_str,str_tif]);
BufferSize=6; %%the size of buffer area. you can change the parameters to improve the effect
weights=[0.95 0.1]; %%the weights of correcting.you can change the parameters to improve the effect
[OutputData_II,~]=Whole_Image_Cloud_Removal_by_Curvature(I,I_mask,J,I_mask*0, BufferSize,weights);
%%% output the result to specified path
imwrite((OutputData_II),[str_InputPath,'\',I_str,'-Result_of_Curvature',str_tif]);

%% =============================== data 2: 2_FalseColorComposition (both I0 and J0 are Cloudy)=============================== 
%%% the Whole Image Cloud Removal by Curvature
str_InputPath='..\data\2_FalseColorComposition';
I=imread([str_InputPath,'\',I_str,'_Clouds',str_tif]);
J=imread([str_InputPath,'\',J_str,'_Clouds',str_tif]);
I_mask=imread([str_InputPath,'\',Imask_str,str_tif]);
J_mask=imread([str_InputPath,'\',Jmask_str,str_tif]);
BufferSize=6; %%the size of buffer area. you can change the parameters to improve the effect
weights=[0.99 0.03]; %%the weights of correcting.you can change the parameters to improve the effect
[OutputData_II,OutputData_JJ]=Whole_Image_Cloud_Removal_by_Curvature(I,I_mask,J,J_mask, BufferSize,weights);
%%% output the result to specified path
imwrite((OutputData_II),[str_InputPath,'\',I_str,'-Result_of_Curvature',str_tif]);
imwrite((OutputData_JJ),[str_InputPath,'\',J_str,'-Result_of_Curvature',str_tif]);
