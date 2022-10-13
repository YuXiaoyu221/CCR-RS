function OutputData2= FindBestBoundry(inputData00,DataMask00, inputData10,BoundryType)
%%%��������Ĥ�ļ�DataMask0���б߽�������ҵ�����Ӱ��ɢ�Ȳ����С���µı߽��
%%%�������inputData0��DataMask0��inputData1����Ϊuint8���ͣ�����ҲΪuint8��
%%%��������ǵ����߽���mask��unit8��
%%%BoundryType����1��2��3��1��ʾ��������ֵ����������ѱ߽�㣬2��ʾ����ԭӰ���ݶ�ֵ����ѱ߽�㣬3��ʾ����ɢ��ֵ����������ѱ߽��
BufferSize=3;
%% ��ǰ���䣬��ֹ����Ӱ��߽�����
[nh,nw,b]=size(inputData00);
inputData0=zeros(nh+20,nw+20,b)-1;inputData0(11:end-10,11:end-10,:)=inputData00;
inputData1=zeros(nh+20,nw+20,b)-1;inputData1(11:end-10,11:end-10,:)=inputData10;
DataMask0=zeros(nh+20,nw+20,1);DataMask0(11:end-10,11:end-10,:)=DataMask00;
%% ��������Ӱ��ɢ�Ȳ���
dif_ImgDivergence=double(inputData0);
if(BoundryType==3)%��������Ӱ��ɢ�Ȳ���
    inputData0=imfilter(double(inputData0), [0,1,0;1,-4,1;0,1,0]);
    inputData1=imfilter(double(inputData1), [0,1,0;1,-4,1;0,1,0]);
    dif_ImgDivergence=abs(inputData0-inputData1);
    dif_ImgDivergence=sum(dif_ImgDivergence,3);
    DataMask0=imfilter(DataMask0,[0,1,0;1,1,1;0,1,0]);%���������mask��������1���߽��ɢ�Ȳ��ɿ�
elseif (BoundryType==2)%ԭӰ��ɢ��
    inputData0=abs(imfilter(double(inputData0), [0,1,0;1,-4,1;0,1,0]));
    dif_ImgDivergence=sum(inputData0,3);
    DataMask0=imfilter(DataMask0,[0,1,0;1,1,1;0,1,0]);%���������mask��������1���߽���ݶȲ��ɿ�
elseif(BoundryType==1)%��������Ӱ������ֵ����
    dif_ImgDivergence=abs(double(inputData0)-double(inputData1));
    dif_ImgDivergence=sum(dif_ImgDivergence,3);
else
    disp('��ȷ��BoundryTypeΪ1��2��3��')
end
% abc=dif_ImgDivergence;
% figure,imshow(DataMask0)
%% ��ȡԭʼ�߽��
OriBoundryMask=imfilter(DataMask0,[0,1,0;1,1,1;0,1,0]);%���������mask��������1
OriBoundryMask(OriBoundryMask>0)=1;
OriBoundryMask(DataMask0>0)=0;
%%  ������Ĥ�ļ�������
DataMask=imfilter(DataMask0,ones(BufferSize*2+1,BufferSize*2+1));%���������mask��������1
DataMask(DataMask>0)=1;
DataMask(DataMask0>0)=0;
dif_ImgDivergence(DataMask==0)=-1;
%% ȷ��mask�ļ��������ڵı߽��
[Array_H,Array_W]=find(DataMask>0);%���������������ص�
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
index_Boundry=index_Boundry(1:num_Boundry,:);%�������߽��
%% ��mask�ļ��Ļ���������ɢ�Ȳ����С�ı߽�㣬����ȷ���߽�
for iter=1:BufferSize-1
    %��ԭ�������ı߽����£������º�߽���ɢ��ֵС��֮ǰ�ģ�����и��£�������
    temp_Divergence=ones(5,5);
    index_NewBoundry=zeros(length(Array_W),2);
    num_NewBoundry=0;
    for i=1:num_Boundry
        %�Ե�ǰ�߽��Ϊ���ĵ�������
        temp_Divergence(2:4,2:4)=dif_ImgDivergence(index_Boundry(i,1)-1:index_Boundry(i,1)+1,index_Boundry(i,2)-1:index_Boundry(i,2)+1);
        s0=temp_Divergence(3,3);
        %���������ڶ����ı߽����и��¡����ڸ��º�ԭ�߽����Ϊ�߽�㣬���������ı߽���⣬ֻ�������������ı߽�㲻���ܼ��������ı߽��
        %�ʽ����ĵ����Ϊ-1���ж�������������û�����ӵı߽�㣬����ɢ�Ⱥͼ��ɡ�
        temp_Divergence(3,3)=-1;
        s1=0;n1=0;
        jj=[1,0;-1,0;0,1;0,-1];
        if(OriBoundryMask(index_Boundry(i,1),index_Boundry(i,2))==1)%����õ㱾������ԭ�߽�ֵ
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
        
        if(s1<=s0) %��������ɢ��ֵ<ԭɢ��ֵ����Ա߽���и���
            dif_ImgDivergence(index_Boundry(i,1),index_Boundry(i,2))=-1;
            DataMask(index_Boundry(i,1),index_Boundry(i,2))=0;
            if(n1==0) %%�ж������ĵ�Ϊ�Ǳ�Ҫ�㣬ֱ�ӽ���ɾ��
                continue;
            end
            index_NewBoundry(num_NewBoundry+1:num_NewBoundry+n1,:)=temp_Boundry;
            num_NewBoundry=num_NewBoundry+n1;
        else %��������ɢ��ֵ>ԭɢ��ֵ���򲻽��и���
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
%%�±߽�ֻ����ԭ�߽������ӣ����ܼ�С����֤һ��
DataMask0(DataMask0~=0)=255;
DataMask(DataMask~=0)=255;
aa=DataMask-DataMask0;
% bb=zeros(size(aa));bb(aa<0)=255;figure,imshow(bb)
if(length(find(aa<0))>0)
    disp('�±߽�ֻ����ԭ�߽������ӣ����ܼ�С������');
end
end