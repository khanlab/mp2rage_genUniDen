
function mp2rage_genUniDen(MP2RAGE_filenameUNI, MP2RAGE_filenameINV1, MP2RAGE_filenameINV2,MP2RAGE_uniden_output_filename,varargin)
% Usage: mp2rage_genUniDen(MP2RAGE_filenameUNI, MP2RAGE_filenameINV1,
% MP2RAGE_filenameINV2,MP2RAGE_uniden_output_filename, multiplyingFactor
% (default=6, increase up to 10 for more noise suppression)

% adapted from RobustCombination function from Jose Marques, https://github.com/JosePMarques/MP2RAGE-related-scripts
% this function shows one possible implementation of the methods suggested
% in http://journals.plos.org/plosone/article?id=10.1371/journal.pone.0099676

if (nargin>4)
    chosenFactor=varargin{1};
else
    chosenFactor=6;
end


disp(sprintf('using multiplying factor of %g',chosenFactor));

%% defines relevant functions

MP2RAGErobustfunc  =@(INV1,INV2,beta)(conj(INV1).*INV2-beta)./(INV1.^2+INV2.^2+2*beta);

rootsquares_pos  =@(a,b,c)(-b+sqrt(b.^2 -4 *a.*c))./(2*a);
rootsquares_neg  =@(a,b,c)(-b-sqrt(b.^2 -4 *a.*c))./(2*a);


%% load Data

MP2RAGEimg=load_untouch_nii(MP2RAGE_filenameUNI);
INV1img=load_untouch_nii(MP2RAGE_filenameINV1);
INV2img=load_untouch_nii(MP2RAGE_filenameINV2);

if and(min(MP2RAGEimg.img(:))>=0,max(MP2RAGEimg.img(:))>=0.51)
    % converts MP2RAGE to -0.5 to 0.5 scale - assumes that it is getting only
    % positive values
    MP2RAGEimg.img=(double(MP2RAGEimg.img)- max(double(MP2RAGEimg.img(:)))/2)./max(double(MP2RAGEimg.img(:)));
    integerformat=1;
    
else
    integerformat=0;
end


%% computes correct INV1 dataset  
INV2img.img=double(INV2img.img);

%gives the correct polarity to INV1;
INV1img.img=sign(MP2RAGEimg.img).*double(INV1img.img);

%
% because the MP2RAGE INV1 and INV2 is a summ of squares data, while the
% MP2RAGEimg is a phase sensitive coil combination.. some more maths has to
% be performed to get a better INV1 estimate which here is done by assuming
% both INV2 is closer to a real phase sensitive combination


INV1pos=rootsquares_pos(-MP2RAGEimg.img,INV2img.img,-INV2img.img.^2.*MP2RAGEimg.img);
INV1neg=rootsquares_neg(-MP2RAGEimg.img,INV2img.img,-INV2img.img.^2.*MP2RAGEimg.img);


INV1final=INV1img.img;
INV1final(abs(INV1img.img-INV1pos)> abs(INV1img.img-INV1neg))=INV1neg(abs(INV1img.img-INV1pos)>abs(INV1img.img-INV1neg));
INV1final(abs(INV1img.img-INV1pos)<=abs(INV1img.img-INV1neg))=INV1pos(abs(INV1img.img-INV1pos)<=abs(INV1img.img-INV1neg));



%% visualizing the data

pos=round(3/5*size(INV1final));
visualize=0;
if visualize
    figureJ(200)
    subplot(411)
    Orthoview(INV1pos,pos,[-200 200])
    title('positive root')
    
    subplot(412)
    Orthoview(INV1neg,pos,[-200 200])
    title('negative root')
    
    subplot(413)
    Orthoview(INV1img.img,pos,[-200 200])
    title('Phase Corrected Sum of Squares  root')
    
    subplot(414)
    Orthoview(INV1final,pos,[-200 200])
    title('INV1 final')
    
end

clear INV1img
clear INV1pos
clear INV1neg

%% lambda calculation

% usually the multiplicative factor shouldn't be greater then 10, but that
% is not the ase when the image is bias field corrected, in which case the
% noise estimated at the edge of the imagemight not be such a good measure

   multiplyingFactor=chosenFactor;
    
    noiselevel=multiplyingFactor*mean(mean(mean(INV2img.img(1:end,end-10:end,end-10:end))));
    
    
    % MP2RAGEimgRobustScanner=MP2RAGErobustfunc(INV1img.img,INV2img.img,noiselevel.^2);
    MP2RAGEimgRobustPhaseSensitive=MP2RAGErobustfunc(INV1final,INV2img.img,noiselevel.^2);
    
     % Robust Image view
    range =[-0.5 0.40];
    if visualize
        
    subplot(211)
    Orthoview(MP2RAGEimg.img,pos,range),title('MP2RAGE UNI_Image')
    
    % subplot(312)
    % Orthoview(MP2RAGEimgRobustScanner,pos,range),title('MP2RAGE Robust Scanner')
    
    subplot(212)
    Orthoview(MP2RAGEimgRobustPhaseSensitive,pos,range),title('MP2RAGE Robust')
    ylabel(['noise level = ',num2str(multiplyingFactor)])
    
    end
    
    
    
%filenameOUT=strrep(MP2RAGE.filenameUNI,'_acq-UNI_',sprintf('_acq-UNIDEN%d_',multiplyingFactor));

        if integerformat==0
            MP2RAGEimg.img=MP2RAGEimgRobustPhaseSensitive;
            save_untouch_nii(MP2RAGEimg,MP2RAGE_uniden_output_filename);
        else
            MP2RAGEimg.img=round(4095*(MP2RAGEimgRobustPhaseSensitive+0.5));
            save_untouch_nii(MP2RAGEimg,MP2RAGE_uniden_output_filename);
            
        end
        

            




end

