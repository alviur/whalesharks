%==========================================================================
% function - extractPattern
%
%
%==========================================================================




function pattern=extractPattern(binImage)


    % Transformacion a escala de grises
    if(size(binImage,3)>1)
        binImage=rgb2gray(binImage);
    end

    % Correccion formato JPG
    A=find(binImage>=100);
    B=find(binImage<100);
    binImage(A)=255;
    binImage(B)=0;
    
    
    binImage=logical(binImage);% Binariza la imagen

    
    % Extraccion centroides regiones
    centrosGT  = regionprops(binImage,'centroid');
    centrosOriginal=zeros(size(centrosGT,1),2);
    
    for i=1:size(centrosGT,1)
       
        pattern(i,:)=centrosGT(i).Centroid;
        
    end
    
    pattern=round(pattern);% Redondeo



end