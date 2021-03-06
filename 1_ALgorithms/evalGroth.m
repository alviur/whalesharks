%==========================================================================
% Evaluacion del algoritmo de identificacion de Groth
%
%
%==========================================================================



    contador=1;
    individuo0='ind';

    PathIMG = '/home/lex/2_SISTEMIC/3_MRF/1_Lagartos/2_Dataset/GT/ind1/DSCN1551.JPG';%Imagen a comparar
    Igt=imread(PathIMG);

    Igt=imresize(Igt,[300 300]);   

    centrosOriginal=extractPattern(Igt);
    
    comparacion=zeros(20,2);



    contador=1;
   individuo0='ind';
   
    for j=1:20
        
      %  /usr/local/MATLAB/MATLAB_Production_Server/R2015a/bin
            numeroToString = int2str(j);
            individuo=horzcat(individuo0,numeroToString,'/');/usr/local/MATLAB/MATLAB_Production_Server/R2015a/bin

            % Paths a archivos
       
            PathResul2='/home/lex/2_SISTEMIC/3_MRF/1_Lagartos/2_Dataset/GT_rotado/';%MS 

            extension='*.JPG';

            % Añade path del folder           
            PathResul3=horzcat(PathResul2,individuo);

            %Extrae cantidad de archivos en folder
            archives=dir(horzcat(PathResul3,extension));
            %horzcat(PathIMG2,extension)

            for i=1:1
                % Navega por cada archivo

                
                PathResul=horzcat(PathResul3,archives(i).name);
                
                
                % Carga imagenes                

                Iraw=imread(PathResul);
                Iraw=imresize(Iraw,[300 300]);
                
                 
                centrosPra=extractPattern(Iraw);                

                
                %voting=grothAlgorithm(centrosOriginal,centrosPra,0.01);
                voting=spotMatchingAlgorithm(centrosOriginal,centrosPra,0.001,(10*2*pi)/360,0.85);
                
                comparacion(j,1)=voting;
                
                
            end
        
    end