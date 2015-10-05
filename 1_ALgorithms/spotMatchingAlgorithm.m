%==========================================================================
%  Spot Matching algorithm
%  Author: Alexander Gomez 
% -----------------------------------------------------------------------
%  Based on:    
%  Arzoumanian, Z., J. Holmberg, and B. Norman. "An astronomical pattern‐
%  matching algorithm for computer‐aided identification of whale sharks
%  Rhincodon typus." Journal of Applied Ecology 42, no. 6 (2005): 999-1011.
% -----------------------------------------------------------------------
%   function [v1,mlog]=spotMatchingAlgorithm(A,B,epsilon)
%
%   A and B are n x 2 matrices with n points wirh x and y position 
%   epsilon: Tolerance factor Groth algorithm
%==========================================================================


function [v1]=spotMatchingAlgorithm(A,B,epsilon,thetaMax,sMax)


%     PathIMG = '/home/lex/2_SISTEMIC/3_MRF/1_Lagartos/2_Dataset/GT/ind1/DSCN1551.JPG';%Imagen a comparar
%     Igt=imread(PathIMG);
%     Igt=imresize(Igt,[300 300]); 
%     A=extractPattern(Igt);
% 
%     PathIMG2 = '/home/lex/2_SISTEMIC/3_MRF/1_Lagartos/2_Dataset/GT/ind1/DSCN1551.JPG';%Imagen a comparar
%     Igt2=imread(PathIMG2);
%     Igt2=imresize(Igt2,[300 300]); 
%     B=extractPattern(Igt2);
%     thetaMax=1;
%     sMax=0.85;
%     epsilon=0.01;

    %% Puntos de entrada

   % A=[1,7;9,2;3,2;10,11;8,13;12,19;3,2];
    %B=[1,9;6,2;3,1;12,11;4,13;13,19;6,2];
   % epsilon=0.01;%parametro de tolerancia para incertidumbre en medidas
   % thetaMax:rotacion relativa entre patrones de puntos

    % Escalamiento


for p=1:1   
    
    if(p>1)
     
    [row,col]=find( pointTable1(:,3));   
    A=A(row,:) 
    [row,col]=find( pointTable2(:,3));
    B=B(row,:) ; 
        
    end

        %% Calculo de triangulos

        % Delaunay triangulation

        tri1 = delaunay(A(:,1),A(:,2));
        tri2 = delaunay(B(:,1),B(:,2));


        % Graficas
        %subplot(1,2,1)
        %triplot(tri1,A(:,1),A(:,2))% grafica triangulos
        %subplot(1,2,2)
        %triplot(tri2,B(:,1),B(:,2))% grafica triangulos

        % Organizacion triangulos segun Groth
        tri1=setGrothTriangles(tri1,A);
        tri2=setGrothTriangles(tri2,B);

        % Normalizacion de coordenadas
        A = (A - min(min(A)))/(max(max(A)) - min(min(A)));
        B = (B - min(min(B)))/(max(max(B)) - min(min(B)));

        % Filtro de triangulos 


        %% Calculo de descriptores

        if(size(A,1)>0 && size(B,1)>0)
           % size(A,1)
           % size(B,1)
            [r1,r2,r3,R,C,F,tr,tc,M,orientacion,theta,s] = spotMatchingDescriptors(A,tri1,epsilon);
            [r1b,r2b,r3b,Rb,Cb,Fb,trb,tcb,Mb,orientacionb,thetab,sb] = spotMatchingDescriptors(B,tri2,epsilon);       



        %% Matching

            % tablas

            triangleTab=ones(size(tri1,1),size(tri2,1),4)*300000;%Tabla con indices de triangulos iguales

            pointTable1=[A zeros(size(A,1),1)]; %tabla puntos 1
            pointTable2=[B zeros(size(B,1),1)]; % tabla puntos 2
            %magTable=zeros(1,1);

            % Contadores
             contFalseMat=0;
             contTrueMat=0;
             Match1Ant=10000;
             Match2Ant=10000;
             magTable=zeros(1,1);
             flag=0;


            % Llenado de la tabla de matching de triangulos

            for t=1:20


                for i=1:size(tri1,1)


                    for j=1:size(tri2,1)
                        

                        cuadratureDiff=(((R(i)-Rb(j))^2)/(tr(i)+trb(j)) + ((C(i)-Cb(j))^2/(tc(i)+tcb(j))) + ...
                                        (theta(i)-thetab(j))^2/(thetaMax^2));
                        if((((R(i)-Rb(j))^2)<tr(i) + trb(j))  && ((C(i)-Cb(j))^2<tc(i) + tcb(j) && ...
                                (R(i)-Rb(j))^2)<min(triangleTab(i,:,2)) && (C(i)-Cb(j))^2< min(triangleTab(i,:,3)) && ...
                                theta(i)-thetab(j)<thetaMax )
                            


                            [a,c]=min(triangleTab(i,:,2));
                            [a,c2]=min(triangleTab(:,j,3));


                             triangleTab(i,c,2)=300000;
                             triangleTab(i,c,3)=300000;
                             triangleTab(i,c,1)=300000;

                             triangleTab(i,j,2)=(R(i)-Rb(j))^2;
                             triangleTab(i,j,3)=(C(i)-Cb(j))^2;


                            if( (triangleTab(i,j,2)<tr(i) + trb(j))  &&  triangleTab(i,j,3)<tc(i) + tcb(j) && ...
                                   cuadratureDiff<min(triangleTab(i,:,4)) )

                               
                               triangleTab(i,j,4)=cuadratureDiff;
                               
                                % Calculo de orientacion de los triangulos
                                orientTri1= ((A(tri1(i,2),1))- (A(tri1(i,1),1)))*((A(tri1(i,3),2))- (A(tri1(i,1),2)))-...
                                            ((A(tri1(i,2),2))- (A(tri1(i,1),2)))*((A(tri1(i,3),1))- (A(tri1(i,1),1)));

                                orientTri2= ((B(tri2(j,2),1))- (B(tri2(j,1),1)))*((B(tri2(j,3),2))- (B(tri2(j,1),2)))-...
                                            ((B(tri2(j,2),2))- (B(tri2(j,1),2)))*((B(tri2(j,3),1))- (B(tri2(j,1),1)));

                                if((orientTri1*orientTri2>=0) && R(i)<8 && Rb(j)<8 && sb(j)<sMax && s(i)<sMax) 

                                    if(t>1)

                                        if(M(i)-Mb(j)>=media-abs(desviacion) && M(i)-Mb(j)<=media+abs(desviacion))

                                            triangleTab(i,j,2)=(R(i)-Rb(j))^2;
                                            triangleTab(i,j,3)=(C(i)-Cb(j))^2;

                                            triangleTab(i,j,1)=1; % triangulos cumplen condiciones de matching
                                            magTable=[magTable  M(i)-Mb(j)]; % Calculo factor de magnificacion


                                            pointTable1(tri1(c2,1),3)= pointTable1(tri1(i,1),3)-1;
                                            pointTable1(tri1(c2,2),3)= pointTable1(tri1(i,2),3)-1;
                                            pointTable1(tri1(c2,3),3)= pointTable1(tri1(i,3),3)-1;

                                            pointTable1(tri1(i,1),3)= pointTable1(tri1(i,1),3)+1;
                                            pointTable1(tri1(i,2),3)= pointTable1(tri1(i,2),3)+1;
                                            pointTable1(tri1(i,3),3)= pointTable1(tri1(i,3),3)+1;


                                            pointTable2(tri2(c,1),3)= pointTable2(tri2(j,1),3)-1;
                                            pointTable2(tri2(c,2),3)= pointTable2(tri2(j,2),3)-1;
                                            pointTable2(tri2(c,3),3)= pointTable2(tri2(j,3),3)-1;

                                            pointTable2(tri2(j,1),3)= pointTable2(tri2(j,1),3)+1;
                                            pointTable2(tri2(j,2),3)= pointTable2(tri2(j,2),3)+1;
                                            pointTable2(tri2(j,3),3)= pointTable2(tri2(j,3),3)+1;

                                            contTrueMat=contTrueMat+1;

                                        else
                                           contFalseMat=contFalseMat +1; 

                                        end

                                    else



                                        triangleTab(i,j,2)=(R(i)-Rb(j))^2;
                                        triangleTab(i,j,3)=(C(i)-Cb(j))^2;

                                        triangleTab(i,j,1)=1; % triangulos cumplen condiciones de matching

                                        if(flag==0)
                                            magTable=[M(i)-Mb(j)]; % Calculo factor de magnificacion
                                            flag=1;

                                        else
                                            magTable=[magTable,M(i)-Mb(j)]; % Calculo factor de magnificacion

                                        end


                                        pointTable1(tri1(c2,1),3)= pointTable1(tri1(i,1),3)-1;
                                        pointTable1(tri1(c2,2),3)= pointTable1(tri1(i,2),3)-1;
                                        pointTable1(tri1(c2,3),3)= pointTable1(tri1(i,3),3)-1;


                                        pointTable1(tri1(i,1),3)= pointTable1(tri1(i,1),3)+1;
                                        pointTable1(tri1(i,2),3)= pointTable1(tri1(i,2),3)+1;
                                        pointTable1(tri1(i,3),3)= pointTable1(tri1(i,3),3)+1;

                                        pointTable2(tri2(c,1),3)= pointTable2(tri2(j,1),3)-1;
                                        pointTable2(tri2(c,2),3)= pointTable2(tri2(j,2),3)-1;
                                        pointTable2(tri2(c,3),3)= pointTable2(tri2(j,3),3)-1;

                                        pointTable2(tri2(j,1),3)= pointTable2(tri2(j,1),3)+1;
                                        pointTable2(tri2(j,2),3)= pointTable2(tri2(j,2),3)+1;
                                        pointTable2(tri2(j,3),3)= pointTable2(tri2(j,3),3)+1;

                                        contTrueMat=contTrueMat+1;

                                    end


                                else
                                    contFalseMat=contFalseMat +1;

                                end



                            end


                        end

             end

             Match1Ant=10000;
             Match2Ant=10000;



         end

         %calculo media y desviacion filtro magnificacion

         media=mean(magTable);
         desviacion=std(magTable);


         aux=triangleTab;     
         aux(aux==300000)=0;
         contTrueMat=sum(sum(aux(:,:,1)));
         contFalseMat=sum(sum(not(aux(:,:,1))));
         mt=abs(contTrueMat-contFalseMat);
         mf=contTrueMat+contFalseMat-mt;

         if(contFalseMat>contTrueMat)

                 desviacion=desviacion;            


         elseif (mf<0.5*mt)

                desviacion=3*desviacion;

         else

                   desviacion=2*des4viacion;                  

         end

         if(t<20)
             % Reinicio de tablas y variables
             %Tablas
            triangleTab=ones(size(tri1,1),size(tri2,1),4)*300000;%Tabla con indices de triangulos iguales
            pointTable1=[A zeros(size(A,1),1)]; %tabla puntos 1
            pointTable2=[B zeros(size(B,1),1)]; % tabla puntos 2
            magTable=zeros(1,1);
            % Contadores
            contFalseMat=0;
            contTrueMat=0;
            Match1Ant=10000;
            Match2Ant=10000;
            flag=0;
         end



            end 


            triangleTab(triangleTab==300000)=0;
            %% Voting

            v1= ( sum((pointTable2(:,3)>1))*(sum(sum(triangleTab(:,:,1))))/size(tri2,1)) ;     
            



        else
            v1=0;

        end
end 
    
end