%==========================================================================
%  Groth algorithm
%  Authors: Alexander Gomez - German Diez
% -----------------------------------------------------------------------
%  Based on:    
%  Groth, Edward J. "A pattern-matching algorithm for two-dimensional
%  coordinate lists." The astronomical journal 91 (1986): 1244-1248.
% -----------------------------------------------------------------------
%   function [v1,mlog]=grothAlgorithm(A,B,epsilon)
%
%   A and B are n x 2 matrices with n points wirh x and y position 
%   epsilon: Tolerance factor Groth algorithm
%==========================================================================


function [v1,mlog]=grothAlgorithm(A,B,epsilon)

    %% Puntos de entrada

   % A=[1,7;9,2;3,2;10,11;8,13;12,19;3,2];
    %B=[1,9;6,2;3,1;12,11;4,13;13,19;6,2];
    %epsilon=0.01;%parametro de tolerancia para incertidumbre en medidas

    % Escalamiento


    %% Calculo de triangulos

    % Delaunay triangulation

    tri1 = delaunay(A(:,1),A(:,2));
    tri2 = delaunay(B(:,1),B(:,2));

    %triplot(tri1,A(:,1),A(:,2))% grafica triangulos
    %figure
    %triplot(tri2,B(:,1),B(:,2))% grafica triangulos


    % Filtro de triangulos 


    %% Calculo de descriptores
    
    if(size(A,1)>0 && size(B,1)>0)
       % size(A,1)
       % size(B,1)
        [r1,r2,r3,R,C,F,tr,tc,M,orientacion] = grothDescriptors(A,tri1,epsilon);
        [r1b,r2b,r3b,Rb,Cb,Fb,trb,tcb,Mb,orientacionb] = grothDescriptors(B,tri2,epsilon);
        
  
    


        %% Matching

        % tablas

        triangleTab=zeros(size(tri1,1),size(tri2,1),3);%Tabla con indices de triangulos iguales

        pointTable1=[A zeros(size(A,1),1)]; %tabla puntos 1
        pointTable2=[B zeros(size(B,1),1)]; % tabla puntos 2
        magTable=zeros(1,1);
        
        % Contadores
         contFalseMat=0;
         contTrueMat=0;


        % Llenado de la tabla de matching de trigngulos

        for i=1:size(tri1,1)


            for j=1:size(tri2,1)

                if(((R(i)-Rb(j))^2<tr(i)^2 + trb(j)^2)  && ((C(i)-Cb(j))^2<tc(i)^2 + tcb(j)^2))


                    if( (triangleTab(i,j,2)<tr(i)^2 + trb(j)^2)  &&  triangleTab(i,j,3)<tc(i)^2 + tcb(j)^2)
                        
                        % Calculo de orientacion de los triangulos
                        orientTri1= ((A(tri1(i,2),1))- (A(tri1(i,1),1)))*((A(tri1(i,3),2))- (A(tri1(i,1),2)))-...
                                    ((A(tri1(i,2),2))- (A(tri1(i,1),2)))*((A(tri1(i,3),1))- (A(tri1(i,1),1)));
                        
                        orientTri2= ((B(tri2(i,2),1))- (B(tri2(i,1),1)))*((B(tri2(i,3),2))- (B(tri2(i,1),2)))-...
                                    ((B(tri2(i,2),2))- (B(tri2(i,1),2)))*((B(tri2(i,3),1))- (B(tri2(i,1),1)));
                                
                        if((orientTri1*orientTri2>=0))        

                            triangleTab(i,j,2)=tr(i)^2 + trb(j)^2;
                            triangleTab(i,j,3)=tc(i)^2 + tcb(j)^2;

                            triangleTab(i,j,1)=1; % triangulos cumplen condiciones de matching
                            magTable=[magTable  M(i)-Mb(j)]; % Calculo factor de magnificacion
                            pointTable1(tri1(i,1),3)= pointTable1(tri1(i,1),3)+1;
                            pointTable1(tri1(i,2),3)= pointTable1(tri1(i,1),3)+1;
                            pointTable1(tri1(i,3),3)= pointTable1(tri1(i,1),3)+1;

                            pointTable2(tri2(i,1),3)= pointTable2(tri2(i,1),3)+1;
                            pointTable2(tri2(i,2),3)= pointTable2(tri2(i,1),3)+1;
                            pointTable2(tri2(i,3),3)= pointTable2(tri2(i,1),3)+1;
                            
                            contTrueMat=contTrueMat+1;
                            
                            
                        else
                            contFalseMat=contFalseMat +1;
                        
                        end
                        
                        

                    end


                end

            end


        end

        % Voting

        v1=(sum((pointTable1(:,3)>1)))  + (sum((pointTable2(:,3)>1))) ;
        
        size(M);
         size(Mb);
        %mag=M-Mb;
       

        mt=abs(contTrueMat-contFalseMat);
        mf=contTrueMat+contFalseMat-mt;
        
         if(mf>mt)
             desviation=1
             
         end    
        
         if (0.1*mt>mf)
                 
            desviation=3;
            
             else
             
               desviation=2;                  
        
         end
        
        mlog=sum(magTable>-desviation & magTable<desviation);
        
        
    else
        v1=0;
        mlog=0;
        
    end
    
    
end
