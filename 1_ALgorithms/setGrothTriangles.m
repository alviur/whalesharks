%==========================================================================
%  SetTriangles
%  Authors: Alexander Gomez - German Diez
% -----------------------------------------------------------------------
%  Based on:    
%  Groth, Edward J. "A pattern-matching algorithm for two-dimensional
%  coordinate lists." The astronomical journal 91 (1986): 1244-1248.
% -----------------------------------------------------------------------
%   function triCorrected=setGrothTriangles(tri,A)
%
%   A : n x 2 matrix with n points wirh x and y position 
%   tri : Delaunay triangulation of A
%   triCorrected: triangles with Groth criterium
%==========================================================================

function triCorrected=setGrothTriangles(tri1,A)

    for i=1:size(tri1,1)
        
        % Calculo longitudes cada lado
        l1=((A(tri1(i,2),1)-A(tri1(i,1),1))^2 + (A(tri1(i,2),2)-A(tri1(i,1),2))^2)^(0.5);
        l2=((A(tri1(i,3),1)-A(tri1(i,2),1))^2 + (A(tri1(i,3),2)-A(tri1(i,2),2))^2)^(0.5);
        l3=((A(tri1(i,1),1)-A(tri1(i,3),1))^2 + (A(tri1(i,1),2)-A(tri1(i,3),2))^2)^(0.5);
        
        % Guarda puntos para intercambiar        
        r1=tri1(i,1);
        r2=tri1(i,2);
        r3=tri1(i,3);
        
        
        if(l1>l2 && l1>l3)
            
            tri1(i,3)=r2;
            tri1(i,2)=r3;
            
            if(l2<l3)
                
                tri1(i,1)=r2;
                tri1(i,3)=r1;
                tri1(i,2)=r3;               
                
            end
            
        elseif(l1<l2 && l1<l3)
            
            if(l2>l3)
                tri1(i,1)=r2;
                tri1(i,2)=r1;                
                
            end
            
        elseif(l1>l2 && l1<l3)
            
                tri1(i,1)=r3;
                tri1(i,3)=r1; 
                
                
         elseif(l1<l2 && l1>l3)

                tri1(i,1)=r2;
                tri1(i,2)=r3;     
                tri1(i,3)=r1; 
            
        end
        
        
    
    end


        
    triCorrected=tri1;


end