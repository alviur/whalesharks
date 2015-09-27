%==========================================================================
%  Groth Descriptors
%  Authors: Alexander Gomez - German Diez
% -----------------------------------------------------------------------
%  Based on:    
%  Groth, Edward J. "A pattern-matching algorithm for two-dimensional
%  coordinate lists." The astronomical journal 91 (1986): 1244-1248.
% -----------------------------------------------------------------------
%   [r1,r2,r3,R,C,F,tr,tc,M,orientacion] = grothDescriptors(A,tri,epsilon )
%
%   A : n x 2 matrix with n points wirh x and y position 
%   tri : Delaunay triangulation of A
%   epsilon: Tolerance factor Groth algorithm
%==========================================================================

function [r1,r2,r3,R,C,F,tr,tc,M,orientacion] = grothDescriptors(A,tri,epsilon )
%grothDescriptors Calcula las caracteristicas para algoritmo de Groth


r2=((A(tri(:,2),1)-A(tri(:,1),1)).^2) + ((A(tri(:,2),2)-A(tri(:,1),2)).^2);
r3=((A(tri(:,3),1)-A(tri(:,1),1)).^2) + ((A(tri(:,3),2)-A(tri(:,1),2)).^2);

% Ratio
R=r3./r2;

% Coseno del angulo en vertex 1
C=((A(tri(:,3),1)-A(tri(:,1),1)).*(A(tri(:,2),1)-A(tri(:,1),1))+(A(tri(:,2),2)...
    -A(tri(:,1),2)).*(A(tri(:,3),2)-A(tri(:,1),2)))./(times(r3,r2));

F=(epsilon^2)*((1./r3.^2)-(C./(r3.*r2))+(1./r2.^2));

% Tolerancias

tr=(2.*R.^2.*F).^(0.5);
tc=(2.*(1-C.^2).*F + 3.*(C.^2).*(F.^2)).^(0.5);

% Factor de magnificacion

r1=((A(tri(:,2),1)-A(tri(:,3),1)).^2) + ((A(tri(:,2),2)-A(tri(:,3),2)).^2);
M=log10( r2+r3+r1);


% Calculo de orientaicion

orientacion=zeros(size(tri,1));
for i=1:size(tri,1)     
    
    orientacion(i) = calcOrientation([A(tri(i,1),:);A(tri(i,2),:);A(tri(i,3),:)]);
       
end

end

