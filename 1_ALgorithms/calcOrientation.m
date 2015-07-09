function [ orientacion ] = calcOrientation(tri )
%calcOrientation calcula la orientacion de un triangulo de entrada

    orden=sort(tri(:,1));%Ordena segun posiciones de x

    indice=find(tri(:,1)==orden(2));
    indice2=find(tri(:,1)==orden(1));


    if(tri(indice,2)>tri(indice2,2) )%compara la ubicacion del 2do punto

       orientacion=1;


    else
       orientacion=2;


    end



end

