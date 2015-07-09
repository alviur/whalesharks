%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Algoritmo de Groth
%  Autores: Alexander Gomez - German Diez
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%clear all, close all, clc


%% Puntos de entrada

A=[1,7;9,2;3,2;10,11;8,13;12,19;3,2];
epsilon=2;%parametro de tolerancia para incertidumbre en medidas

% Escalamiento


%% Calculo de triangulos

% Delaunay triangulation

tri = delaunay(A(:,1),A(:,2));

triplot(tri,A(:,1),A(:,2))

% Filtro de triangulos 


%% Calculo de descriptores
[r1,r2,r3,R,C,F,tr,tc,M,orientacion] = grothDescriptors(A,tri,epsilon);

%% Matching



