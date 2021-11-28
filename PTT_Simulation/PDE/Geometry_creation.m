%% Creation of geometry for PDE simulation (heat distribution)
% Author: Mauricio Cespedes Tenorio
% Date: Nov. 27th, 2021
% Copyright: Laboratorio de Investigacion en Ingenieria Biomedica, UCR. 2021
% Description: This code is used to create the geometry used for the PDE simulation
% to emulate heat distribution in the tissue. This code was based on the example 
% 'Sphere in Cube' from the MATLAB PDE toolbox. For further reference, see
% (https://www.mathworks.com/help/pdf_doc/pde/pde.pdf). For this case, only a cube with a 
% mesh of 0.2 mm max in space was used. Therefore, no actual geometry was used. 
%
clear;clc;
% 1) Creation of the PDE's that we need to solve the problem:
%Main PDE
model = createpde(1);

%% 2) Creation of Geometry:
% 2.i)Definition of the mesh size that is going to be used
[xg, yg, zg] = meshgrid(0:0.15:3, -3:0.15:3, -3:0.15:3);
xg = xg(:);
yg = yg(:);
zg = zg(:);
K = convhull(xg,yg,zg);
nodes = [xg';yg';zg'];
elements = K';
% 2.vii) Creation of a geometry with a spherical hole in the middle
g = geometryFromMesh(model,nodes,elements);
model.Geometry = g;
figure('Position',[10,10,800,400]);
subplot(1,2,1)
pdegplot(model,'FaceAlpha',0.25,'CellLabel','on')
title('Geometry with Cell Labels')
subplot(1,2,2)
pdegplot(model,'FaceAlpha',0.25,'FaceLabel','on')
title('Geometry with Face Labels')
mesh=generateMesh(model,'Hmax',0.2,'Hmin', 0.05);
figure
pdemesh(mesh)