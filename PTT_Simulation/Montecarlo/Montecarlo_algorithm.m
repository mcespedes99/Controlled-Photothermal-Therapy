%% Montecarlo algorithm to simulate light propagation in tissue
% Author: Mauricio Cespedes Tenorio
% Date: Nov. 18th, 2021
% Copyright: Laboratorio de Investigacion en Ingenieria Biomedica, UCR. 2021
% Description: This code is an implementation of Montecarlo algorithm in MATLAB.
% To execute it, it is necessary to have the file 'tejido_piel.vol' within the same
% directory, which can be generated through the code 'Geometry_generation.py'.
% This code was based on the example code "netgentest3d.m" from ValoMC repository.
%
% Note: This is example requires a powerful computer to run.

%% Import the NetGen file


clear all;

% Import geometry from vol file.
if(exist('tejido_piel.vol', 'file') ~= 2)
    error('Could not find the mesh data file. Please run Geometry_generation.py');
end

[vmcmesh regions region_names boundaries boundary_names] = importNetGenMesh('tejido_piel.vol', false);

%% Find indices from the file. 
% This is how MATLAB knows the different figures or layers that are contained within the geometry.
muscle = cell2mat(regions(find(strcmp(region_names,'muscle'))));
adipose = cell2mat(regions(find(strcmp(region_names,'adipose'))));
skin = cell2mat(regions(find(strcmp(region_names,'skin'))));
skin_2 = regions(find(strcmp(region_names,'skin')));
sphere = cell2mat(regions(find(strcmp(region_names,'sphere'))));


%% Tissue optical coefficients:
% First, general absorption coeff are given to all the geometry. I chose muscle coeff since
% it is the type of tissue with the biggest volume within the geometry
vmcmedium.absorption_coefficient(:) = 0.0214336; %[1/mm]
vmcmedium.scattering_coefficient(:) = 8.2729053; %[1/mm]
vmcmedium.scattering_anisotropy(:) = 0.933957;                            
vmcmedium.refractive_index(:) = 1.37;           

% Use the indices to assign optical coefficients to each type of tissue.
vmcmedium.absorption_coefficient(muscle) = 0.0214336; %[1/mm]   
vmcmedium.scattering_coefficient(muscle) = 8.2729053; %[1/mm]
vmcmedium.scattering_anisotropy(muscle) = 0.933957;                            
vmcmedium.refractive_index(muscle) = 1.37; 

vmcmedium.absorption_coefficient(skin) = 0.0207638; %[1/mm]
vmcmedium.scattering_coefficient(skin) = 5.1004766; %[1/mm]
vmcmedium.scattering_anisotropy(skin) = 0.715;                            
vmcmedium.refractive_index(skin) = 1.3773113; 

vmcmedium.absorption_coefficient(adipose) = 0.0083523; %[1/mm]
vmcmedium.scattering_coefficient(adipose) = 3.7088950; %[1/mm]
vmcmedium.scattering_anisotropy(adipose) = 0.715;                            
vmcmedium.refractive_index(adipose) = 1.44; 

vmcmedium.absorption_coefficient(sphere) = 0.0895438;  %[1/mm]
vmcmedium.scattering_coefficient(sphere) = 16.7833342; %[1/mm]
vmcmedium.scattering_anisotropy(sphere) = 0.933957;                            
vmcmedium.refractive_index(sphere) = 1.37; 


%% Find boundary elements
% A circular domain for the light source was meshed (circle r = 5.0 at the
% face of the cube in the skin layer.
% Command to find a circle of r=5mm in the face of the cube.
% The actual radius used in the command was 5.2 to be able to find the correct part of the mesh
% containing the circle of r=5mm. 
lightsource1 = findBoundaries(vmcmesh, 'direction', [0 0 0 ], [-5 0 0], 5.2);
vmcboundary.lightsource(lightsource1) = {'direct'};

%% Plotting of the mesh containing the skin and the sphere:
figure
hold on

trimesh(vmcmesh.BH,vmcmesh.r(:,1),vmcmesh.r(:,2),vmcmesh.r(:,3),'facecolor', 'r','FaceAlpha',0.2);

% Highlight the location for the lightsource for the plot
trimesh(vmcmesh.BH(lightsource1,:),vmcmesh.r(:,1), vmcmesh.r(:,2),vmcmesh.r(:,3),'facecolor', 'b');
% Show the sphere
tetramesh(vmcmesh.H(sphere,:), vmcmesh.r);

title('Laser and tumor location');
xlabel('x [mm]');
ylabel('y [mm]');
zlabel('z [mm]');
view(-36,16);

hold off

%% Plot adipose layer
figure
hold on

trimesh(vmcmesh.BH,vmcmesh.r(:,1),vmcmesh.r(:,2),vmcmesh.r(:,3),'facecolor', 'r','FaceAlpha',0.2);

%Show the adipose tissue
tetramesh(vmcmesh.H(adipose,:), vmcmesh.r);

title('Adipose tissue layer');
xlabel('x [mm]');
ylabel('y [mm]');
zlabel('z [mm]');
view(-36,16);

hold off

%% Plot skin layer
figure
hold on

trimesh(vmcmesh.BH,vmcmesh.r(:,1),vmcmesh.r(:,2),vmcmesh.r(:,3),'facecolor', 'r','FaceAlpha',0.2);

%Show the adipose tissue
tetramesh(vmcmesh.H(skin,:), vmcmesh.r);

title('Skin layer');
xlabel('x [mm]');
ylabel('y [mm]');
zlabel('z [mm]');
view(-36,16);

hold off


%% Run the simulation. 
% Recommendation: first corroborate that the mesh and geometry are as desired since repeating
% the simulation is difficult due to its execution time.
options.photon_count=1e8;
solution = ValoMC(vmcmesh, vmcmedium, vmcboundary,options);

%% Visualize the solution
figure
hold on
halfspace_elements = findElements(vmcmesh, 'halfspace', [0 0 0], [0 1 0]);
tetramesh(vmcmesh.H(halfspace_elements,:), vmcmesh.r, solution.element_fluence(halfspace_elements));
view(-10,10);
%
hold
xlabel('x [mm]');
ylabel('y [mm]');
zlabel('z [mm]');

cc = colorbar;                       
c.Label.String = 'Fluence [(W/mm^2)]';
