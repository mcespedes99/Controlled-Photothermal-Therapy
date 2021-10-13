%% Netgen 3D example
% This example demonstrates how to setup a 3D simulation on a geometry
% built with Netgen and run it. The geometry is a sphere within a
% cube. The sphere has a different refractive index than the cube. The
% mesh is created so that there is a circular domain within one
% surface of the cube.  You can view the Netgen Python source code
% <netgen_sphere_in_box.py here>
%
% Note: This is example requires a powerful computer to run.

%% Import the NetGen file


clear all;

if(exist('tejido_piel.vol', 'file') ~= 2)
    error('Could not find the mesh data file. Please run netgen netgen_sphere_in_box.py');
end

[vmcmesh regions region_names boundaries boundary_names] = importNetGenMesh('tejido_piel.vol', false);

%% Find indices from the file
muscle = cell2mat(regions(find(strcmp(region_names,'muscle'))));
adipose = cell2mat(regions(find(strcmp(region_names,'adipose'))));
skin = cell2mat(regions(find(strcmp(region_names,'skin'))));
skin_2 = regions(find(strcmp(region_names,'skin')));
sphere = cell2mat(regions(find(strcmp(region_names,'sphere'))));


%% Set optical coefficients
vmcmedium.absorption_coefficient(:) = 0.0214336;   
vmcmedium.scattering_coefficient(:) = 8.2729053;  
vmcmedium.scattering_anisotropy(:) = 0.933957;                            
vmcmedium.refractive_index(:) = 1.37;           

% Use the indices
vmcmedium.absorption_coefficient(muscle) = 0.0214336;   
vmcmedium.scattering_coefficient(muscle) = 8.2729053;  
vmcmedium.scattering_anisotropy(muscle) = 0.933957;                            
vmcmedium.refractive_index(muscle) = 1.37; 

vmcmedium.absorption_coefficient(skin) = 0.0207638;   
vmcmedium.scattering_coefficient(skin) = 5.1004766;  
vmcmedium.scattering_anisotropy(skin) = 0.715;                            
vmcmedium.refractive_index(skin) = 1.3773113; 

vmcmedium.absorption_coefficient(adipose) = 0.0083523;   
vmcmedium.scattering_coefficient(adipose) = 3.7088950;  
vmcmedium.scattering_anisotropy(adipose) = 0.715;                            
vmcmedium.refractive_index(adipose) = 1.44; 

vmcmedium.absorption_coefficient(sphere) = 0.0895438;   
vmcmedium.scattering_coefficient(sphere) = 16.7833342;  
vmcmedium.scattering_anisotropy(sphere) = 0.933957;                            
vmcmedium.refractive_index(sphere) = 1.37; 


%% Find boundary elements
% A circular domain for the light source was meshed (circle r = 1.0 at the
% face of the cube whose normal points to [-1 0 0]) but it is not contained
% as a separate boundary condition. We can use findBoundaries to find it
% manually.

lightsource1 = findBoundaries(vmcmesh, 'direction', [0 0 0 ], [-5 0 0], 5.2);
vmcboundary.lightsource(lightsource1) = {'direct'};
%lightsource2 = findBoundaries(vmcmesh, 'direction', [0 0 0 ], [-2.5 0 0], 0.4);
%vmcboundary.lightsource_position(lightsource1)=lightsource1;

%% Ploteo esfera + láser
figure
hold on

trimesh(vmcmesh.BH,vmcmesh.r(:,1),vmcmesh.r(:,2),vmcmesh.r(:,3),'facecolor', 'r','FaceAlpha',0.2);

% Highlight the location for the lightsource for the plot
trimesh(vmcmesh.BH(lightsource1,:),vmcmesh.r(:,1), vmcmesh.r(:,2),vmcmesh.r(:,3),'facecolor', 'b');
% Show the sphere
tetramesh(vmcmesh.H(sphere,:), vmcmesh.r);

title('Laser and tumor location');
xlabel('x [cm]');
ylabel('y [cm]');
zlabel('z [cm]');
view(-36,16);

hold off

%Plot adipose layer
figure
hold on

trimesh(vmcmesh.BH,vmcmesh.r(:,1),vmcmesh.r(:,2),vmcmesh.r(:,3),'facecolor', 'r','FaceAlpha',0.2);

%Show the adipose tissue
tetramesh(vmcmesh.H(adipose,:), vmcmesh.r);

title('Adipose tissue layer');
xlabel('x [cm]');
ylabel('y [cm]');
zlabel('z [cm]');
view(-36,16);

hold off

%Plot skin layer
figure
hold on

trimesh(vmcmesh.BH,vmcmesh.r(:,1),vmcmesh.r(:,2),vmcmesh.r(:,3),'facecolor', 'r','FaceAlpha',0.2);

%Show the adipose tissue
tetramesh(vmcmesh.H(skin,:), vmcmesh.r);

title('Skin layer');
xlabel('x [cm]');
ylabel('y [cm]');
zlabel('z [cm]');
view(-36,16);

hold off


%% Run the simulation
%options.photon_count=1e8;
%solution = ValoMC(vmcmesh, vmcmedium, vmcboundary,options);

%% Visualize the solution
% Visualizing large tetrahedral meshes is often cumbersome. Alternative,
% less power consuming option is to use exportX3D to export solution to X3D
% format and view the file using e.g. meshlab. See 'help exportX3D' for 
% more details.

figure
hold on
halfspace_elements = findElements(vmcmesh, 'halfspace', [0 0 0], [0 1 0]);
tetramesh(vmcmesh.H(halfspace_elements,:), vmcmesh.r, solution.element_fluence(halfspace_elements));
view(-10,10);
%
hold
xlabel('x [cm]');
ylabel('y [cm]');
zlabel('z [cm]');

cc = colorbar;                       
c.Label.String = 'Fluence [(W/cm^2)]';
