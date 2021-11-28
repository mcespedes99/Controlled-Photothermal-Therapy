"""
Author: Mauricio Cespedes Tenorio
Date: Nov. 18th, 2021
Copyright: Laboratorio de Investigacion en Ingenieria Biomedica, UCR. 2021
Description: This code is used for the creation of the 3D geometry required for 
Montecarlo simulation. The geometry consists on a brick composed by three layers
of tissue with a spherical tumor inside of it. This code was based on the example code
"netgen_sphere_in_box.py" from ValoMC repository.
"""

#The Netgen/NGSolve package is used, which can be downloaded from: https://ngsolve.org/downloads
from netgen.csg import *
from ngsolve import Mesh
from ngsolve.solve import Draw

geo = CSGeometry()

# Create a box (layer) of muscle. This uses two opposite corners of the cube as reference
# with coordinates (x,y,z)
muscle = OrthoBrick(Pnt(4,-30,-30),Pnt(30,30,30))
# Create a box (layer) of adipose tissue
adipose = OrthoBrick(Pnt(2,-30,-30),Pnt(4,30,30))
# Create a box (layer) of skin
skin_o = OrthoBrick(Pnt(0,-30,-30),Pnt(2,30,30))
# Add a sphere into the middle (the maxh refers to the mesh max size)
sphere = Sphere(Pnt(0,0,0),5).maxh(2)

# Remesh a face of the cube so that it contains
# a cylindrical inclusion to simulate where the laser beam would go through
infcylinder = Cylinder(Pnt(-5, 0, 0), Pnt(0, 0, 0), 5); # infinite cylinder along x = -5 ... 0

# Create planes that cut the cylinder to make it finite.
# Two planes are required per layer to create intersections between 3D objects. 
plane1 = infcylinder*Plane (Pnt(0,0,0), Vec(-1,0,0) ); #One point of the plane and its normal vector are given as inputs
plane2 = infcylinder*Plane (Pnt(2,0,0), Vec(1,0,0) );
plane3 = infcylinder*Plane (Pnt(2,0,0), Vec(-1,0,0) );
plane4 = infcylinder*Plane (Pnt(4,0,0), Vec(1,0,0) );
plane5 = infcylinder*Plane (Pnt(4,0,0), Vec(-1,0,0) );
plane6 = infcylinder*Plane (Pnt(30,0,0), Vec(1,0,0) );
plane7 = infcylinder*Plane (Pnt(5,0,0), Vec(1,0,0) );

#Creation of the intersections between the infinite cylinder with each of the layers of tissue.
cylinder1 = (plane1*plane2*infcylinder);
cylinder2 = (plane3*plane4*infcylinder);
cylinder3 = (plane5*plane6*infcylinder);

# Creation of the official layers of tissue with all the intersections.
# 1. Skin: union of the cylinder1 with the skin brick - (minus) its intersection with the sphere.
sphere = (sphere*plane1*plane7);
skin = skin_o+cylinder1;
skin = skin-sphere;
# Command to give a nave to this part of the geometry (later used in MATLAB):
skin.mat("skin");
geo.Add(skin);

# 2. Subcutaneous fat: union of the cylinder2 with the adipose 
#brick - (minus) its intersection with the sphere.
adipose = adipose+cylinder2;
adipose = adipose-sphere;
adipose.mat("adipose");
geo.Add(adipose);

# 3. Muscle: union of the cylinder3 with the muscle brick - (minus) its 
# intersection with the sphere.
muscle = muscle+cylinder3;
muscle = muscle-sphere;
muscle.mat("muscle");
geo.Add(muscle);

# 4. Addition of the sphere to the geometry.
sphere.mat("sphere")
geo.Add(sphere);

# Creation of the mesh with a max size of 2 mm:
ngmesh = geo.GenerateMesh(maxh=2)

# Exporting the geometry to *.vol file.
ngmesh.Save("tejido_1mm_afuera.vol")