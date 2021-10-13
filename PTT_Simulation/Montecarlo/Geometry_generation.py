from netgen.csg import *
from ngsolve import Mesh
from ngsolve.solve import Draw

geo = CSGeometry()

# Create a box as muscle
muscle = OrthoBrick(Pnt(4,-30,-30),Pnt(30,30,30))
# Create a box as adipose tissue
adipose = OrthoBrick(Pnt(2,-30,-30),Pnt(4,30,30))
# Create a box as skin
skin_o = OrthoBrick(Pnt(0,-30,-30),Pnt(2,30,30))
# Add a sphere into the middle
sphere = Sphere(Pnt(0,0,0),5).maxh(2)

# Remesh a face of the cube so that it contains
# a cylindrical inclusion

infcylinder = Cylinder(Pnt(-5, 0, 0), Pnt(0, 0, 0), 5); # infinite cylinder along x = -5 ... 5

# Create two planes that cut the cylinder to make it finite
plane1 = infcylinder*Plane (Pnt(0,0,0), Vec(-1,0,0) );
plane2 = infcylinder*Plane (Pnt(2,0,0), Vec(1,0,0) );
plane3 = infcylinder*Plane (Pnt(2,0,0), Vec(-1,0,0) );
plane4 = infcylinder*Plane (Pnt(4,0,0), Vec(1,0,0) );
plane5 = infcylinder*Plane (Pnt(4,0,0), Vec(-1,0,0) );
plane6 = infcylinder*Plane (Pnt(30,0,0), Vec(1,0,0) );
plane7 = infcylinder*Plane (Pnt(5,0,0), Vec(1,0,0) );

cylinder1 = (plane1*plane2*infcylinder);
cylinder2 = (plane3*plane4*infcylinder);
cylinder3 = (plane5*plane6*infcylinder);

# Boundary conditions can be set using e.g.
# geo.Add(brick.bc("lightsource"));. See NetGen documentation for more details.
# However, it is currently difficult apply them to a domain of a face
# hence not done in this example
sphere = (sphere*plane1*plane7);
skin = skin_o+cylinder1;
skin = skin-sphere;
skin.mat("skin");

geo.Add(skin);

adipose = adipose+cylinder2;
adipose = adipose-sphere;
adipose.mat("adipose");
geo.Add(adipose);

muscle = muscle+cylinder3;
muscle = muscle-sphere;
muscle.mat("muscle");
geo.Add(muscle);

sphere.mat("sphere")
geo.Add(sphere);

ngmesh = geo.GenerateMesh(maxh=2)
ngmesh.Save("tejido_1mm_afuera.vol")
