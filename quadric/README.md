Quadric package

This is a set of functions that work with quadric surfaces, rays, points, and their relations.

A quadric surface may defined in several different ways. In these routines, I adopt the forms:

Implicit:
       S(x,y,z) =  Ax^2 + Bx^2 + Cx^2 + 
                   2Dxy + 2Exz + 2Fyz +
                   2Gx + 2Hy + 2Iz + K == 0

 	Note that the order of the cross-terms is xy, xz, yz

Matrix:
       [A D E G;
        D B F H;
        E F C I;
        G H I K]

Vector:
       [A B C D E F G H I J K]


The directories are:
convert		- Convert quadrics between these forms, and from angles to rays  geodetics	- Points and paths on ellipsoidal surfacesplot		- Plot quadric surfacesprimitives	- Return primitive quadric surfaces (sphere, paraboloid, hyperboloid)properties	- Obtain the properties of a given quadric surfacerelations	- Spatial relationships between rays, surfaces, and pointstransform	- Affine transforms of surfaces, and rays.