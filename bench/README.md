## bench

These routines for create, alter, and ray-trace an optical system.

An optical system is a set of surfaces with which light rays interact. In this code, the `opticalSystem` variable is an m x 19 matrix, with each row describing one of the m surfaces. The 19 values in a row specify (in order):

- `S` - 1x10 vector that defines a quadric surface (e.g., ellipsoid, hyperboloid). More information on quadric surfaces can be found in the `quadric` directory.
- `side` - Scalar taking the value -1 or 1 that defines which of the two possible points of intersection of a ray with a quadric should be used as the refractive surface.
- `bb` - 1x6 vector defining the bounding box (in the eye coordinate system) within which the refractive surface is present.
- `must` - Scalar taking the value of 0 or 1, where 1 indicates that the ray must intersect the surface if the ray trace is to continue to the next surface.
- `n` - Scalar which provides the relative refractive index of the surface. A negative value for the refractive index indicates that the ray trace routine should reflect, instead of refract, the ray when it intersects the surface.

The first row of the opticalSystem corresponds to the initial conditions of the ray. Thus, the refractive index value given in the first row specifies the index of the medium in which the ray arises. The other values for the first row are typically left as nans and ignored.

Ray tracing is conventionally performed from "right-to-left". Because we are interested in describing the appearance of an eye as viewed from a camera, the optical system is constructed by default in the order retina --> camera. An opticalSystem variable is only valid for ray-tracing in one direction at a time. The function `reverseSystemDirection` may be used to switch a system from the retina --> camera direction to the camera --> retina direction, and vice-a-versa.

The opticalSystem matrix may be stored in an opticalSystem structure, which contains as well a set of labels for the surfaces, and information regarding how the surfaces should be plotted.

The routine `findPupilRay` identifies a ray that leaves a location on the border of the aperture stop of the iris in a rotated eye and arrives at the nodal point of a camera that is observing the eye. This routine is used to account for the refractive effects of the cornea and any artificial lenses to generate the entrance pupil.

The routine `findGlintRay` identifies a ray that leaves a light source adjacent to the camera, is reflected from the tear film in a rotated eye, and arrives at the nodal point of the camera. This routine is used to identify the location of any glints in an image of the eye.

The contents of the directory are:

- assembleOpticalSystem.m - Given an eye structure (SEE: `modelEyeParameters.m`) returns a opticalSystem. The key `surfaceSetName` specifies which opticalSystem to generate (e.g., stopToCamera)- addBiconvexLens.m, addContactLens.m, addSpectacleLens.m, addStopAfter.m - Adds these surfaces to an optical system- findGlintRay.m, findPupilRay.m - Implements "inverse" ray-tracing. The compiled, MEX versions of these routines are stored in the `bin` directory.- reverseSystemDirection.m - Swaps the valid direction for ray-tracing of an opticalSystem variable.