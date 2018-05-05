# geometry
A set of functions borrowed from other sources:
 * `createLine3d`, `intersectLineSphere` - by David Legland and the geom3d library (https://www.mathworks.com/matlabcentral/fileexchange/24484-geom3d). Small modifications to the functions to allow compilation using codegen.
 * `EllipsoidPlaneIntersection` - by Sebahattin Bektas (https://www.mathworks.com/matlabcentral/fileexchange/52958-intersection-ellipsoid-and-a-plane)
 * `AtoG` - by Hui Ma (https://www.mathworks.com/matlabcentral/fileexchange/32105-conversion-of-conics-parameters). Required by `EllipsoidPlaneIntersection`.
