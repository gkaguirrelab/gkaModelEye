project

These routines create a projection of the model eye in the image plane of a camera. The main routine (projectModelEye.m) is at the top level. The stages of the projection are inside the "stages" directory, and are (in sequence): 

- addPseudoTorsion.m- addStopPerimeter.m- addFullEyeModel.m- refractEyePoints.m- addGlint.m- applyEyeRotation.m- projectToImagePlane.m- applyRadialLensDistortion.m
- obtainImagePlaneEllipse.m- obtainGlintCoord.m