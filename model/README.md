## model

The code that implements the model eye. The directories are organized as follows:

- eye: The primary routines `modelEyeParameters` and `createSceneGeometry` are here, along with sub-directories that define the model properties for different "species" of eyes.
- project: Create a projection of the model eye to the image plane of a camera observing the eye
- bench: Routines that create, alter, validate, and ray-trace an optical system
- bin: Compiled versions of some of the ray-tracing routines
- quadric: functions for the manipulation of quadric surfaces, rays, and planes. Includes the function `rayTraceQuadrics`
- ellipse: functions for fitting ellipses to sets of points, and transforming between different forms of ellipse parameters. Includes the routine `eyePoseEllipseFit` which fits an ellipse to points on the perimeter of the entrance pupil in terms of the translation and rotation of a model eye.
- coord: transform between "eye" and "world" coordinate systems, and perform rotations of points and rays within these systems.
- calc: functions that derive optical properties from optical systems and eye models.
- plot: create depictions of eyes, optical systems, and scene elements