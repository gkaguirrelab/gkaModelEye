# gkaModelEye
These routines implement a ray-traced model eye in MATLAB routines. A primary application of the model is to describe the entrance pupil in the image plane for a rotated eye. The entrance pupil is described by the parameters of an ellipse fit to the pupil perimeter, and those parameters are given in "transparent" form.

The model is described in:

	GK Aguirre (2018) The Entrance Pupil of the Human Eye. bioRxiv.

These routines are support model-based eye tracking with transparentTrack: https://github.com/gkaguirrelab/transparentTrack

The function `pupilProjection_fwd` implements the forward model of this projection. Inputs to this routine are:
 * `eyePose` which is a vector that describes dynamic aspects of the eye, specifically rotation in degrees of azimuth, elevation, and torsion, and the radius of the pupil aperture in mm.
 * `sceneGeometry` which is a structure that describes static aspects of the scene, including the properties and position of a pinhole camera model and the eye. The sceneGeometry structure is generated by the function `createSceneGeometry`.

The function `pupilProjection_inv` implements a search over eyePose parameters and executions of the forward model to find the eyePose values that best describe an observed entrance pupil ellipse.

The forward model of the appearance of the pupil and iris accounts for the refractive properties of the cornea (and any artificial lenses between the eye and the camera). The routine `virtualImageFunc.m` calculates the effect of refraction, making use of calls to `rayTraceEllipsoids.m`. An improvement in the execution time of the forward model can be achieved by compiling the ray tracing routines. To do so, issue the command `compileVirtualImageFunc` at the MATLAB console. A compiled MEX file version of `virtualImageFunc` will be placed in the bin directory of this repository if it is not already present.

To install and configure gkaModelEye, first install toolboxToolbox (tBtB), which provides for declarative dependency management for Matlab (https://github.com/ToolboxHub/ToolboxToolbox). Once tBtB is installed, transparentTrack (and all its dependencies) will be installed and readied for use with the command `tbUse('gkaModelEye')`. If you do not wish to use tBtB, add the quadfit toolbox to your path (https://www.mathworks.com/matlabcentral/fileexchange/45356-fitting-quadratic-curves-and-surfaces). Additionally, to automatically run all examples, the ExampleTest toolbox is needed (https://github.com/isetbio/ExampleTestToolbox.git)

A hierarchy of the functions is as follows:
```
    pupilProjection_inv
            |
            V
    pupilProjection_fwd  <--  createSceneGeometry
            |                   |-- modelEyeParameters
            V                   |     '--returnRefractiveIndex
    virtualImageFunc            |
            |                   |-- addContactLens
            V                   '-- addSpectacleLens
    rayTraceEllipsoids
```

Most functions have associated examples in the header comments. This command issued in the MATLAB console will test all examples:
```
	[names,status] = RunExamples(fullfile(userpath(),'toolboxes','gkaModelEye'))
```

A good place to get started is to try rendering the model eye for a few views and examining the parameters of the pupil ellipse:
```
    sceneGeometry=createSceneGeometry();
    eyePose = [-30 -5 0 3];
    renderEyePose(eyePose, sceneGeometry);
    pupilEllipse = pupilProjection_fwd(eyePose,sceneGeometry);
```
