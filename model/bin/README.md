## bin

This directory holds compiled MEX function versions of the ray tracing routines `findPupilRay`, `findGlintRay` and `findNodalRay`, and a compiled version of the Levente Hunyadi function `ellipsefit_robust`.

If you need to compile the MEX function for your system, issue the commands `compileInverseRayTrace` and `compileEllipseFit` in the MATLAB console.

The compiled function is placed by default in the directory that contains `compileInverseRayTrace`, which should be:
```
gkaModelEye/bin
```

Pre-compiled mex files are provided for Apple Silicon, and the ray trace function are available for Apple Intel. It's easy to compile your own versions for your computer.