# bin
This directory holds compiled MEX function versions of virtualImageFunc.

If you need to compile the mex function for your system, issue the command `compileVirtualImageFunc` in the MATLAB console.

The compiled function is placed by default in the directory:
```
fullfile(userpath(),'toolboxes','transparentTrack/code/bin')
```

If your have transparentTrack in some other location, you can specify the bin directory in the compile call by using:
```
compileVirtualImageFunc('functionDirPath','myDirPath')
``` 