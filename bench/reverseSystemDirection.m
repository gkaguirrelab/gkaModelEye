function opticalSystemOut = reverseSystemDirection( opticalSystemIn )
% Reverse the valid direction of ray tracing for an optical system
%
% Syntax:
%  opticalSystemOut = reverseSystemDirection( opticalSystemIn )
%
% Description:
%   The implementation of optical systems and ray tracing in this code
%   results in only one direction of ray tracing being available for a
%   given opticalSystem variable. This routine reverses an optical system
%   from the "eyeToCamera" to the "cameraToEye" orientation, or
%   vice-a-vera.
%
% Inputs:
%   opticalSystemIn       - Struct or matrix. If struct, must have the
%                           fields {opticalSystem, surfaceLabels,
%                           surfaceColors}. The matrix form is mx19. See
%                           "assembleOpticalSystem.m" for details. 
%
% Optional key/values pairs:
%   none
%
% Outputs:
%   opticalSystemOut      - Struct or matrix. Will match the type of
%                           opticalSystemIn
%
% Examples:
%{
    eye = modelEyeParameters();
    opticalSystemIn=assembleOpticalSystem(eye,'surfaceSetName','retinaToCamera','opticalSystemNumRows',[]);
    opticalSystemOut = reverseSystemDirection( opticalSystemIn );
    calcSystemDirection(opticalSystemOut);
%}


%% Prepeare the optical system
% If we were supplied a struct, extract the system components
if isstruct(opticalSystemIn)
    opticalSystemMatrix = opticalSystemIn.opticalSystem;
    surfaceLabels = opticalSystemIn.surfaceLabels;
    surfaceColors = opticalSystemIn.surfaceColors;
else
    opticalSystemMatrix = opticalSystemIn;
end

% Strip the optical system of nan rows. We will add these back before
% returning
numRows = size(opticalSystemMatrix,1);
opticalSystemMatrix = opticalSystemMatrix(sum(isnan(opticalSystemMatrix),2)~=size(opticalSystemMatrix,2),:);

% Check that the optical system is valid
systemDirection = calcSystemDirection(opticalSystemMatrix);
if ~contains(systemDirection,{'eyeToCamera','cameraToEye'})
    errorString = ['Not a valid opticalSystem: ' systemDirection];
    error('reverseSystemDirection:invalidSystemMatrix',errorString);
end

% Copy the input to output
opticalSystemMatrixOut = opticalSystemMatrix;

% Add a nan row to the end of the matrix
opticalSystemMatrixOut = [opticalSystemMatrixOut; nan(1,19)];

% Shift the refractive index column down one position
opticalSystemMatrixOut(2:end,19) = opticalSystemMatrixOut(1:end-1,19);

% Remove the initial "nan" row
opticalSystemMatrixOut = opticalSystemMatrixOut(2:end,:);

% Flip the sign of the "side" column
opticalSystemMatrixOut(:,11) = opticalSystemMatrixOut(:,11) * (-1);

% Reverse the matrix along the column direction
opticalSystemMatrixOut = flipud(opticalSystemMatrixOut);


%% Pad the optical system matrix
% The number of rows in the optical system matrix is set to a fixed value
% so that the compiled ray-tracing routines can expect a constant size for
% the input variables. The nan rows are stripped out at the time of ray
% tracing.
if numRows > size(opticalSystemMatrixOut,1)
    opticalSystemMatrixOut = [opticalSystemMatrixOut; ...
        nan(numRows-size(opticalSystemMatrixOut,1),19)];
end


%% Handle the opticalSystem structure on output
if isstruct(opticalSystemIn)
    opticalSystemOut.opticalSystem = opticalSystemMatrixOut;
    opticalSystemOut.surfaceColors = flipud([surfaceColors(2:end); {[nan nan nan]}]);
    opticalSystemOut.surfaceLabels = flipud([surfaceLabels(2:end); {'start'}]);
else
    opticalSystemOut = opticalSystemMatrixOut;
end
        
        
end % assembleOpticalSystem
