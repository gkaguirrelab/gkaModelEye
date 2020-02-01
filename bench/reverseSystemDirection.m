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
%   opticalSystemIn       - An mx19 matrix, where m is set by the key value
%                           opticalSystemNumRows. Each row contains the 
%                           values:
%                               [S side bb must n]
%                           where:
%                               S     - 1x10 quadric surface vector
%                               side  - Scalar taking the value -1 or 1
%                                       that defines which of the two
%                                       points of intersection on the
%                                       quadric should be used as the
%                                       refractive surface.
%                               bb    - 1x6 vector defining the bounding
%                                       box within which the refractive
%                                       surface is present.
%                               must  - Scalar taking the value of 0 or 1,
%                                       where 1 indicates that the ray must
%                                       intersect the surface. If the ray
%                                       misses a required surface, the
%                                       routine exits with nans for the
%                                       outputRay.
%                               n     - Refractive index of the surface.
%                           The first row corresponds to the initial
%                           conditions of the ray. Thus, the refractive
%                           index value given in the first row specifies
%                           the index of the medium in which the ray
%                           arises. The other values for the first row are
%                           ignored. The matrix may have rows of all nans.
%                           These are used to define a fixed sized input
%                           variable for compiled code. They are removed
%                           from the matrix and have no effect.%
% Optional key/values pairs:
%   none
%
% Outputs:
%   opticalSystemOut        An mx19 matrix. This is the opticalSystemIn,
%                           now reversed in direction.
%
% Examples:
%{
    eye = modelEyeParameters();
    opticalSystemIn=assembleOpticalSystem(eye,'surfaceSetName','retinaToCamera','opticalSystemNumRows',[]);
    opticalSystemOut = reverseSystemDirection( opticalSystemIn );
    calcSystemDirection(opticalSystemOut)
%}

%% input parser
p = inputParser; p.KeepUnmatched = true;

% Required inputs
p.addRequired('opticalSystemIn',@isnumeric);

% parse
p.parse(opticalSystemIn)

% Strip the optical system of nan rows. We will add these back before
% returning
numRows = size(opticalSystemIn,1);
opticalSystemIn = opticalSystemIn(sum(isnan(opticalSystemIn),2)~=size(opticalSystemIn,2),:);

% Check that the optical system is valid
systemDirection = calcSystemDirection(opticalSystemIn);
if ~any(strcmp(systemDirection,{'eyeToCamera','cameraToEye'}))
    errorString = ['Not a valid opticalSystem: ' systemDirection];
    error('reverseSystemDirection:invalidSystemMatrix',errorString);
end

% Copy the input to output
opticalSystemOut = opticalSystemIn;

% Add a nan row to the end of the matrix
opticalSystemOut = [opticalSystemOut; nan(1,19)];

% Shift the refractive index column down one position
opticalSystemOut(2:end,19) = opticalSystemOut(1:end-1,19);

% Remove the initial "nan" row
opticalSystemOut = opticalSystemOut(2:end,:);

% Flip the sign of the "side" column
opticalSystemOut(:,11) = opticalSystemOut(:,11) * (-1);

% Reverse the matrix along the column direction
opticalSystemOut = flipud(opticalSystemOut);


%% Pad the optical system matrix
% The number of rows in the optical system matrix is set to a fixed value
% so that the compiled ray-tracing routines can expect a constant size for
% the input variables. The nan rows are stripped out at the time of ray
% tracing.
if numRows > size(opticalSystemOut,1)
    opticalSystemOut = [opticalSystemOut; ...
        nan(numRows-size(opticalSystemOut,1),19)];
end

end % assembleOpticalSystem
