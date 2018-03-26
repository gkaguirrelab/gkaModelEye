function sceneGeometry = createSceneGeometry(varargin)
% Create and return a sceneGeometry structure
%
% Syntax:
%  sceneGeometry = createSceneGeometry()
%
% Description:
%   Using default values and passed key/value pairs, this routine creates a
%   sceneGeometry structure, with fields the describe a camera, an eye,
%   corrective lenses, and the geometric relationship between them. The
%   fields are:
%
%   intrinsicCameraMatrix - This matrix has the form:
%
%       [fx  s x0]
%       [0  fy y0]
%       [0   0  1]
%
%   where fx and fy are the focal lengths of the camera in the x and y
%   image dimensions, s is the axis skew, and x0, y0 define the principle
%   offset point. For a camera sensor with square pixels, fx = fy. Ideally,
%   skew should be zero. The principle offset point should be in the center
%   of the sensor. Units are traditionally in pixels. These values can be
%   empirically measured for a camera using a calibration approach
%   (https://www.mathworks.com/help/vision/ref/cameramatrix.html). Note
%   that the values in the matrix returned by the Matlab camera calibration
%   routine must be re-arranged to correspond to the X Y image dimension
%   specifications we use here.
%
%   radialDistortionVector - A two element vector of the form:
%
%       [k1 k2]
%
%   that models the radial distortion introduced the lens. This is an
%   empirically measured property of the camera system.
%
%   extrinsicTranslationVector - a vector of the form:
%
%       [x]
%       [y]
%       [z]
%
%   with the values specifying the location (horizontal, vertical, and
%   depth, respectively) of the principle offset point of the camera in mm
%   relative to the scene coordinate system. We define the origin of the
%   scene coordinate system to be x=0, y=0 along the optical axis of the
%   eye, and z=0 to be the apex of the corneal surface.
%
%   extrinsicRotationMatrix - A 3x3 identity matrix:
%
%       [1  0  0]
%       [0  1  0]
%       [0  0  1]
%
%   The projection of pupil circles from scene to image is invariant to
%   rotations of the camera matrix, so these valyes should not require
%   adjustment.
%
%   primaryPosition - A 1x3 vector of:
%
%       [eyeAzimuth, eyeElevation, eyeTorsion]
%
%   that specifies the rotation angles (in head fixed axes) for which the
%   eye is in primary position (as defined by Listing's Law). These values
%   may also be used to define the position at which the subject is
%   fixating the origin point of a stimulus array.
%
%   constraintTolerance - A scalar. The inverse projection from ellipse on
%   the image plane to eye params (azimuth, elevation) imposes a constraint
%   on how well the solution must match the shape of the ellipse (defined
%   by ellipse eccentricity and theta) and the area of the ellipse. This
%   constraint is expressed as a proportion of error, relative to either
%   the largest possible area in ellipse shape or an error in area equal to
%   the area of the ellipse itself (i.e., unity). If the constraint is made
%   too stringent, then in an effort to perfectly match the shape of the
%   ellipse, error will increase in matching the position of the ellipse
%   center. It should be noted that even in a noise-free simulation, it is
%   not possible to perfectly match ellipse shape and area while matching
%   ellipse center position, as the shape of the projection of the pupil
%   upon the image plane deviates from perfectly elliptical due to
%   perspective effects. We find that a value in the range 0.01 - 0.03
%   provides an acceptable compromise in empirical data.
%
%   eye -  A sub-structure returned by the function modelEyeParameters. The
%   parameters define the anatomical properties of the eye, including the
%   size and shape of the anterior and posterior chamber. These parameters
%   are adjusted for the measured spherical refractive error of the subject
%   and (optionally) measured axial length. Unmatched key-value pairs
%   passed to createSceneGeometry are passed to modelEyeParameters.
%
%   spectacleLens -  [optional] A set of fields that define a "single
%   vision" corrective lens. Only modeling of spherical correction is
%   supported.
%
%   contactLens -  [optional] A set of fields that define a corrective
%   contact lens. Only modeling of spherical correction is supported.
%
%   opticalSystem - An mx3 matrix that assembles information regarding the
%   m surfaces of the cornea and any corrective lenses into a format needed
%   for ray tracing.
%
% Inputs:
%   none
%
% Optional key/value pairs
%  'sceneGeometryFileName' - Full path to file 
%  'intrinsicCameraMatrix' - 3x3 matrix
%  'radialDistortionVector' - 1x2 vector of radial distortion parameters
%  'extrinsicRotationMatrix' - 3x3 matrix
%  'extrinsicTranslationVector' - 3x1 vector
%  'primaryPosition'      - 1x3 vector
%  'constraintTolerance'  - Scalar. Range 0-1. Typical value 0.01 - 0.03
%  'contactLens'          - Scalar or 1x2 vector, with values for the lens
%                           refraction in diopters, and (optionally) the
%                           index of refraction of the lens material. If
%                           left empty, no contact lens is added to the
%                           model.
%  'spectacleLens'        - Scalar, 1x2, or 1x3 vector, with values for the
%                           lens refraction in diopters, (optionally) the
%                           index of refraction of the lens material, and
%                           (optinally) the vertex distance in mm. If left
%                           empty, no contact lens is added to the model.
%  'medium'               - String, options include:
%                           {'air','water','vacuum'}. This sets the index
%                           of refraction of the medium between the eye an
%                           the camera.
%  'spectralDomain'       - String, options include {'vis','nir'}.
%                           This is the light domain within which imaging
%                           is being performed. The refractive indices vary
%                           based upon this choice.
%                           
%
% Outputs
%	sceneGeometry         - A structure that contains the components of the
%                           projection model.
%
% Examples:
%{
    % Create a scene geometry file for a myopic eye wearing a contact lens.
    % The key-value sphericalAmetropia is passed to modelEyeParameters
    % within the routine
    sceneGeometry = createSceneGeometry('sphericalAmetropia',-2,'contactLens',-2);
%}
%{
    % Create a scene geometry file for a hyperopic eye wearing spectacles
    % that provide appropriate correction when underwater.
    sceneGeometry = createSceneGeometry('sphericalAmetropia',-2,'spectacleLens',2,'medium','water');
    % Plot a figure that traces a ray arising from the optical axis at the
    % pupil plane, departing at 15 degrees.    
    figureFlag.zLim = [-10 20]; figureFlag.hLim = [-10 10];
    rayTraceCenteredSphericalSurfaces([sceneGeometry.eye.pupilCenter(1) 0], deg2rad(15), sceneGeometry.opticalSystem,figureFlag);
%}


%% input parser
p = inputParser; p.KeepUnmatched = true;

% Optional analysis params
p.addParameter('sceneGeometryFileName','', @(x)(isempty(x) | ischar(x)));
p.addParameter('intrinsicCameraMatrix',[2600 0 320; 0 2600 240; 0 0 1],@isnumeric);
p.addParameter('radialDistortionVector',[0 0],@isnumeric);
p.addParameter('extrinsicTranslationVector',[0; 0; 120],@isnumeric);
p.addParameter('extrinsicRotationMatrix',[1 0 0; 0 1 0; 0 0 1],@isnumeric);
p.addParameter('primaryPosition',[0 0 0],@isnumeric);
p.addParameter('constraintTolerance',0.02,@isscalar);
p.addParameter('contactLens',[], @(x)(isempty(x) | isnumeric(x)));
p.addParameter('spectacleLens',[], @(x)(isempty(x) | isnumeric(x)));
p.addParameter('medium','air',@ischar);
p.addParameter('spectralDomain','nir',@ischar);

% parse
p.parse(varargin{:})


%% assemble the sceneGeometry structure
% Values defined locally
sceneGeometry.intrinsicCameraMatrix = p.Results.intrinsicCameraMatrix;
sceneGeometry.radialDistortionVector = p.Results.radialDistortionVector;
sceneGeometry.extrinsicTranslationVector = p.Results.extrinsicTranslationVector;
sceneGeometry.extrinsicRotationMatrix = p.Results.extrinsicRotationMatrix;
sceneGeometry.primaryPosition = p.Results.primaryPosition;
sceneGeometry.constraintTolerance = p.Results.constraintTolerance;

% Values returned by the modelEyeParameters() routine
sceneGeometry.eye = modelEyeParameters('spectralDomain',p.Results.spectralDomain,varargin{:});

% Generate the opticalSystem matrix through the cornea. Note that we model
% the corneal surfaces as spheres, with a radius equal to the radius of
% curvature value R. To create a more accurate model, we would need to
% update rayTraceCenteredSphericalSurfaces() to model aspherical surfaces.
mediumRefractiveIndex = returnRefractiveIndex( p.Results.medium, p.Results.spectralDomain );

% The center of the cornea front surface is at a position equal to its
% radius of curvature, thus placing the apex of the front corneal surface
% at a z position of zero. The back surface is shifted back to produce
% the appropriate corneal thickness.
cornealThickness = -sceneGeometry.eye.corneaBackSurfaceCenter(1)-sceneGeometry.eye.corneaBackSurfaceRadii(1);
opticalSystem = [nan, nan, sceneGeometry.eye.aqueousRefractiveIndex; ...
    -sceneGeometry.eye.corneaBackSurfaceR-cornealThickness, -sceneGeometry.eye.corneaBackSurfaceR, sceneGeometry.eye.corneaRefractiveIndex; ...
    -sceneGeometry.eye.corneaFrontSurfaceR, -sceneGeometry.eye.corneaFrontSurfaceR, mediumRefractiveIndex];

% Add a contact lens if requested
if ~isempty(p.Results.contactLens)
    switch length(p.Results.contactLens)
        case 1
            lensRefractiveIndex=returnRefractiveIndex( 'hydrogel', p.Results.spectralDomain );
            [opticalSystem, pOutFun] = addContactLens(opticalSystem, p.Results.contactLens, 'lensRefractiveIndex', lensRefractiveIndex);
        case 2
            [opticalSystem, pOutFun] = addContactLens(opticalSystem, p.Results.contactLens(1), 'lensRefractiveIndex', p.Results.contactLens(2));
        otherwise
            error('The key-value pair contactLens is limited to two elements: [refractionDiopters, refractionIndex]');
    end
    sceneGeometry.contactLens = pOutFun.Results;
end

% Add a spectacle lens if requested
if ~isempty(p.Results.spectacleLens)
    switch length(p.Results.spectacleLens)
        case 1
            lensRefractiveIndex=returnRefractiveIndex( 'polycarbonate', p.Results.spectralDomain );
            [opticalSystem, pOutFun] = addSpectacleLens(opticalSystem, p.Results.spectacleLens, 'lensRefractiveIndex', lensRefractiveIndex);
        case 2
            [opticalSystem, pOutFun] = addSpectacleLens(opticalSystem, p.Results.spectacleLens(1), 'lensRefractiveIndex', p.Results.spectacleLens(2));
        case 3
            [opticalSystem, pOutFun] = addSpectacleLens(opticalSystem, p.Results.spectacleLens(1), 'lensRefractiveIndex', p.Results.spectacleLens(2),'lensVertexDistance', p.Results.spectacleLens(3));
        otherwise
            error('The key-value pair spectacleLens is limited to three elements: [refractionDiopters, refractionIndex, vertexDistance]');
    end
    sceneGeometry.spectacleLens = pOutFun.Results;
end

% Store the optical system
sceneGeometry.opticalSystem = opticalSystem;

% Save the meta data
sceneGeometry.meta.createSceneGeometry = p.Results;


%% Save the sceneGeometry file
if ~isempty(p.Results.sceneGeometryFileName)
    save(p.Results.sceneGeometryFileName,'sceneGeometry');
end

end % createSceneGeometry

