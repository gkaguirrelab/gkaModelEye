function sceneGeometry = createSceneGeometry(varargin)
% Create and return a sceneGeometry structure
%
% Syntax:
%  sceneGeometry = createSceneGeometry()
%
% Description:
%   Using default values and passed key/value pairs, this routine creates a
%   sceneGeometry structure, with fields that describe a camera, an eye,
%   corrective lenses, and the geometric relationship between them. The
%   fields are:
%
%  'cameraIntrinsic' - A structure that defines the properties of a pinhole
%       camera. Sub-fields:
%
%      'matrix' - A 3x3 matrix of the form:
%
%               [fx  s x0; 0  fy y0; 0   0  1]
%
%           where fx and fy are the focal lengths of the camera in the x
%           and y image dimensions, s is the axis skew, and x0, y0 define
%           the principle offset point. For a camera sensor with square
%           pixels, fx = fy. Ideally, skew should be zero. The principle
%           offset point is usually in the center of the sensor. Units are
%           traditionally in pixels. These values can be empirically
%           measured for a camera using a resectioning approach
%           (https://www.mathworks.com/help/vision/ref/cameramatrix.html).
%           Note that the values in the matrix returned by the MATLAB
%           camera resectioning routine must be re-arranged to correspond
%           to the X Y image dimension specifications used here.
%
%      'radialDistortionVector' - A 1x2 vector that models the radial 
%           distortion introduced by the lens. This is an empirically 
%           measured property of the camera system.
%
%      'sensorResolution' - A 1x2 vector that provides the dimension of the
%           camera image in pixels, along the X and Y dimensions,
%           respectively.
%
%  'cameraPosition' - A structure that defines the spatial position of the
%       camera w.r.t. to world coordinates. This is used to assemble the
%       camera extrinsic matrix. Sub-fields:
%
%      'translation' - A 3x1 vector of the form [horizontal; vertical;
%           depth] in units of mm. Specifies the location of the nodal
%           point of the camera relative to the world coordinate system. We
%           define the origin of the world coordinate system to be x=0, y=0
%           along the optical axis of the eye with zero rotation, and z=0
%           to be the apex of the corneal surface.
%
%      'torsion' - Scalar in units of degrees that specifies the torsional
%           rotation of the camera relative to the z-axis of the world
%           coordinate space. Because eye rotations are set have a value of
%           zero when the camera axis is aligned with the pupil axis of the
%           eye, the camera rotation around the X and Y axes of the
%           coordinate system are not used. Rotation about the Z axis
%           (torsion) is meaningful, as the model eye has different
%           movement properties in the azimuthal and elevational
%           directions, and has a non circular exit pupil.
%
%       'glintSourceRelative' - A 3xn vector of the form [horizontal;
%           vertical; depth] in units of mm, with n equal to the number of
%           light sources. Specifies the relative location of an active
%           light source of a camera relative to the translation camera
%           position. This is the source of light for the modeled glint.
%
%  'screenPosition' - A structure that defines the spatial position of a
%           screen that the eye is fixating upon. Sub-fields:
%
%      'distance' - A scalar that specifies the distance (in mm) of the 
%           center of the screen relative to the corneal apex when the
%           optical axis of the eye is aligned with the center of the
%           screen.
%
%      'dimensions' - 1x2 vector in units of mm that provides the width and
%           height of the screen
%
%      'resolutions' - 1x2 vector in units of pixels for the width and
%           height.
%
%      'poseRegParams' - Structure. Contains fields that may be used to
%           calculate gaze position on a screen (in degrees of visual
%           angle) based upon the pose (rotation) of the eye. The sub-
%           field "meta" contains the formulat to be used, which is of the
%           form:
%
%               f = R * [azi; ele] + t
%
%      'vecRegParams' - Structure. Contains fields that may be used to
%           calculate gaze position on a screen (in degrees of visual
%           angle) based upon the vector between the center of the pupil
%           and a glint. The sub- field "meta" contains the formula to be
%           used, which is of the form:
%
%               f = s * R * [Xp - Xg; -1*(Yp-Yg)] + t
%
%           where Xp,Yp and Xg,Yg are the X and Y positions of the pupil
%           center and the glint center. The sign of the Y dimension is
%           reversed due to the sign convention of positive vertical
%           rotation of the eye producing a negative change in Y position.
%
%
%  'eye' - A structure that is returned by the function modelEyeParameters.
%       The parameters define the anatomical properties of the eye. These
%       parameters are adjusted for the measured spherical refractive error
%       of the subject or measured axial length. Unmatched key-value pairs
%       passed to createSceneGeometry are passed to modelEyeParameters.
%
%  'refraction' - A structure with sub-fields for ray-tracing through sets
%       of optical surfaces. Standard fields are 'retinaToStop',
%       'stopToMedium', and 'retinaToMedium'. Each subset has the fields:
%
%         'opticalSystem' - An mx19 matrix, where m is the number of
%                           surfaces in the model, including the initial
%                           state of the ray. The matrix may have rows of
%                           all nans. These are used to define a fixed
%                           sized input variable for compiled code. They
%                           are removed from the matrix and have no effect.
%                           Further details in assembleOpticalSystem.
%         'surfaceLabels' - A cell array of strings or character vectors
%                           that identify each of the optical surfaces
%         'surfaceColors' - A cell array of 3x1 vectors that provide the
%                           color specification for plotting each surface
%                           of the optical system.
%         'magnification' - An optional field that is only populated if
%                           the sceneGeometry includes artificial lenses 
%                           (e.g., spectacles or contacts). This field
%                           provides the angular magnification of the world
%                           as experienced by the eye through the lens.
%
%  'meta' - A structure that contains information regarding the creation
%       and modification of the sceneGeometry.
%
% Inputs:
%   none
%
% Optional key/value pairs
%  'sceneGeometryFileName' - Full path to file
%  'intrinsicCameraMatrix' - 3x3 matrix
%  'radialDistortionVector' - 1x2 vector of radial distortion parameters
%  'cameraTranslation'    - 3x1 vector
%  'cameraRotation'       - 1x3 vector
%  'constraintTolerance'  - Scalar. Range 0-1. Typical value 0.01 - 0.10
%  'contactLens'          - Scalar or 1x2 vector, with values for the lens
%                           refraction in diopters, and (optionally) the
%                           index of refraction of the lens material. If
%                           left empty, no contact lens is added to the
%                           model.
%  'spectacleLens'        - Scalar, 1x2, 1x3, or 1x4 vector, with values for the
%                           lens refraction in diopters, (optionally) the
%                           index of refraction of the lens material,
%                           (optinally) the vertex distance in mm, and
%                           (optionally) the base curvature. If left empty,
%                           no spectacle is added to the model.
%  'cameraMedium'         - String, options include:
%                           {'air','water','vacuum'}. This sets the index
%                           of refraction of the medium between the eye and
%                           the camera.
%  'spectralDomain'       - String or numerical scalar. This is the
%                           wavelength domain within which imaging is being
%                           performed. The refractive indices vary based
%                           upon this choice. Either a wavelength (in nm)
%                           may be provided, or one of the char vectors:
%                           {'vis','nir'}.
%
% Outputs
%	sceneGeometry         - A structure.
%
% Examples:
%{
    % Create a sceneGeometry file for a myopic eye wearing a contact lens.
    % The key-value sphericalAmetropia is passed to modelEyeParameters
    % within the routine
    sceneGeometry = createSceneGeometry('sphericalAmetropia',-2,'contactLens',-2);
%}



%% input parser
p = inputParser; p.KeepUnmatched = true; p.PartialMatching = false;

% Optional analysis params
p.addParameter('sceneGeometryFileName','', @(x)(isempty(x) | ischar(x)));
p.addParameter('intrinsicCameraMatrix',[2600 0 320; 0 2600 240; 0 0 1],@isnumeric);
p.addParameter('sensorResolution',[640, 480],@isnumeric);
p.addParameter('radialDistortionVector',[0, 0],@isnumeric);
p.addParameter('cameraTranslation',[0; 0; 120],@isnumeric);
p.addParameter('cameraGlintSourceRelative',[-14; 0; 0],@isnumeric);
p.addParameter('cameraTorsion',0,@isnumeric);
p.addParameter('screenDistance',1065,@isnumeric);
p.addParameter('screenDimensions',[697.347,392.257],@isnumeric);
p.addParameter('screenResolutions',[1920,1080],@isnumeric);
p.addParameter('poseRegParams',[],@isstruct);
p.addParameter('vecRegParams',[],@isstruct);
p.addParameter('eye',[],@isstruct);
p.addParameter('surfaceSetName',{'retinaToCamera','cameraToRetina','retinaToStop','stopToMedium','retinaToMedium','mediumToRetina','mediumToCamera','cameraToMedium','glint'},@iscell);
p.addParameter('contactLens',[], @(x)(isempty(x) | isnumeric(x)));
p.addParameter('spectacleLens',[], @(x)(isempty(x) | isnumeric(x)));
p.addParameter('cameraMedium','air',@ischar);
p.addParameter('spectralDomain','nir',@(x)(ischar(x) || isnumeric(x)));

% parse
p.parse(varargin{:})

% Obtain the unmatched varargin
vararginUnmatched = namedargs2cell(p.Unmatched);

%% cameraIntrinsic
sceneGeometry.cameraIntrinsic.matrix = p.Results.intrinsicCameraMatrix;
sceneGeometry.cameraIntrinsic.radialDistortion = p.Results.radialDistortionVector;
sceneGeometry.cameraIntrinsic.sensorResolution = p.Results.sensorResolution;


%% cameraPosition
sceneGeometry.cameraPosition.translation = p.Results.cameraTranslation;
sceneGeometry.cameraPosition.torsion = p.Results.cameraTorsion;
sceneGeometry.cameraPosition.glintSourceRelative = p.Results.cameraGlintSourceRelative;


%% screenPosition
sceneGeometry.screenPosition.screenDistance = p.Results.screenDistance;
sceneGeometry.screenPosition.dimensions = p.Results.screenDimensions;
sceneGeometry.screenPosition.resolutions = p.Results.screenResolutions;
if isempty(p.Results.poseRegParams)
    sceneGeometry.screenPosition.poseRegParams.R = [1 0; 0 1];
    sceneGeometry.screenPosition.poseRegParams.t = [0; 0];
    sceneGeometry.screenPosition.poseRegParams.s = 1;
    sceneGeometry.screenPosition.poseRegParams.M = [1 0 0; 0 1 0; 0 0 1];
    sceneGeometry.screenPosition.poseRegParams.theta = 0;
    sceneGeometry.screenPosition.poseRegParams.meta = 'f = R * [azi; ele] + t';
    
else
    sceneGeometry.screenPosition.poseRegParams = p.Results.poseRegParams;
end
sceneGeometry.screenPosition.vecRegParams = p.Results.vecRegParams;


%% eye
if isempty(p.Results.eye)
    sceneGeometry.eye = modelEyeParameters('spectralDomain',p.Results.spectralDomain,vararginUnmatched{:});
else
    sceneGeometry.eye = p.Results.eye;
end

%% refraction
for ii = 1:length(p.Results.surfaceSetName)
    [opticalSystem, surfaceLabels, surfaceColors, magnification] = ...
        assembleOpticalSystem( sceneGeometry.eye, ...
        'surfaceSetName', p.Results.surfaceSetName{ii}, ...
        'cameraMedium', p.Results.cameraMedium, ...
        'contactLens', p.Results.contactLens, ...
        'spectacleLens', p.Results.spectacleLens, ...
        vararginUnmatched{:});
    sceneGeometry.refraction.(p.Results.surfaceSetName{ii}).opticalSystem = opticalSystem;
    sceneGeometry.refraction.(p.Results.surfaceSetName{ii}).surfaceLabels = surfaceLabels;
    sceneGeometry.refraction.(p.Results.surfaceSetName{ii}).surfaceColors = surfaceColors;
    if ~isempty(magnification)
        sceneGeometry.refraction.(p.Results.surfaceSetName{ii}).magnification = magnification;
    end
end


%% meta
sceneGeometry.meta.createSceneGeometry = p.Results;
sceneGeometry.meta.createSceneGeometry.varargin = varargin;

%% Save the sceneGeometry file
if ~isempty(p.Results.sceneGeometryFileName)
    save(p.Results.sceneGeometryFileName,'sceneGeometry');
end

end % createSceneGeometry

