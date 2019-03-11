function [opticalSystem, surfaceLabels, surfaceColors] = assembleOpticalSystem( eye, varargin )
% Create a matrix to be used for ray tracing a set of optical surfaces
%
% Syntax:
%  [opticalSystem, surfaceLabels, surfaceColors] = assembleOpticalSystem( eye )
%
% Description:
%   This routine assembles a matrix that contains the parameters that
%   define a set of optical surfaces for ray tracing. In addition to an eye
%   structure, the routine requires specification of the particular set of
%   surfaces to be assembled. By default, the surfaces that define the
%   pupil to camera path are returned.
%
% Inputs:
%   eye                   - A model eye structure.
%
% Optional key/values pairs:
%  'surfaceSetName'       - A string or char vector that from the set
%                           {'retinaToStop','stopToCamera',
%                           'retinaToCamera'}
%  'cameraMedium'         - String, options include:
%                           {'air','water','vacuum'}. This sets the index
%                           of refraction of the medium between the eye and
%                           the camera.
%  'contactLens'          - Scalar or 1x2 vector, with values for the lens
%                           refraction in diopters, and (optionally) the
%                           index of refraction of the lens material. If
%                           left empty, no contact lens is added to the
%                           model.
%  'spectacleLens'        - Scalar, 1x2, or 1x3 vector, with values for the
%                           lens refraction in diopters, (optionally) the
%                           index of refraction of the lens material, and
%                           (optinally) the vertex distance in mm. If left
%                           empty, no spectacle is added to the model.
%  'opticalSystemNumRows' - Scalar. The optical system is defined with a
%                           fixed number of rows. This is done so that the
%                           compiled ray trace routines can expect input
%                           parameters of fixed size. If there are m
%                           surfaces defined in the surfaceSet, then the
%                           optical system will have 100-m rows of nans.
%
% Outputs:
%   opticalSystem         - An mx19 matrix, where m is set by the key value
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
%                           from the matrix and have no effect.
%	surfaceLabels         - A cell array of strings or character vectors
%                           that identify each of the optical surfaces
%   surfaceColors         - A cell array of 3x1 vectors that provide the
%                           color specification for plotting each surface
%                           of the optical system.
%

%% input parser
p = inputParser; p.KeepUnmatched = true;

p.addRequired('eye',@isstruct);

% Optional analysis params
p.addParameter('surfaceSetName','stopToCamera',@ischar);
p.addParameter('cameraMedium','air',@ischar);
p.addParameter('contactLens',[], @(x)(isempty(x) | isnumeric(x)));
p.addParameter('spectacleLens',[], @(x)(isempty(x) | isnumeric(x)));
p.addParameter('opticalSystemNumRows',100,@isnumeric);

% parse
p.parse(eye, varargin{:})

% Get the refractive index of the medium in which the camera resides
mediumRefractiveIndex = returnRefractiveIndex( p.Results.cameraMedium, eye.meta.spectralDomain );

% Initialize an empty optical system matrix
opticalSystem=[];

switch p.Results.surfaceSetName
    case 'retinaToCamera'
        
        % Start on the retina with vitreous refractive index
        opticalSystem = [opticalSystem; ...
            [eye.retina.S eye.retina.side eye.retina.boundingBox eye.retina.mustIntersect returnRefractiveIndex( 'vitreous', eye.meta.spectralDomain )]];
        
        % Add the lens
        opticalSystem = [opticalSystem; ...
            [eye.lens.S eye.lens.side eye.lens.boundingBox eye.lens.mustIntersect [eye.lens.index; returnRefractiveIndex( 'aqueous', eye.meta.spectralDomain )]]];
        
        % Add the cornea
        opticalSystem = [opticalSystem; ...
            [eye.cornea.S eye.cornea.side eye.cornea.boundingBox eye.cornea.mustIntersect [eye.cornea.index; mediumRefractiveIndex]]];
        
        % Assemble the labels
        surfaceLabels = [eye.retina.label; eye.lens.label; eye.cornea.label];
        
        % Assemble the surface plot colors
        surfaceColors = [eye.retina.plot.color; eye.lens.plot.color; eye.cornea.plot.color];
        
        % Add a contact lens if requested
        if ~isempty(p.Results.contactLens)
            switch length(p.Results.contactLens)
                case 1
                    lensRefractiveIndex=returnRefractiveIndex( 'hydrogel', eye.meta.spectralDomain );
                    opticalSystem = addContactLens(opticalSystem, p.Results.contactLens, 'lensRefractiveIndex', lensRefractiveIndex);
                case 2
                    opticalSystem = addContactLens(opticalSystem, p.Results.contactLens(1), 'lensRefractiveIndex', p.Results.contactLens(2));
                otherwise
                    error('The key-value pair contactLens is limited to two elements: [refractionDiopters, refractionIndex]');
            end
            surfaceLabels = [surfaceLabels; {'contactLens'}];
            surfaceColors = [surfaceColors; {[.5 .5 .5]}];
        end
        
        % Add a spectacle lens if requested
        if ~isempty(p.Results.spectacleLens)
            switch length(p.Results.spectacleLens)
                case 1
                    lensRefractiveIndex=returnRefractiveIndex( 'polycarbonate', eye.meta.spectralDomain );
                    opticalSystem = addSpectacleLens(opticalSystem, p.Results.spectacleLens, 'lensRefractiveIndex', lensRefractiveIndex);
                case 2
                    opticalSystem = addSpectacleLens(opticalSystem, p.Results.spectacleLens, 'lensRefractiveIndex', p.Results.spectacleLens(2));
                case 3
                    opticalSystem = addSpectacleLens(opticalSystem, p.Results.spectacleLens, 'lensRefractiveIndex', p.Results.spectacleLens(2),'lensVertexDistance', p.Results.spectacleLens(3));
                otherwise
                    error('The key-value pair spectacleLens is limited to three elements: [refractionDiopters, refractionIndex, vertexDistance]');
            end
            surfaceLabels = [surfaceLabels; {'spectacleLensBack'}; {'spectacleLensFront'}];
            surfaceColors = [surfaceColors; {[.5 .5 .5]}; {[.5 .5 .5]}];
        end
        
    case 'retinaToStop'
        
        % Start in the retina
        opticalSystem = [opticalSystem; ...
            [eye.retina.S eye.retina.side eye.retina.boundingBox eye.retina.mustIntersect returnRefractiveIndex( 'vitreous', eye.meta.spectralDomain )]];
        
        % Add the lens, ending in the aqueous medium
        opticalSystem = [opticalSystem; ...
            [eye.lens.S eye.lens.side eye.lens.boundingBox eye.lens.mustIntersect [eye.lens.index; returnRefractiveIndex( 'aqueous', eye.meta.spectralDomain )]]];
        
        % Assemble the labels
        surfaceLabels = [eye.retina.label; eye.lens.label];
        
        % Assemble the surface plot colors
        surfaceColors = [eye.retina.plot.color; eye.lens.plot.color];
        
    case 'stopToCamera'
        
        % Start in the aqueous
        opticalSystem = [opticalSystem; ...
            [nan(1,10) nan nan(1,6) nan returnRefractiveIndex( 'aqueous', eye.meta.spectralDomain )]];
        
        % Add the cornea
        opticalSystem = [opticalSystem; ...
            [eye.cornea.S eye.cornea.side eye.cornea.boundingBox eye.cornea.mustIntersect [eye.cornea.index; mediumRefractiveIndex]]];
        
        % Assemble the labels
        surfaceLabels = [{'anteriorChamber'}; eye.cornea.label];
        
        % Assemble the surface plot colors
        surfaceColors = [{[nan nan nan]}; eye.cornea.plot.color];
        
        % Add a contact lens if requested
        if ~isempty(p.Results.contactLens)
            switch length(p.Results.contactLens)
                case 1
                    lensRefractiveIndex=returnRefractiveIndex( 'hydrogel', eye.meta.spectralDomain );
                    opticalSystem = addContactLens(opticalSystem, p.Results.contactLens, 'lensRefractiveIndex', lensRefractiveIndex);
                case 2
                    opticalSystem = addContactLens(opticalSystem, p.Results.contactLens(1), 'lensRefractiveIndex', p.Results.contactLens(2));
                otherwise
                    error('The key-value pair contactLens is limited to two elements: [refractionDiopters, refractionIndex]');
            end
            surfaceLabels = [surfaceLabels; {'contactLens'}];
            surfaceColors = [surfaceColors; {[.5 .5 .5]}];
        end
        
        % Add a spectacle lens if requested
        if ~isempty(p.Results.spectacleLens)
            switch length(p.Results.spectacleLens)
                case 1
                    lensRefractiveIndex=returnRefractiveIndex( 'polycarbonate', eye.meta.spectralDomain );
                    opticalSystem = addSpectacleLens(opticalSystem, p.Results.spectacleLens, 'lensRefractiveIndex', lensRefractiveIndex);
                case 2
                    opticalSystem = addSpectacleLens(opticalSystem, p.Results.spectacleLens, 'lensRefractiveIndex', p.Results.spectacleLens(2));
                case 3
                    opticalSystem = addSpectacleLens(opticalSystem, p.Results.spectacleLens, 'lensRefractiveIndex', p.Results.spectacleLens(2),'lensVertexDistance', p.Results.spectacleLens(3));
                otherwise
                    error('The key-value pair spectacleLens is limited to three elements: [refractionDiopters, refractionIndex, vertexDistance]');
            end
            surfaceLabels = [surfaceLabels; {'spectacleLensBack'}; {'spectacleLensFront'}];
            surfaceColors = [surfaceColors; {[.5 .5 .5]}; {[.5 .5 .5]}];
        end
        
    case 'cameraToStop'
        
        % Start in the camera medium
        opticalSystem = [opticalSystem; ...
            [nan(1,10) nan nan(1,6) nan mediumRefractiveIndex]];
        
        % Add the cornea, ending in the aqueous refractive index
        opticalSystem = [opticalSystem; ...
            [flipud(eye.cornea.S) flipud(eye.cornea.side)*(-1) flipud(eye.cornea.boundingBox) flipud(eye.cornea.mustIntersect) [flipud(eye.cornea.index); returnRefractiveIndex( 'aqueous', eye.meta.spectralDomain )]]];
        
        % Assemble the labels
        surfaceLabels = [{'camera'}; flipud(eye.cornea.label)];
        
        % Assemble the surface plot colors
        surfaceColors = [{[nan nan nan]}; flipud(eye.cornea.plot.color)];
        
    case 'stopToRetina'
        
        % Start in the aqueous
        opticalSystem = [opticalSystem; ...
            [nan(1,10) nan nan(1,6) nan returnRefractiveIndex( 'aqueous', eye.meta.spectralDomain )]];

        % Add the lens, ending in the vitreous medium
       opticalSystem = [opticalSystem; ...
          [flipud(eye.lens.S) flipud(eye.lens.side)*(-1) flipud(eye.lens.boundingBox) flipud(eye.lens.mustIntersect) [flipud(eye.lens.index); returnRefractiveIndex( 'vitreous', eye.meta.spectralDomain )]]];

        % Add the retina. We assign the vitreous refractive index so that
        % the final unit vector direction of the ray is unchanged. The ray
        % stops here.
       opticalSystem = [opticalSystem; ...
          [flipud(eye.retina.S) flipud(eye.retina.side)*(-1) flipud(eye.retina.boundingBox) flipud(eye.retina.mustIntersect) returnRefractiveIndex( 'vitreous', eye.meta.spectralDomain )]];

        % Assemble the labels
        surfaceLabels = [{'anteriorChamber'}; flipud(eye.lens.label); flipud(eye.retina.label)];
        
        % Assemble the surface plot colors
        surfaceColors = [{[nan nan nan]}; flipud(eye.lens.plot.color); flipud(eye.retina.plot.color)];

    case 'cameraToRetina'
        
        % Start in the camera medium
        opticalSystem = [opticalSystem; ...
            [nan(1,10) nan nan(1,6) nan mediumRefractiveIndex]];
        
        % Add the cornea, ending in the aqueous refractive index
        opticalSystem = [opticalSystem; ...
            [flipud(eye.cornea.S) flipud(eye.cornea.side)*(-1) flipud(eye.cornea.boundingBox) flipud(eye.cornea.mustIntersect) [flipud(eye.cornea.index); returnRefractiveIndex( 'aqueous', eye.meta.spectralDomain )]]];

        % Add the lens, ending in the vitreous medium
       opticalSystem = [opticalSystem; ...
          [flipud(eye.lens.S) flipud(eye.lens.side)*(-1) flipud(eye.lens.boundingBox) flipud(eye.lens.mustIntersect) [flipud(eye.lens.index); returnRefractiveIndex( 'vitreous', eye.meta.spectralDomain )]]];

        % Add the retina. We assign the vitreous refractive index so that
        % the final unit vector direction of the ray is unchanged. The ray
        % stops here.
       opticalSystem = [opticalSystem; ...
          [flipud(eye.retina.S) flipud(eye.retina.side)*(-1) flipud(eye.retina.boundingBox) flipud(eye.retina.mustIntersect) returnRefractiveIndex( 'vitreous', eye.meta.spectralDomain )]];

        % Assemble the labels
        surfaceLabels = [{'camera'}; flipud(eye.cornea.label); flipud(eye.lens.label); flipud(eye.retina.label)];
        
        % Assemble the surface plot colors
        surfaceColors = [{[nan nan nan]}; flipud(eye.cornea.plot.color); flipud(eye.lens.plot.color); flipud(eye.retina.plot.color)];        
        
    otherwise
        error('Unrecognized surfaceSetName');
        
end % switch statement


%% Pad the optical system matrix
% The number of rows in the optical system matrix is set to a fixed value
% so that the compiled ray-tracing routines can expect a constant size for
% the input variables. The nan rows are stripped out at the time of ray
% tracing.
osRowLength = size(opticalSystem,2);
opticalSystem = [opticalSystem; ...
    nan(p.Results.opticalSystemNumRows-size(opticalSystem,1),osRowLength)];

end % assembleOpticalSystem
