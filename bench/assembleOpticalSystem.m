function [opticalSystem, surfaceLabels, surfaceColors, magnification] = assembleOpticalSystem( eye, varargin )
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
%   pupil to medium path are returned.
%
%   The optical system of the eye (and any contact lens) is separate from
%   the optical system of the medium to the camera. This is to allow
%   separate behavior for eye rotation of the two systems.
%
% Inputs:
%   eye                   - A model eye structure.
%
% Optional key/values pairs:
%  'surfaceSetName'       - A string or char vector that from the set:
%                           {'retinaToMedium','retinaToStop',
%                            'stopToMedium','mediumToStop','stopToRetina',
%                            'mediumToRetina','mediumToCamera', ...
%                            'cameraToMedium','retinaToCamera', ...
%                            'cameraToRetina','glint'}
%  'cameraMedium'         - String, options include:
%                           {'air','water','vacuum'}. This sets the index
%                           of refraction of the medium between the eye and
%                           the medium.
%  'contactLens'          - Scalar or 1x2 vector, with values for the lens
%                           refraction in diopters, and (optionally) the
%                           index of refraction of the lens material. If
%                           left empty, no contact lens is added to the
%                           model.
%  'spectacleLens'        - Scalar, 1x2, 1x3, 1x4 vector, with values for 
%                           the lens refraction in diopters, (optionally)
%                           the index of refraction of the lens material,
%                           (optinally) the vertex distance in mm, and
%                           (optionally) the base curve in doiopters.
%                           Typically, the spectacleLens is the only item
%                           within the mediumToCamera or cameratToMedium
%                           systems.
%                           empty, no spectacle is added to the model.
%  'opticalSystemNumRows' - Scalar. The optical system is defined with a
%                           fixed number of rows. This is done so that the
%                           compiled ray trace routines can expect input
%                           parameters of fixed size. If there are m
%                           surfaces defined in the surfaceSet, then the
%                           optical system will have 100-m rows of nans. If
%                           the key-value is set to empty, no nan rows will
%                           be added.
%  'skipMagCalc'          - During calculation of magnification, this
%                           routine is called recursively. This flag is set
%                           to prevent a loop.
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
%   magnification         - Scalar. The magnification of the visual world
%                           experienced by the eye as a consequence of
%                           artificial lenses. Only calculated for the
%                           'retinaToMedium' system set.
%

%% input parser
p = inputParser; p.KeepUnmatched = true;

p.addRequired('eye',@isstruct);

% Optional analysis params
p.addParameter('surfaceSetName','stopToMedium',@ischar);
p.addParameter('cameraMedium','air',@ischar);
p.addParameter('contactLens',[], @(x)(isempty(x) | isnumeric(x)));
p.addParameter('spectacleLens',[], @(x)(isempty(x) | isnumeric(x)));
p.addParameter('opticalSystemNumRows',100,@isnumeric);
p.addParameter('skipMagCalc',false, @islogical);

% parse
p.parse(eye, varargin{:})

% Get the refractive index of the medium in which the medium resides
mediumRefractiveIndex = returnRefractiveIndex( p.Results.cameraMedium, eye.meta.spectralDomain );

% This variable is set to empty by default
magnification = struct();

% The optical system is always assembled in the eyeToMedium direction, but
% the reversed version is returned if that is what was requested
switch p.Results.surfaceSetName
    
    case {'retinaToCamera','cameraToRetina'}
        
        % We start in the vitreous chamber. Assign this refractive index
        opticalSystem = initializeOpticalSystem(returnRefractiveIndex( 'vitreous', eye.meta.spectralDomain ));
        
        % Add the vitreous chamber surface. As this has the same refractive
        % index as the first line of the optical system, this surface does
        % not induce any refraction.
        opticalSystem = [opticalSystem; ...
            [eye.retina.S eye.retina.side eye.retina.boundingBox eye.retina.mustIntersect returnRefractiveIndex( 'vitreous', eye.meta.spectralDomain )]];
        
        % Add the lens
        opticalSystem = [opticalSystem; ...
            [eye.lens.S eye.lens.side eye.lens.boundingBox eye.lens.mustIntersect [eye.lens.index; returnRefractiveIndex( 'aqueous', eye.meta.spectralDomain )]]];
        
        % Add the cornea
        opticalSystem = [opticalSystem; ...
            [eye.cornea.S eye.cornea.side eye.cornea.boundingBox eye.cornea.mustIntersect [eye.cornea.index; mediumRefractiveIndex]]];
        
        % Assemble the labels
        surfaceLabels = [{'vitreous'}; eye.retina.label; eye.lens.label; eye.cornea.label];
        
        % Assemble the surface plot colors
        surfaceColors = [{[nan nan nan]}; eye.retina.plot.color; eye.lens.plot.color; eye.cornea.plot.color];
        
        % Add a contact lens if requested
        if ~isempty(p.Results.contactLens)
            switch length(p.Results.contactLens)
                case 1
                    lensRefractiveIndex=returnRefractiveIndex( 'hydrogel', eye.meta.spectralDomain );
                    opticalSystem = addContactLens(opticalSystem, p.Results.contactLens, 'lensRefractiveIndex', lensRefractiveIndex,'cornealRotation',eye.cornea.rotation);
                case 2
                    opticalSystem = addContactLens(opticalSystem, p.Results.contactLens(1), 'lensRefractiveIndex', p.Results.contactLens(2),'cornealRotation',eye.cornea.rotation);
                otherwise
                    error('The key-value pair contactLens is limited to two elements: [refractionDiopters, refractionIndex]');
            end
            surfaceLabels = [surfaceLabels; {'contactLens'}; {'tearfilm'}];
            surfaceColors = [surfaceColors; {[.5 .5 .5]}; {'blue'}];
            
            % Calculate the magnification produced by this lens
            if ~p.Results.skipMagCalc
                magnification.contact = calcAngularMagnification(eye,'contactLens',p.Results.contactLens);
            end
        end
        
        % Add a spectacle lens if requested
        if ~isempty(p.Results.spectacleLens)
            switch length(p.Results.spectacleLens)
                case 1
                    lensRefractiveIndex=returnRefractiveIndex( 'polycarbonate', eye.meta.spectralDomain );
                    opticalSystem = addSpectacleLens(opticalSystem, p.Results.spectacleLens, 'lensRefractiveIndex', lensRefractiveIndex);
                case 2
                    opticalSystem = addSpectacleLens(opticalSystem, p.Results.spectacleLens(1), 'lensRefractiveIndex', p.Results.spectacleLens(2));
                case 3
                    opticalSystem = addSpectacleLens(opticalSystem, p.Results.spectacleLens(1), 'lensRefractiveIndex', p.Results.spectacleLens(2),'lensVertexDistance', p.Results.spectacleLens(3));
                case 4
                    opticalSystem = addSpectacleLens(opticalSystem, p.Results.spectacleLens(1), 'lensRefractiveIndex', p.Results.spectacleLens(2),'lensVertexDistance', p.Results.spectacleLens(3), 'baseCurve', p.Results.spectacleLens(4));
                otherwise
                    error('The key-value pair spectacleLens is limited to four elements: [refractionDiopters, refractionIndex, vertexDistance, baseCurve]');
            end
            surfaceLabels = [surfaceLabels; {'spectacleLensBack'}; {'spectacleLensFront'}];
            surfaceColors = [surfaceColors; {[.5 .5 .5]}; {[.5 .5 .5]}];

            % Calculate the magnification produced by the this lens
            if ~p.Results.skipMagCalc
                magnification.spectacle = calcAngularMagnification(eye,'spectacleLens',p.Results.spectacleLens);
            end
        end
        
        % Reverse the system if needed
        if strcmp(p.Results.surfaceSetName,'cameraToRetina')
            opticalSystem = reverseSystemDirection(opticalSystem);
            surfaceColors = flipud([surfaceColors(2:end); {[nan nan nan]}]);
            surfaceLabels = flipud([surfaceLabels(2:end); {'medium'}]);
        end
        
    case {'retinaToMedium','mediumToRetina'}
        
        % We start in the vitreous chamber. Assign this refractive index
        opticalSystem = initializeOpticalSystem(returnRefractiveIndex( 'vitreous', eye.meta.spectralDomain ));
        
        % Add the vitreous chamber surface. As this has the same refractive
        % index as the first line of the optical system, this surface does
        % not induce any refraction.
        opticalSystem = [opticalSystem; ...
            [eye.retina.S eye.retina.side eye.retina.boundingBox eye.retina.mustIntersect returnRefractiveIndex( 'vitreous', eye.meta.spectralDomain )]];
        
        % Add the lens
        opticalSystem = [opticalSystem; ...
            [eye.lens.S eye.lens.side eye.lens.boundingBox eye.lens.mustIntersect [eye.lens.index; returnRefractiveIndex( 'aqueous', eye.meta.spectralDomain )]]];
        
        % Add the cornea
        opticalSystem = [opticalSystem; ...
            [eye.cornea.S eye.cornea.side eye.cornea.boundingBox eye.cornea.mustIntersect [eye.cornea.index; mediumRefractiveIndex]]];
        
        % Assemble the labels
        surfaceLabels = [{'vitreous'}; eye.retina.label; eye.lens.label; eye.cornea.label];
        
        % Assemble the surface plot colors
        surfaceColors = [{[nan nan nan]}; eye.retina.plot.color; eye.lens.plot.color; eye.cornea.plot.color];
        
        % Add a contact lens if requested
        if ~isempty(p.Results.contactLens)
            switch length(p.Results.contactLens)
                case 1
                    lensRefractiveIndex=returnRefractiveIndex( 'hydrogel', eye.meta.spectralDomain );
                    opticalSystem = addContactLens(opticalSystem, p.Results.contactLens, 'lensRefractiveIndex', lensRefractiveIndex,'cornealRotation',eye.cornea.rotation);
                case 2
                    opticalSystem = addContactLens(opticalSystem, p.Results.contactLens(1), 'lensRefractiveIndex', p.Results.contactLens(2),'cornealRotation',eye.cornea.rotation);
                otherwise
                    error('The key-value pair contactLens is limited to two elements: [refractionDiopters, refractionIndex]');
            end
            surfaceLabels = [surfaceLabels; {'contactLens'}; {'tearfilm'}];
            surfaceColors = [surfaceColors; {[.5 .5 .5]}; {'blue'}];
            
            % Calculate the magnification produced by this lens
            if ~p.Results.skipMagCalc
                magnification.contact = calcAngularMagnification(eye,'contactLens',p.Results.contactLens);
            end
        end
                
        % Reverse the system if needed
        if strcmp(p.Results.surfaceSetName,'mediumToRetina')
            opticalSystem = reverseSystemDirection(opticalSystem);
            surfaceColors = flipud([surfaceColors(2:end); {[nan nan nan]}]);
            surfaceLabels = flipud([surfaceLabels(2:end); {'medium'}]);
        end
        
        
    case {'retinaToStop','stopToRetina'}
        
        % We start in the vitreous chamber. Assign this refractive index
        opticalSystem = initializeOpticalSystem(returnRefractiveIndex( 'vitreous', eye.meta.spectralDomain ));
        
        % Add the vitreous chamber surface. As this has the same refractive
        % index as the first line of the optical system, this surface does
        % not induce any refraction.
        
        % Start in the retina
        opticalSystem = [opticalSystem; ...
            [eye.retina.S eye.retina.side eye.retina.boundingBox eye.retina.mustIntersect returnRefractiveIndex( 'vitreous', eye.meta.spectralDomain )]];
        
        % Add the lens, ending in the aqueous medium
        opticalSystem = [opticalSystem; ...
            [eye.lens.S eye.lens.side eye.lens.boundingBox eye.lens.mustIntersect [eye.lens.index; returnRefractiveIndex( 'aqueous', eye.meta.spectralDomain )]]];
        
        % Assemble the labels
        surfaceLabels = [{'vitreous'}; eye.retina.label; eye.lens.label];
        
        % Assemble the surface plot colors
        surfaceColors = [{[nan nan nan]}; eye.retina.plot.color; eye.lens.plot.color];
        
        % Reverse the system if needed
        if strcmp(p.Results.surfaceSetName,'stopToRetina')
            opticalSystem = reverseSystemDirection(opticalSystem);
            surfaceColors = flipud([surfaceColors(2:end); {[nan nan nan]}]);
            surfaceLabels = flipud([surfaceLabels(2:end); {'aqueous'}]);
        end
        
        
    case {'stopToMedium','mediumToStop'}
        
        % We start in the aqueous. Assign this refractive index
        opticalSystem = initializeOpticalSystem(returnRefractiveIndex( 'aqueous', eye.meta.spectralDomain ));
        
        % Add the cornea
        opticalSystem = [opticalSystem; ...
            [eye.cornea.S eye.cornea.side eye.cornea.boundingBox eye.cornea.mustIntersect [eye.cornea.index; mediumRefractiveIndex]]];
        
        % Assemble the labels
        surfaceLabels = [{'aqueous'}; eye.cornea.label];
        
        % Assemble the surface plot colors
        surfaceColors = [{[nan nan nan]}; eye.cornea.plot.color];
        
        % Add a contact lens if requested
        if ~isempty(p.Results.contactLens)
            switch length(p.Results.contactLens)
                case 1
                    lensRefractiveIndex=returnRefractiveIndex( 'hydrogel', eye.meta.spectralDomain );
                    opticalSystem = addContactLens(opticalSystem, p.Results.contactLens, 'lensRefractiveIndex', lensRefractiveIndex,'cornealRotation',eye.cornea.rotation);
                case 2
                    opticalSystem = addContactLens(opticalSystem, p.Results.contactLens(1), 'lensRefractiveIndex', p.Results.contactLens(2),'cornealRotation',eye.cornea.rotation);
                otherwise
                    error('The key-value pair contactLens is limited to two elements: [refractionDiopters, refractionIndex]');
            end
            surfaceLabels = [surfaceLabels; {'contactLens'}; {'tearfilm'}];
            surfaceColors = [surfaceColors; {[.5 .5 .5]}; {'blue'}];
        end
        
        % Reverse the system if needed
        if strcmp(p.Results.surfaceSetName,'mediumToStop')
            opticalSystem = reverseSystemDirection(opticalSystem);
            surfaceColors = flipud([surfaceColors(2:end); {[nan nan nan]}]);
            surfaceLabels = flipud([surfaceLabels(2:end); {'medium'}]);
        end
        

	case {'mediumToCamera','cameraToMedium'}
        
        % We start in the medium. Assign this refractive index
        opticalSystem = initializeOpticalSystem(returnRefractiveIndex( p.Results.cameraMedium, eye.meta.spectralDomain ));

        % Assemble the labels
        surfaceLabels = [{'cameraMedium'}];
        
        % Assemble the surface plot colors
        surfaceColors = [{[nan nan nan]}];

        % Add a spectacle lens if requested
        if ~isempty(p.Results.spectacleLens)
            switch length(p.Results.spectacleLens)
                case 1
                    lensRefractiveIndex=returnRefractiveIndex( 'polycarbonate', eye.meta.spectralDomain );
                    opticalSystem = addSpectacleLens(opticalSystem, p.Results.spectacleLens, 'lensRefractiveIndex', lensRefractiveIndex);
                case 2
                    opticalSystem = addSpectacleLens(opticalSystem, p.Results.spectacleLens(1), 'lensRefractiveIndex', p.Results.spectacleLens(2));
                case 3
                    opticalSystem = addSpectacleLens(opticalSystem, p.Results.spectacleLens(1), 'lensRefractiveIndex', p.Results.spectacleLens(2),'lensVertexDistance', p.Results.spectacleLens(3));
                case 4
                    opticalSystem = addSpectacleLens(opticalSystem, p.Results.spectacleLens(1), 'lensRefractiveIndex', p.Results.spectacleLens(2),'lensVertexDistance', p.Results.spectacleLens(3), 'baseCurve', p.Results.spectacleLens(4));
                otherwise
                    error('The key-value pair spectacleLens is limited to four elements: [refractionDiopters, refractionIndex, vertexDistance, baseCurve]');
            end
            surfaceLabels = [surfaceLabels; {'spectacleLensBack'}; {'spectacleLensFront'}];
            surfaceColors = [surfaceColors; {[.5 .5 .5]}; {[.5 .5 .5]}];

            % Calculate the magnification produced by the this lens
            if ~p.Results.skipMagCalc
                magnification.spectacle = calcAngularMagnification(eye,'spectacleLens',p.Results.spectacleLens);
            end
            
        end
        
        % Reverse the system if needed
        if strcmp(p.Results.surfaceSetName,'cameraToMedium')
            opticalSystem = reverseSystemDirection(opticalSystem);
            surfaceColors = flipud([surfaceColors(2:end); {[nan nan nan]}]);
            surfaceLabels = flipud([surfaceLabels(2:end); {'medium'}]);
        end
        
        
    case {'glint'}
        
        % First assemble the path tearfilm --> medium
        
        % We start in the tearfilm
        opticalSystem = initializeOpticalSystem(returnRefractiveIndex( 'tears', eye.meta.spectralDomain ));
        
        % Add the tear surface
        opticalSystem = [opticalSystem; ...
            eye.cornea.S(3,:) eye.cornea.side(3) eye.cornea.boundingBox(3,:) 1 mediumRefractiveIndex];
        
        % Assemble the labels
        surfaceLabels = [eye.cornea.label(3); eye.cornea.label(3)];
        
        % Assemble the surface plot colors
        surfaceColors = [eye.cornea.plot.color(3); eye.cornea.plot.color(3)];
        
        % Add a contact lens if requested. We are going to model the glint
        % as taking place at the corneal tear film, regardless of the
        % presence of the contact lens. This might be introducing some
        % error, and I'll probably want to return to test if having the
        % glint off the contact lens tear film is more accurate.
        if ~isempty(p.Results.contactLens)
            switch length(p.Results.contactLens)
                case 1
                    lensRefractiveIndex=returnRefractiveIndex( 'hydrogel', eye.meta.spectralDomain );
                    opticalSystem = addContactLens(opticalSystem, p.Results.contactLens, 'lensRefractiveIndex', lensRefractiveIndex,'cornealRotation',eye.cornea.rotation);
                case 2
                    opticalSystem = addContactLens(opticalSystem, p.Results.contactLens(1), 'lensRefractiveIndex', p.Results.contactLens(2),'cornealRotation',eye.cornea.rotation);
                otherwise
                    error('The key-value pair contactLens is limited to two elements: [refractionDiopters, refractionIndex]');
            end
            surfaceLabels = [surfaceLabels; {'contactLens'}; {'tearfilm'}];
            surfaceColors = [surfaceColors; {[.5 .5 .5]}; {'blue'}];
        end

        % Reverse the system to give us a path from the medium to the
        % corneal tear film
        opticalSystem = reverseSystemDirection(opticalSystem);
        surfaceColors = flipud([surfaceColors(2:end); {[nan nan nan]}]);
        surfaceLabels = flipud([surfaceLabels(2:end); {'cameraMedium'}]);
        
        % Set the refractive index of the corneal tear film to be that of
        % the prior surface, with a negative value. This results in a ray
        % that is reflected and does not undergo refraction
        opticalSystem(end,19) = -opticalSystem(end-1,19);
        
        % If the optical system has more than two rows, that means we have
        % surfaces in between the medium and the tear film. In this case,
        % we have to add the surfaces heading back out.
        if size(opticalSystem,1)>2
            opticalSystemR = reverseSystemDirection(opticalSystem);
            surfaceColorsR = flipud([surfaceColors(2:end); {[nan nan nan]}]);
            surfaceLabelsR = flipud([surfaceLabels(2:end); eye.cornea.label(3)]);            
            
            opticalSystem = [opticalSystem; opticalSystemR(3:end,:)];
            surfaceColors = [surfaceColors; surfaceColorsR(3:end)];
            surfaceLabels = [surfaceLabels; surfaceLabelsR(3:end)];
        end
        
    otherwise
        error('Unrecognized surfaceSetName');
        
end % switch statement


%% Pad the optical system matrix
% The number of rows in the optical system matrix is set to a fixed value
% so that the compiled ray-tracing routines can expect a constant size for
% the input variables. The nan rows are stripped out at the time of ray
% tracing.
if ~isempty(p.Results.opticalSystemNumRows)
    osRowLength = size(opticalSystem,2);
    opticalSystem = [opticalSystem; ...
        nan(p.Results.opticalSystemNumRows-size(opticalSystem,1),osRowLength)];
end

end % assembleOpticalSystem
