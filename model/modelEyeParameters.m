function eye = modelEyeParameters( varargin )
% Return the parameters of a model eye
%
% Syntax:
%  eye = modelEyeParameters()
%
% Description:
%   This routine returns the parameters of a model eye used in the
%   sceneGeometry routines.
%
%   The parameters returned by this routine correspond to the eyeWorld
%   coordinate space used in pupilProjection_fwd, which is relative to the
%   optical / pupillary axis, with the apex of the cornea set as zero in
%   depth. The space has the dimensions [depth, horizontal, vertical];
%   negative values of depth are towards the back of the eye. The model
%   assumes the optical and pupil axis of the eye are algined.
%
% Inputs:
%   none
%
% Optional key/value pairs:
%  'sphericalAmetropia'   - Scalar, in units of diopters. Several
%                           parameters of the model eye are adjusted by the
%                           spherical refractive error of the eye. A
%                           negative number is the correction that would be
%                           used for a myopic person.
%  'axialLength'          - Scalar. This is the axial length along the 
%                           optical axis. When set, this fixes the axial
%                           length of the eye to the passed value in
%                           millimeters. As the modeled anterior chamber
%                           depth is not variable, this change is enforced
%                           on the vitreous chamber. The remaining
%                           dimensions of the vitreous chamber are scaled
%                           to fit the proportions predicted by the
%                           Atchison model for the specified degree of
%                           ametropia.
%  'eyeLaterality'        - A text string that specifies which eye (left,
%                           right) to model. Allowed values (in any case)
%                           are {'left','right','L','R','OS','OD'}
%  'species'              - A text string that specifies the species to be
%                           modeled. Supported values (in any case) are
%                           {'human','dog'}
%  'spectralDomain'       - String, options include {'vis','nir'}.
%                           This is the wavelength domain within which
%                           imaging is being performed. The refractive
%                           indices vary based upon this choice.
%  'visualAxisDegRetina'  - 1x3 vector. This is the position of the fovea 
%                           w.r.t. to optical axis in degrees of retina.
%                           The values are [azimuth, elevation, torsion].
%                           Used in model development.
%  'opticDiscAxisDegRegina'  - 1x3 vector. This is the position of the  
%                           optic disc w.r.t. to optical axis in degrees of
%                           retina. The values are [azimuth, elevation,
%                           torsion]. Used in model development.
%
% Outputs:
%   eye                   - A structure with fields that contain the values
%                           for the model eye.
%
% Examples:
%{
    % Default parameters, corresponding to an emmetropic, right, human eye
    eye = modelEyeParameters();
%}
%{
    % Parameters for an myopic (-3), left, human eye
    eye = modelEyeParameters('sphericalAmetropia',-3,'eyeLaterality','left');
%}


%% input parser
p = inputParser; p.KeepUnmatched = true;

% Optional
p.addParameter('sphericalAmetropia',0,@isscalar);
p.addParameter('axialLength',[],@(x)(isempty(x) || isscalar(x)));
p.addParameter('cornealAxis',[],@(x)(isempty(x) || isnumeric(x)));
p.addParameter('eyeLaterality','Right',@ischar);
p.addParameter('species','Human',@ischar);
p.addParameter('spectralDomain','nir',@ischar);
p.addParameter('skipEyeAxes',false,@islogical);

% parse
p.parse(varargin{:})

% Interpret the passed laterality
switch p.Results.eyeLaterality
    case {'right','RIGHT','Right','R','r','od','OD'}
        eyeLaterality = 'Right';
    case {'left','LEFT','Left','L','l','os','OS'}
        eyeLaterality = 'Left';
    otherwise
        error('Please specify a valid eye laterality for the model eye');
end

% Create an empty eye struct
eye = struct();

% Meta data regarding the units of the model
eye.meta.p = p.Results;
eye.meta.units = 'mm';
eye.meta.coordinates = 'eyeWorld';
eye.meta.dimensions = {'depth (axial)' 'horizontal' 'vertical'};
eye.meta.eyeLaterality = eyeLaterality;
eye.meta.sphericalAmetropia = p.Results.sphericalAmetropia;
eye.meta.species = p.Results.species;
eye.meta.spectralDomain = p.Results.spectralDomain;

% Switch parameters at the top level by species
switch eye.meta.species

    %% Human eye
    case {'human','Human','HUMAN'}
                
        % Cornea
        eye.cornea = human.cornea(eye, p.Results.cornealAxis);

        % Iris
        eye.iris = human.iris(eye);

        % Pupil
        eye.pupil = human.pupil(eye);

        % retina
        eye.retina = human.retina(eye);

        % Lens
        eye.lens = human.lens(eye);

        % Rotation centers
        eye.rotationCenters = human.rotationCenters(eye);

        % Axes
        if ~p.Results.skipEyeAxes
           eye.axes = human.axes(eye);
        end
        
        %% Refractive indices
        % Obtain refractive index values for this spectral domain.
        eye.index.vitreous = returnRefractiveIndex( 'vitreous', p.Results.spectralDomain );
        eye.index.aqueous = returnRefractiveIndex( 'aqueous', p.Results.spectralDomain );

        
    %% Dog eye
    case {'dog','Dog','canine','Canine'}
        
        % Unless othewise stated, values taken from:
        %	Coile, D. C., and L. P. O'Keefe. "Schematic eyes for domestic
        %	animals." Ophthalmic and Physiological Optics 8.2 (1988):
        %	215-219.
        %
        % and
        %   Mutti, Donald O., Karla Zadnik, and Christopher J. Murphy.
        %   "Naturally occurring vitreous chamber-based myopia in the
        %   Labrador retriever." Investigative ophthalmology & visual
        %   science 40.7 (1999): 1577-1584.
        %
        % Values are given for an emmetropic canine eye.

        %% Cornea front surface
        % I cannot find a value for the asphericity, so am using the human
        % value
        corneaFrontR = 8.375;
        corneaFrontQ = -0.15;        
        a = corneaFrontR / ( corneaFrontQ + 1 );
        b = corneaFrontR * sqrt(1/(corneaFrontQ+1)) ;
        eye.cornea.front.radii(1) = a;
        eye.cornea.front.radii(2:3) = b;
        
        % We set the axial apex of the corneal front surface at position
        % [0, 0, 0]
        eye.cornea.front.center = [-eye.cornea.front.radii(1) 0 0];
        
        %% Cornea back surface
        % Asphericity is the human value.
        corneaBackR = 8;
        corneaBackQ = -0.275;
        
        % Compute the radii of the ellipsoid
        a = corneaBackR / ( corneaBackQ + 1 );
        b = corneaBackR * sqrt(1/(corneaBackQ+1)) ;
        eye.cornea.back.radii(1) = a;
        eye.cornea.back.radii(2:3) = b;
        
        % The thickness of the canine cornea is given as 0.587 mm by:
        %   Alario, Anthony F., and Christopher G. Pirie. "Central corneal
        %   thickness measurements in normal dogs: a comparison between
        %   ultrasound pachymetry and optical coherence tomography."
        %   Veterinary ophthalmology 17.3 (2014): 207-211.
        %
        % The center of the cornea circle for the back surface is
        % positioned to provide this thickness  between
        % the front and back surface of the cornea at the apex. 
        eye.cornea.back.center = [-0.587-eye.cornea.back.radii(1) 0 0];
        
        % Lacking information otherwise, the cornea is aligned with the
        % optical axis
        eye.cornea.axis = [0 0 0];

        
        %% Iris
        % An anterior chamber depth of 4.29 mm is given in Table 3 of:
        %
        %   Thomasy, Sara M., et al. "Species differences in the geometry
        %   of the anterior segment differentially affect anterior chamber
        %   cell scoring systems in laboratory animals." Journal of Ocular
        %   Pharmacology and Therapeutics 32.1 (2016): 28-37.
        %
        % Adding corneal thickness gives an iris center depth of -4.877.
        % Apparently the iris plane is tilted substantially in the dog, so
        % some estimate of this will be needed. Using human value for
        % thickness.
        eye.iris.thickness = 0.15;
        eye.iris.radius = 7;
       	eye.iris.center = [-4.877 0 0];

        
        %% Pupil
        % Centered on iris and optical axis.
        eye.pupil.center = [-4.877 0 0];
        
        % We assume that the canine actual pupil is circular
        eye.pupil.eccenParams = []; 
        eye.pupil.eccenFcnString = sprintf('@(x) 0'); 
        eye.pupil.thetas = [0  0];
        

        %% vitreous chamber
        eye.retina.radii = [ 8.25 8.25 8.25];
        
        % This is the human value; Need to do the computation for the dog.
        retinaApexDepth = 3.25;

        if isempty(p.Results.axialLength)
            eye.axialLength = retinaApexDepth + eye.vitreousChamber.radii(1)*2;
        else
            % If a specific axial length was passed (perhaps obtained by
            % measurement using the IOL Master apparatus), set the model
            % eye to have this length, and scale the other dimensions of
            % the vitreous chamber to maintain the specified ametropia. We
            % adjust the axial length for the component of the anterior
            % chamber that contibutes to length (vitreousChamberApexDepth)
            scaleFactor = (p.Results.axialLength - vitreousChamberApexDepth) / (eye.vitreousChamberRadii(1)*2);
            eye.vitreousChamber.radii = eye.vitreousChamber.radii .* scaleFactor;
            eye.axialLength = p.Results.axialLength;
        end

        % Set the depth of the center of the vitreous chamber
        eye.vitreousChamber.center = ...
            [(-4.2 - eye.vitreousChamber.radii(1)) 0 0];

        % I have not yet implemented fovea and optic disc positioning in
        % the canine eye, so set these to nan for now
        eye.vitreousChamber.fovea = [nan nan nan];
        eye.vitreousChamber.opticDisc = [nan nan nan];
        
        %% Lens
        % The lens parameters are included to support an illustration of a
        % complete eye model. The front and back surfaces of the lens are
        % modeled as hyperbolas. This simplified model does not model the
        % gradient in refractive index across the extent of the lens, and
        % therefore does not support ray tracing. All values taken from
        % Atchison 2006.
        % To convert R and Q to radii of a hyperbola:
        %   R = b^2/a
        %	Q = (a^2 / b^2) + 1
        % Therefore, given R and Q, we can obtain a and b, which correspond
        % to the radii of the ellipsoid model, with a corresponding to the
        % axial dimension, and b to the horizontal and verical dimensions.
        % Checking my algebra here:
        %{
            syms a b R Q
            eqn1 = R == a^2/b;
            eqn2 = Q == (a^2 / b^2) + 1;
            solution = solve([eqn1, eqn2]);
            solution.a
            solution.b
        %}
        eye.lens.front.R = 6.945;
        eye.lens.front.Q = -5;
        a = eye.lens.front.R * sqrt(abs( 1 / (eye.lens.front.Q - 1 ) )) * sign(eye.lens.front.Q);
        b = eye.lens.front.R / (eye.lens.front.Q - 1 );
        eye.lens.front.radii(1) = b;
        eye.lens.front.radii(2:3) = a;
        eye.lens.front.center = [eye.pupil.center(1)-eye.lens.front.radii(1) 0 0];
        
        eye.lens.back.R = -6.520;
        eye.lens.back.Q = -2;
        a = eye.lens.back.R * sqrt(abs( 1 / (eye.lens.back.Q - 1 ) )) * sign(eye.lens.back.Q);
        b = eye.lens.back.R / (eye.lens.back.Q - 1 );
        eye.lens.back.radii(1) = b;
        eye.lens.back.radii(2:3) = a;
        eye.lens.back.center = [eye.pupil.center(1)-3.6-eye.lens.back.radii(1) 0 0];
        
        % I specify the location of a single nodal point to support
        % calculation of the visual axis. The nodal point is placed at the
        % depth of the anterior principal point given by 
        eye.lens.nodalPoint = [-8.907 0 0];        
        
        
        %% Axes
        % For the canine eye, I treat the optical and visual axes as
        % aligned. I don't have any information regarding the position of
        % the optic disc.
        eye.axes.optical.degRetina = [0 0 0];
        eye.axes.visual.degRetina = [0 0 0];
        eye.axes.opticalDisc.degRetina = [0 0 0];
        

        %% Rotation centers
        % For lack of better information, these are placed on the optical
        % axis, 10 mm posterior to the corneal apex.
        eye.rotationCenters.azi = [-10 0 0];
        eye.rotationCenters.ele = [-10 0 0];
        eye.rotationCenters.tor = [0 0 0];

        
        %% Refractive indices
        % Using the human values for now
        % Obtain refractive index values for this spectral domain.
        eye.index.cornea = returnRefractiveIndex( 'cornea', p.Results.spectralDomain );
        eye.index.aqueous = returnRefractiveIndex( 'aqueous', p.Results.spectralDomain );
        eye.index.lens = returnRefractiveIndex( 'lens', p.Results.spectralDomain );

        
    otherwise
        error('Please specify a valid species for the eye model');
end



end % function


