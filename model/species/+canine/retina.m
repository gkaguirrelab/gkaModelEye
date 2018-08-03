function [ output_args ] = untitled( input_args )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
        
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

end

