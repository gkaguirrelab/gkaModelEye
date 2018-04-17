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
%   optic / pupil axis, with the apex of the cornea set as zero in depth.
%   The space has the dimensions [depth, horizontal, vertical]; negative
%   values of depth are towards the back of the eye. The model assumes the
%   optical and pupil axis of the eye are algined.
%
% Inputs:
%   none
%
% Optional key/value pairs:
%  'sphericalAmetropia'   - Scalar, in units of diopters. The
%                           dimensions of the posterior chamber of the eye
%                           (and to a lesser extent the curvature of the
%                           cornea) change with the observed refractive
%                           error of the subject. This value is the
%                           spherical refractive correction for the
%                           subject. A negative number is the correction
%                           that would be used for a myopic person.
%  'axialLength'          - Scalar. When set, this fixes the axial length
%                           of the eye to the passed value in millimeters.
%                           As the modeled anterior chamber depth is not
%                           variable, this change is enforced on the
%                           posterior chamber. The remaining dimensions of
%                           the posterior chamber are scaled to fit the
%                           proportions predicted by the Atchison model for
%                           the specified degree of ametropia.
%  'gammaAngle'           - 1x4 vector. This is the angle of the fixation
%                           axis in degrees w.r.t. to optical axis. The
%                           values are [azimuth, elevation]. An eyePose of:
%                             [-gamma(1), -gamma(2), 0, radius]
%                           aligns the visual axis of the eye with the
%                           optical axis of the camera.
%  'eyeLaterality'        - A text string that specifies which eye (left,
%                           right) to model. Allowed values (in any case)
%                           are {'left','right','L','R','OS','OD'}
%  'species'              - A text string that specifies the species to be
%                           modeled. Supported values (in any case) are
%                           {'human'}
%  'spectralDomain'       - String, options include {'vis','nir'}.
%                           This is the wavelength domain within which
%                           imaging is being performed. The refractive
%                           indices vary based upon this choice.
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
p.addParameter('gammaAngle',[],@(x)(isempty(x) || isnumeric(x)));
p.addParameter('eyeLaterality','Right',@ischar);
p.addParameter('species','Human',@ischar);
p.addParameter('spectralDomain','nir',@ischar);

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

% Switch parameters at the top level by species
switch p.Results.species

    %% Human eye
    case {'human','Human','HUMAN'}
                
        %% Cornea
        % We model the cornea as an ellipsoid, taking the "canonical
        % representation" parameters from Table 1 of Navarro 2006:
        %
        %   Navarro, Rafael, Luis González, and José L. Hernández. "Optics
        %   of the average normal cornea from general and canonical
        %   representations of its surface topography." JOSA A 23.2 (2006):
        %   219-232.
        %
        % Their dimensions [a,b,c] correspond to our [p2, p3, p1].
        %
        % The radius of curvature at the vertex of the cornea was found by
        % Atchison to vary as a function of spherical ametropia (Table 1):
        %
        %	Atchison, David A. "Optical models for human myopic eyes."
        %	Vision research 46.14 (2006): 2236-2250.
        %
        % Atchison provides parameters for a radially symmetric ellipsoid
        % in terms of the radius of curvature (R) at the vertex and its
        % asphericity (Q). R varies with spherical ametropia (D):
        %
        %   R = 7.77 + 0.022 * D
        %   Q = -0.15
        % 
        % Because the asphericity of the cornea did not change, the change
        % in R corresponds to an overall scaling of the ellipsoid in all
        % dimensions. We adjust the Navarro values to account for this
        % effect. R and Q are related to the radii of an ellipse along the
        % primary and secondy axes (a, b) by:
        %
        %   R = b^2/a
        %	Q = (b^2 / a^2) - 1
        %
        % when Q < 0. Therefore, given R and Q, we can obtain a and b,
        % which correspond to the radii of the ellipsoid model, with a
        % corresponding to the axial dimension, and b to the horizontal and
        % verical dimensions. Checking my algebra here:
        %{
            syms a b R Q
            eqn1 = R == b^2/a;
            eqn2 = Q == (b^2 / a^2) - 1;
            solution = solve([eqn1, eqn2]);
            solution.a
            solution.b
        %}        
        % We calculate the change in parameters of the Navarro model that
        % would be expected given the Atchison effect for ametropia. 
        %{
            R = @(D) 7.77 + 0.022 .* D;
            Q = -0.15;
            a = @(D) R(D) ./ (Q+1);
            b = @(D) R(D) .* sqrt(1./(Q+1));
            radiiAtchFront = @(D) [a(D) b(D) b(D)];
            % Show that the ametropia correction scales all radii equally
            radiiAtchFront(0)./radiiAtchFront(1)
            % Calculate the proportion change in radius 
            radiusScalerPerD = 1-a(1)/a(0);
            radiiNavFront = [14.26   10.43   10.27];
            radiiNavFrontCorrected = @(D) radiiNavFront.* (D.*radiusScalerPerD+1);
            % Report the ratio of the Atchison and Navarro axial radii
            % for the front surface of the corneal; we use this below.
            atchNavScaler = a(0) ./ radiiNavFront(1)
        %}
        % Atchison finds that the back surface of cornea does not vary by
        % ametropia. Navarro does not provide posterior cornea parameters.
        % Therefore, we scale the parameters provided by Atchison to relate
        % to the axial corneal radius specified by Navarro:
        %{
            R = 6.4;
            Q = -0.275;
            a = R ./ (Q+1);
            b = R .* sqrt(1./(Q+1));
            % Taken from the prior block of code 
            atchNavScaler = 0.6410;
            % Scale the overall back cornea ellipsoid to match Navarro
            radiiNavBack = [a b b]./atchNavScaler;
            % Now scale the relative horizontal and vertical axes so that
            % the relationship between the horizontal (and vertical) radii
            % and the axial radius is of the same proportion to the front
            % surface in the Atchison model
            radiiAtchFront0D = radiiAtchFront(0);
            frontHorizToAxAtch = radiiAtchFront0D(2)/radiiAtchFront0D(1);
            backHorizToAxAtch = b / a;
            radiiNavFront0D = radiiNavFrontCorrected(0);
            frontHorizToAxNav = radiiNavFront0D(2)/radiiNavFront0D(1);
            backHorizToAxNav = radiiNavBack(2)/radiiNavBack(1);
            targetBackHorizToAxNav = backHorizToAxAtch / frontHorizToAxAtch * frontHorizToAxNav;
            radiiNavBackCorrected = [a a*targetBackHorizToAxNav a*targetBackHorizToAxNav]./atchNavScaler
        %}
        eye.cornea.front.radii = [14.26    10.43    10.27] .* ...
            ((p.Results.sphericalAmetropia .* -0.0028)+1);
        eye.cornea.back.radii = [13.7716    9.3027    9.3027];
        
        % We set the center of the cornea front surface ellipsoid so that
        % the axial apex is at position [0, 0, 0]
        eye.cornea.front.center = [-eye.cornea.front.radii(1) 0 0];
                
        % The center of the back cornea ellipsoid is positioned so that
        % there is 0.55 mm of corneal thickness between the front and back
        % surface of the cornea at the apex, following Atchison 2006.
        eye.cornea.back.center = [-0.55-eye.cornea.back.radii(1) 0 0];
        
        % Navarro 2006 measured the angle of rotation of the axes of the
        % corneal ellipsoid relative to the keratometric axis, which is the
        % axis that connects a fixation point with the center of curvature
        % of the cornea. We convert those angles here to be relative to the
        % optic axis of the eye. To do so, we first assume that the
        % keratometric axis is equal to the fixation axis [MAY WANT TO ADD
        % A CONVERSION STEP FOR THIS]. Next, we add the Navarro
        % measurements to the gamma angle values that we have for the
        % model.
        %{
            % Navarro values for the displacement of the corneal axis from
            % keratometric axis for the right eye (in degrees)
            keratometricAxisWRTcornealAxis = [2.35 0.85 0.02];
            % assume that the fixation and keratometric axes are equal
            fixationAxisWRTcornealAxis = [2.35 0.85 0.02];
            % specify our gamma angles
            eye = modelEyeParameters();
            fixationAxisWRTopticalAxis = eye.gamma;
            % Now obtain the corneal axes relative to optical axis
            cornealAxisWRTopticalAxis = fixationAxisWRTopticalAxis - fixationAxisWRTcornealAxis            
        %}
        switch eyeLaterality
            case 'Right'
                eye.cornea.axis = [2.6539    1.3017   -0.0200];
            case 'Left'
                eye.cornea.axis = [-2.6539    1.3017   0.0200];
        end

        
        %% Pupil
        % We position the pupil plane at the depth of the anterior point of
        % the lens. The coordinate space of the model eye is defined w.r.t.
        % the center of the pupil, so the p2 and p3 values are zero
        eye.pupil.center = [-3.7 0 0];
        
        % The exit pupil of the eye is elliptical. Further, the
        % eccentricity and theta of the exit pupil ellipse changes with
        % pupil dilation:
        %
        %   Wyatt, Harry J. "The form of the human pupil." Vision Research
        %   35.14 (1995): 2021-2036.
        %
        % Wyatt reported the average ellipse parameters for the entrance
        % pupil (with the visual axis aligned with camera axis) under dim
        % and bright light conditions. We calculate the corresponding
        % parameters of the exit pupil on the optical axis. We then fit a
        % hyperbolic tangent (sigmoidal) function to the the eccentricity
        % of the exit pupil as a function of the exit pupil radius. The
        % theta values observed by Wyatt were close to vertically
        % orientated in the dark, and horizontally oriented in the light,
        % so we round to these values. When the exit pupil eccentricity is
        % below zero, the theta is set to zero (horizontal), and above zero
        % value it is set to pi/2 (vertical). In the forward model, we take
        % the absolute value of the eccentricity returned by the parameters
        % for the exit pupil eccentrivity.
        %{
            % Observed entrance pupil diameters reported in Wyatt 1995.
            entranceRadius = [3.09/2 4.93/2];
            % Wyatt reported an eccentricity of the pupil of 0.21 under
            % dark conditions. We find that using that value produces
            % model results that disagree with Malthur 2013. We have
            % adopted an upper value of 0.15 instead. We also use the 
            % convention of a negative eccentricity for a horizontal major
            % axis and a positive eccentricity for vertical.
            entranceEccen = [-0.12 0.15];
            % Prepare scene geometry and eye pose aligned with visual axis
            sceneGeometry = createSceneGeometry();
            % Fix the exit pupil eccentricity at 0
            sceneGeometry.eye.pupil.eccenFcnString = '@(x) 0';
            sceneGeometry.eye.pupil.thetas = [0, 0];
            % Obtain the pupil area in the image for each entrance radius
            % assuming no ray tracing
            sceneGeometry.refraction = [];
            pupilImage = pupilProjection_fwd([-sceneGeometry.eye.gamma(1), -sceneGeometry.eye.gamma(2), 0, entranceRadius(1)],sceneGeometry);
            exitArea(1) = pupilImage(3);
            pupilImage = pupilProjection_fwd([-sceneGeometry.eye.gamma(1), -sceneGeometry.eye.gamma(2), 0, entranceRadius(2)],sceneGeometry);
            exitArea(2) = pupilImage(3);
            % Add the ray tracing function to the sceneGeometry
            sceneGeometry = createSceneGeometry();
            % Search across exit pupil radii to find the values that match
            % the observed entrance areas.
            myPupilEllipse = @(radius) pupilProjection_fwd([-sceneGeometry.eye.gamma(1), -sceneGeometry.eye.gamma(2), 0, radius],sceneGeometry);
            myArea = @(ellipseParams) ellipseParams(3);
            myObj = @(radius) (myArea(myPupilEllipse(radius))-exitArea(1)).^2;
            exitRadius(1) = fminunc(myObj, entranceRadius(1));
            myObj = @(radius) (myArea(myPupilEllipse(radius))-exitArea(2)).^2;
            exitRadius(2) = fminunc(myObj, entranceRadius(2));
            % Now find the exit pupil eccentricity that produces the
            % observed entrance pupil eccentricity
            place = {'eye' 'pupil' 'eccenFcnString'};
            sceneGeometry.eye.pupil.thetas = [0, 0];
            mySceneGeom = @(eccen) setfield(sceneGeometry,place{:},['@(x) ' num2str(eccen)]);
            myPupilEllipse = @(eccen) pupilProjection_fwd([-sceneGeometry.eye.gamma(1), -sceneGeometry.eye.gamma(2), 0, exitRadius(1)],mySceneGeom(eccen));
            myEccen = @(ellipseParams) ellipseParams(4);
            myObj = @(eccen) 1e4*(myEccen(myPupilEllipse(eccen))-abs(entranceEccen(1))).^2;
            exitEccen(1) = -fminsearch(myObj, 0.1);
            sceneGeometry.eye.pupil.thetas = [pi/2, pi/2];
            mySceneGeom = @(eccen) setfield(sceneGeometry,place{:},['@(x) ' num2str(eccen)]);
            myPupilEllipse = @(eccen) pupilProjection_fwd([-sceneGeometry.eye.gamma(1), -sceneGeometry.eye.gamma(2), 0, exitRadius(2)],mySceneGeom(eccen));
            myEccen = @(ellipseParams) ellipseParams(4);
            myObj = @(eccen) 1e4*(myEccen(myPupilEllipse(eccen))-abs(entranceEccen(2))).^2;
            exitEccen(2) = fminsearch(myObj, 0.2);        
            % We then interpolate the observed values, assuming that the
            % observed values are close to asymptote
            exitRadiusInterp = [exitRadius(1)-.5 exitRadius(1) mean(exitRadius) exitRadius(2) exitRadius(2)+.5];
            exitEccenInterp = [exitEccen(1)/0.96 exitEccen(1) mean(exitEccen) exitEccen(2) exitEccen(2)/0.96];
            % Fit a hand-tuned sigmoidal function
            sigFit = @(scaleX, shiftY, scaleY, x) (tanh((x-mean(exitRadius)).*scaleX)+shiftY)*scaleY;
            fitEccen = fit(exitRadiusInterp',exitEccenInterp',sigFit);
            fprintf('eye.pupil.eccenParams = [-%4.3f %4.3f %4.3f %4.3f];\n',mean(exitRadius),fitEccen.scaleX,fitEccen.shiftY,fitEccen.scaleY);
            % Plot the fit
            figure
            plot(exitRadiusInterp,exitEccenInterp,'kx');
            hold on
            plot(0.5:.1:3,fitEccen(0.5:.1:3),'-r');        
        %}
        % Specify the params and equation that defines the exit pupil
        % ellipse. This can be invoked as a function using str2func.
        eye.pupil.eccenParams = [-1.766 4.716 0.290 0.087]; 
        eye.pupil.eccenFcnString = sprintf('@(x) (tanh((x+%f).*%f)+%f)*%f',eye.pupil.eccenParams(1),eye.pupil.eccenParams(2),eye.pupil.eccenParams(3),eye.pupil.eccenParams(4)); 

        % The theta values of the exit pupil ellipse for eccentricities
        % less than, and greater than, zero. We have structure here to add
        % a bit of tilt from vertical by laterality, but are not currently
        % using it.
        switch eyeLaterality
            case 'Right'
                eye.pupil.thetas = [0  pi/2];
            case 'Left'
                eye.pupil.thetas = [0  pi/2];
        end
        
        
        %% Iris
        % Define the iris radius. One study measured the horizontal visible
        % iris diameter (HVID) in 200 people, and found a mean of 11.8 with
        % a range of 10.2 - 13.0.
        %    PJ Caroline & MP Andrew. "The Effect of Corneal Diameter on
        %    Soft Lens Fitting, Part 1" Contact Lens Spectrum, Issue: April
        %    2002
        %    https://www.clspectrum.com/issues/2002/april-2002/contact-lens-case-reports
        %
        % Bernd Bruckner of the company Appenzeller Kontaktlinsen AG
        % supplied me with a tech report from his company (HVID & SimK
        % study) that measured HVID in 461 people. These data yield a mean
        % iris radius of 5.92 mm, 0.28 SD. The values from the histogram
        % are represented here, along with a Gaussian fit to the
        % distribution
        %{
            counts = [0 2 2 0 0 4 5 12 19 23 36 44 52 41 39 37 43 30 28 12 15 10 4 1 0 2 0];
            HVIDRadiusmm = (10.5:0.1:13.1)/2;
            hvidGaussFit = fit(HVIDRadiusmm', counts', 'gauss1');
            hvidRadiusMean = hvidGaussFit.b1;
            hvidRadiusSD =  hvidGaussFit.c1;
            figure
            plot(HVIDRadiusmm, hvidGaussFit(HVIDRadiusmm), '-r')
            hold on
            plot(HVIDRadiusmm, counts, '*k')
            xlabel('HVID radius in mm')
            ylabel('counts')
        %}
        % The HVID is the refracted iris size. We can use the forward model
        % to find the size of the true iris.
        %{
            sceneGeometry = createSceneGeometry();
            virtualImageFunc = compileVirtualImageFunc(sceneGeometry);
            sceneGeometry.virtualImageFunc = [];
            % Get the area in pixels of a "pupil" that is the same radius
            % as the HVID when there is no ray tracing
            hvidP=pupilProjection_fwd([0 0 0 hvidRadiusMean],sceneGeometry);
            % Restore ray tracing
            sceneGeometry.virtualImageFunc = virtualImageFunc;
            % Set up the objective function
            myArea = @(p) p(3);
            myObj = @(r) (hvidP(3) - myArea(pupilProjection_fwd([0 0 0 r],sceneGeometry)))^2;
            [r,pixelError] = fminsearch(myObj,5.5);
            fprintf('An unrefracted iris radius of %4.2f yields a refracted HVID of %4.2f \n',r,hvidRadiusMean)
        %}
        % We use this true iris size and then subject the iris perimeter
        % points to refraction
        eye.iris.radius = 5.55;
        
        % The iris is shifted slightly temporally and upward with respect
        % to the pupil center:
        %
        %   ...the typical entrance pupil is decentered
        %   approximately 0.15 mm nasally and 0.1 mm inferior to the
        %   geometric center of the visible iris circumference
        %
        % Bennett, Edward S., and Barry A. Weissman, eds. Clinical contact
        % lens practice. Lippincott Williams & Wilkins, 2005, p119
        %
        % We model an eye with zero iris angle, and thus set the depth of
        % the iris plane equal to the pupil plane.
        switch eyeLaterality
            case 'Right'
                eye.iris.center = [-3.7 -0.15 0.1];
            case 'Left'
                eye.iris.center = [-3.7 0.15 0.1];
        end
        
        
        %% Lens
        % Although the lens does not influence the pupil tracking, we
        % include it here to support an illustration of a complete eye
        % model. The front and back surfaces of the lens are modeled as
        % hyperbolas. This simplified model does not model the gradient in
        % refractive index across the extent of the lens, and therefore
        % does not support ray tracing. All values taken from Atchison
        % 2006.
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
        eye.lens.front.R = 11.48;
        eye.lens.front.Q = -5;
        a = eye.lens.front.R * sqrt(abs( 1 / (eye.lens.front.Q - 1 ) )) * sign(eye.lens.front.Q);
        b = eye.lens.front.R / (eye.lens.front.Q - 1 );
        eye.lens.front.radii(1) = b;
        eye.lens.front.radii(2:3) = a;
        eye.lens.front.center = [-3.7-eye.lens.front.radii(1) 0 0];
        
        eye.lens.back.R = -5.9;
        eye.lens.back.Q = -2;
        a = eye.lens.back.R * sqrt(abs( 1 / (eye.lens.back.Q - 1 ) )) * sign(eye.lens.back.Q);
        b = eye.lens.back.R / (eye.lens.back.Q - 1 );
        eye.lens.back.radii(1) = b;
        eye.lens.back.radii(2:3) = a;
        eye.lens.back.center = [-7.3-eye.lens.back.radii(1) 0 0];
        
        
        %% Posterior chamber
        % The posterior chamber of the eye is modeled as an ellipsoid.
        % Atchison 2006 provides radii of curvature and asphericity for the
        % posterior chamber as they vary by spherical ametropia. We perform
        % the calculations here and save only the corresponding radii.
        % Calculated using the formula for a positive Q value. We can
        % compare the posterior chamber radii calculated here to those
        % reported in Atchison 2005, and we find they are very similar:
        %
        %   Atchison, David A., et al. "Shape of the retinal surface in
        %   emmetropia and myopia." Investigative ophthalmology & visual
        %   science 46.8 (2005): 2698-2707.
        %{
            eye = modelEyeParameters();
            fprintf('Atchison emmetropic eye 2005, posterior chamber radii [axial x horizontal x vertical]:\n');
            fprintf('\t10.148\t11.455\t11.365\n');
            fprintf('The current model eye, emmetropic eye:\n');
            fprintf('\t%4.3f \t%4.3f \t%4.3f \n',eye.posteriorChamberRadii(1),eye.posteriorChamberRadii(2),eye.posteriorChamberRadii(3));
        %}
        Rzx = -12.91-0.094*p.Results.sphericalAmetropia;
        Rzy = -12.72+0.004*p.Results.sphericalAmetropia;
        Qzx = 0.27+0.026*p.Results.sphericalAmetropia;
        Qzy = 0.25+0.017*p.Results.sphericalAmetropia;
        eye.posteriorChamber.radii = [ -(Rzx/(Qzx+1)) -(Rzx*sqrt(1/(Qzx+1))) -(Rzy*sqrt(1/(Qzy+1)))];

        % Our model holds the depth of the anterior chamber constant.
        % Atchison found that anterior chamber depth does not vary with
        % spherical ametropia, although this is not a consistent finding:
        %
        %   Hosny, Mohamed, et al. "Relationship between anterior chamber
        %   depth, refractive state, corneal diameter, and axial length."
        %   Journal of Refractive Surgery 16.3 (2000): 336-340.
        %
        % To position the posterior chamber, we need to know the distance
        % between the apex of the anterior chamber and the apex of the
        % posterior chamber. The value for this fixed distance is derived
        % from the Atchison 2006 model eye. The total interior depth of the
        % eye from Table 1 for an emmetrope is:
        %
        %   0.55 + 3.15 + 1.44 + 2.16 + 16.28 = 23.58
        %
        % The axial length of the posterior chamber ellipsoid in an
        % emmetrope is:
        %
        %   10.1654 * 2 = 20.3308
        %
        % Therefore, the apex of the posterior ellipsoid is:
        %   23.5800 - 20.3308 = 3.2492 mm
        % behind the corneal apex.
        posteriorChamberApexDepth = 3.2492;
        
        % Compute and store axial length
        if isempty(p.Results.axialLength)
            eye.axialLength = posteriorChamberApexDepth + eye.posteriorChamber.radii(1)*2;
        else
            % If a specific axial length was passed (perhaps obtained by
            % measurement using the IOL Master apparatus), set the model
            % eye to have this length, and scale the other dimensions of
            % the posterior chamber to maintain the specified ametropia. We
            % adjust the axial length for the component of the anterior
            % chamber that contibutes to length (posteriorChamberApexDepth)
            %
            % GKA to follow up: Axial length is usually measured with the
            % IOL master along the visual (as opposed to optic or
            % pupillary) axis of the eye. May want to correct for this
            % somewhere.
            scaleFactor = (p.Results.axialLength - posteriorChamberApexDepth) / (eye.posteriorChamber.radii(1)*2);
            eye.posteriorChamber.radii = eye.posteriorChamber.radii .* scaleFactor;
            eye.axialLength = p.Results.axialLength;
        end
        
        % Set the depth of the center of the posterior chamber
        eye.posteriorChamber.center = ...
            [(-posteriorChamberApexDepth - eye.posteriorChamber.radii(1)) 0 0];
        

        %% Rotation centers
        % The rotation center of the eye is often treated as a single,
        % fixed point. A typical assumption is that the eye center of
        % rotation in emmetropes is 13.3 mm behind the corneal apex:
        %
        %   Gunter K. vonNoorden, MD; Emilio C. Campos "Binocular Vision
        %   and Ocular Motility Theory and Management of Strabismus"
        %   American Orthoptic Journal 51.1 (2001): 161-162.
        %
        % The source of this value in the cited text is not entirely clear.
        % It appears to be some compromise between the observed centers of
        % rotation that are obtained for azimuthal and elevation rotations. 
        % Measurements by Fry & Hill in 1962 and 1963 find that the
        % center of rotation is slightly nasal to the optical axis of the
        % eye, and differs for horizontal and vertical rotations:
        %
        %   Fry, G. A., and W. W. Hill. "The center of rotation of the
        %   eye." Optometry and Vision Science 39.11 (1962): 581-595.
        %
        %   Fry, Glenn A., and W. W. Hill. "The mechanics of elevating the
        %   eye." Optometry and Vision Science 40.12 (1963): 707-716.
        %
        % This difference in the apparent horizontal and vertical radii of
        % the eye was subsequently confirmed:
        %
        %   Hayami, Takehito, Kazunori Shidoji, and Katsuya Matsunaga. "An
        %   ellipsoidal trajectory model for measuring the line of sight."
        %   Vision research 42.19 (2002): 2287-2293.
        %
        % Fry & Hill report that the average azimuthal center of rotation
        % was 14.8 mm posterior to the corneal apex (14.7 in the
        % emmetropes), and 0.79 mm nasal to visual axis; and the elevation
        % center of rotation was 12.2 mm posterior to the corneal apex
        % (12.0 in the emmetropes) and 0.33 mm superior. These measurements
        % were made relative to the visual axis of the eye. While our model
        % is in optical axis coordinates, the effect of this difference is
        % very small (less than 1/100th of a millimeter).
        % 
        % Note that the Fry & Hill measurements supersede the earlier, Park
        % & Park measurements that claimed substantial translation of the
        % eye during rotation:
        %
        %   Park, Russell Smith, and George E. Park. "The center of ocular
        %   rotation in the horizontal plane." American Journal of
        %   Physiology--Legacy Content 104.3 (1933): 545-552.
        % 
        % The Park & Park result was due to their assumption that all
        % "sight lines" (i.e., rotations of the visual axis of the eye)
        % pass through the same point in space. Fry & Hill that some
        % subjects (2 of 31) show translation of the eye with rotation.
        % Also, there is a small, transient retraction of the eye following
        % a saccade that we do not attempt to model:
        %
        %   Enright, J. T. "The aftermath of horizontal saccades: saccadic
        %   retraction and cyclotorsion." Vision research 26.11 (1986):
        %   1807-1814.
        %
        % We provide three rotation centers, corresponding to the point of
        % rotation for azimuth, elevation, and torsional eye movements. The
        % values differ by eye because of the nasal displacement of the
        % rotation center.
        switch eyeLaterality
            case 'Right'
                eye.rotationCenters.azi = [-14.7 0.79 0];
            case 'Left'
                eye.rotationCenters.azi = [-14.7 -0.79 0];
        end
        eye.rotationCenters.ele = [-12.0 0 -0.33];
        eye.rotationCenters.tor = [0 0 0];
        
        % Spherical ametropia is correlated with the axial length of the
        % eye. We assume here that the center of rotation reflects this
        % change in length. Fry & Hill found that azimuthal rotation depth
        % increased by 0.167 mm for each negative diopter of spherical
        % refraction, and elevation rotation depth by 0.15 mm for each
        % negative diopter. Dick and colleagues (Figure 6) found that for
        % each mm of increase in axial length, the center of rotation
        % increased by 0.5 mm:
        %
        %   Dick, Graham L., Bryan T. Smith, and Peter L. Spanos.
        %   "Axial length and radius of rotation of the eye."
        %   Clinical and Experimental Optometry 73.2 (1990): 43-50.
        %
        % Given that in the Atchison data the axial length of the eye
        % increases by 0.27 mm for each negative diopter of spherical
        % ametropic error, this would imply a lengthening of the radius of
        % eye rotation by 0.14 mm, which is in good agreement with the Fry
        % & Hill observation of 0.15 - 0.167 mm of increase.
        %
        % We scale the azi and ele rotation centers by the ratio of the
        % posterior chamber axial and vertical radii relative to the
        % emmetropic size
        emmetropicPostChamberRadii = [10.1653543307087 11.455772536562 11.3771138695189];
        for dim=1:3
            eye.rotationCenters.azi(dim) = eye.rotationCenters.azi(dim) .* (eye.posteriorChamber.radii(dim)/emmetropicPostChamberRadii(dim));
            eye.rotationCenters.ele(dim) = eye.rotationCenters.ele(dim) .* (eye.posteriorChamber.radii(dim)/emmetropicPostChamberRadii(dim));
            eye.rotationCenters.tor(dim) = eye.rotationCenters.tor(dim) .* (eye.posteriorChamber.radii(dim)/emmetropicPostChamberRadii(dim));
        end

        
        %% Gamma
        % We now calculate gamma, which is the angle (in degrees) between
        % the optical and fixation axes of the eye. A related measurement
        % is kappa, which is the angle between the pupil and visual axes.
        % We use the names and greek letter designations for eye axes from
        % Atchison & Smith:
        %
        %   Atchison, David A., George Smith, and George Smith. "Optics of
        %   the human eye." (2000): 34-35.
        %
        % The optical and pupil axes of the our model eye are aligned. The
        % fixation axis of the eye is the line that connects the fixation
        % point to the center of rotation of the eye. The visual axis of
        % the eye connects the fixation point to the fovea via nodal points
        % through the eye optics. The visual axis lies within 1% of the
        % fixation axis when the fixation target is greater than 50 cm away
        % (Atchison & Smith, 2000, p.37). Because the visual and fixation
        % axes are closely aligned, and because our model equates the pupil
        % and optical axes, we treat empirical measurements of kappa
        % (visual to pupillary axis) as if they were the same as gamma
        % (fixation to optical axis). A related angle (equivalent to kappa
        % and gamma in our model) is alpha, which relates the visual and
        % optical axes.
        % 
        % The visual axis is displaced nasally and inferiorly within the
        % visual field relative to the optical axis. Horizontal kappa is
        % usually defined with positive values being more nasal. We adopt a
        % different convention in which kappa is defined in head-fixed
        % coordinates. Thus, positive values for the right eye, and
        % negative values for the left eye, are more nasal. Positive values
        % for vertical kappa are upward.
        %
        % A source for an estimate of kappa comes from Mathur 2013:
        %
        %	Mathur, Ankit, Julia Gehrmann, and David A. Atchison. "Pupil
        %	shape as viewed along the horizontal visual field." Journal of
        %	vision 13.6 (2013): 3-3.
        %
        % They measured the shape of the entrance pupil as a function of
        % viewing angle relative to the fixation point of the eye. Their
        % data from the right eye is well fit by a horizontal kappa of 5.5
        % degrees (see TEST_Mathur2013.m).
        %
        % While a horizontal kappa of ~5 degrees is a consistent finding,
        % measurements of vertical kappa differ:
        %
        %   Hashemi, Hassan, et al. "Distribution of angle kappa
        %   measurements with Orbscan II in a population-based survey."
        %   Journal of Refractive Surgery 26.12 (2010): 966-971.
        %
        %   Gharaee, Hamid, et al. "Angle kappa measurements: normal values
        %   in healthy iranian population obtained with the Orbscan II."
        %   Iranian Red Crescent Medical Journal 17.1 (2015).
        %
        % We note that there is evidence that the vertical kappa value can
        % vary based upon the subject being in a sittng or supine position.
        % Tscherning measured a downward alpha of 2-3 degrees, although
        % this varied amongst subjects:
        %
        %   Tscherning, Marius Hans Erik. Physiologic Optics: Dioptrics of
        %   the Eye, Functions of the Retina Ocular Movements and Binocular
        %   Vision. Keystone Publishing Company, 1920.
        %
        % Until better evidene is available, we adopt a vertical kappa of
        % -2.15 degrees for the emmetropic model eye, as this value best
        % fits the Mathur 2013 data.
        %
        % Measured kappa has been found to depend upon axial length:
        %
        %   Tabernero, Juan, et al. "Mechanism of compensation of
        %   aberrations in the human eye." JOSA A 24.10 (2007): 3274-3283.
        %
        % Tabernero 2007 report a mean horizontal kappa of 5 degrees in
        % emmetropes, and their Equation 6 expresses kappa (technically
        % alpha, the angle w.r.t. the optical axis) as a function of axial
        % length. Their formula assumes an emmetropic model eye of 24 mm,
        % while the model eye used here has an emmetropic axial length of
        % 23.592. The equation implemented below is adjusted so that an
        % emmetropic eye of 23.5924 mm has a horizontal (nasal directed)
        % kappa of 5.5 degrees and a vertical (inferiorly directed) kappa
        % of -2.15 degrees.
        %
        if isempty(p.Results.gammaAngle)
            switch eyeLaterality
                case 'Right'
                    eye.gamma(1) = atand((15.0924/(eye.axialLength-8.5000))*tand(5.5));
                case 'Left'
                    eye.gamma(1) = -atand((15.0924/(eye.axialLength-8.5000))*tand(5.5));
            end
            eye.gamma(2) = atand((15.0924/(eye.axialLength-8.5000))*tand(2.15));
            eye.gamma(3)=0;
        else
            eye.gamma = p.Results.gammaAngle;
        end
        
        
        %% Refractive indices
        % Obtain refractive index values for this spectral domain.
        eye.index.cornea = returnRefractiveIndex( 'cornea', p.Results.spectralDomain );
        eye.index.aqueous = returnRefractiveIndex( 'aqueous', p.Results.spectralDomain );
        eye.index.lens = returnRefractiveIndex( 'lens', p.Results.spectralDomain );

        
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
        eye.cornea.front.R = 8.375;
        eye.cornea.front.Q = -0.15;        
        a = eye.cornea.front.R / ( eye.cornea.front.Q + 1 );
        b = eye.cornea.front.R * sqrt(1/(eye.cornea.front.Q+1)) ;
        eye.cornea.front.radii(1) = a;
        eye.cornea.front.radii(2:3) = b;
        
        % We set the axial apex of the corneal front surface at position
        % [0, 0, 0]
        eye.cornea.front.center = [-eye.cornea.front.radii(1) 0 0];
        
        %% Cornea back surface
        % Asphericity is the human value.
        eye.cornea.back.R = 8;
        eye.cornea.back.Q = -0.275;
        
        % Compute the radii of the ellipsoid
        a = eye.cornea.back.R / ( eye.cornea.back.Q + 1 );
        b = eye.cornea.back.R * sqrt(1/(eye.cornea.back.Q+1)) ;
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
        

        %% Pupil
        % We position the pupil plane at the depth of the anterior point of
        % the lens. Table 3 of:
        %
        %   Thomasy, Sara M., et al. "Species differences in the geometry
        %   of the anterior segment differentially affect anterior chamber
        %   cell scoring systems in laboratory animals." Journal of Ocular
        %   Pharmacology and Therapeutics 32.1 (2016): 28-37.
        %
        % gives an anterior chamber depth of 4.29 mm. We must then add
        % corneal thickness to properly position the pupil plane.
        eye.pupil.center = [-4.877 0 0];
        
        % We assume that the canine exit pupil is circular
        eye.pupil.eccenParams = []; 
        eye.pupil.eccenFcnString = sprintf('@(x) 0'); 
        % The theta values of the exit pupil ellipse for eccentricities
        % less than and greater than zero.
        eye.pupil.thetas = [0  0];
        
        
        %% Iris
        % Need values for this. Apparently the iris plane is tilted
        % substantially in the dog, so some estimate of this will be
        % needed.
        eye.iris.radius = 7;
       	eye.iris.center = [-4.877 0 0];


        %% Posterior chamber
        eye.posteriorChamber.radii = [ 8.25 8.25 8.25];
        
        % This is the human value; Need to do the computation for the dog.
        posteriorChamberApexDepth = 3.25;

        if isempty(p.Results.axialLength)
            eye.axialLength = posteriorChamberApexDepth + eye.posteriorChamber.radii(1)*2;
        else
            % If a specific axial length was passed (perhaps obtained by
            % measurement using the IOL Master apparatus), set the model
            % eye to have this length, and scale the other dimensions of
            % the posterior chamber to maintain the specified ametropia. We
            % adjust the axial length for the component of the anterior
            % chamber that contibutes to length (posteriorChamberApexDepth)
            scaleFactor = (p.Results.axialLength - posteriorChamberApexDepth) / (eye.posteriorChamberRadii(1)*2);
            eye.posteriorChamber.radii = eye.posteriorChamber.radii .* scaleFactor;
            eye.axialLength = p.Results.axialLength;
        end
        
        % Set the depth of the center of the posterior chamber
        eye.posteriorChamber.center = ...
            [(-4.2 - eye.posteriorChamber.radii(1)) 0 0];
        
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

% Meta data regarding the units of the model
eye.meta.p = p.Results;
eye.meta.units = 'mm';
eye.meta.coordinates = 'eyeWorld';
eye.meta.dimensions = {'depth (axial)' 'horizontal' 'vertical'};
eye.meta.gamma = 'Degrees angle of fixation axis w.r.t. optical axis.';

end % function

