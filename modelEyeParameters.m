function eye = modelEyeParameters( varargin )
% Return the parameters of a model eye
%
% Syntax:
%  eye = modelEyeParameters()
%
% Description:
%   This routine returns the parameters of a model eye used in the
%   sceneGeometry routines. We make use of values derived by Atchison for
%   human eyes:
%
%       Atchison, David A. "Optical models for human myopic eyes." Vision
%       research 46.14 (2006): 2236-2250.
%
%       Atchison, David A., et al. "Shape of the retinal surface in
%       emmetropia and myopia." IOVS 46.8 (2005): 2698-2707.
%
%   Atchison uses the dimensions [x, y, z] corresponding to the width,
%   height, and depth (axial length) of the model eye. The parameters
%   returned by this routine correspond to the eyeWorld coordinate space
%   used in pupilProjection_fwd, which is relative to the pupil axis, with
%   the apex of the cornea set as zero in depth. The space has the
%   dimensions [depth, horizontal, vertical]; negative values of depth are
%   towards the center of the eye. The model assumes the optical and pupil
%   axis of the eye are algined.
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
%  'kappaAngle'           - 1x2 matrix. This is the angle of the visual
%                           axis in degrees w.r.t. to pupil axis. The
%                           values are [azimuth, elevation]. An eyePose of:
%                               [-kappa(1), -kappa(2), 0, radius]
%                           aligns the visual axis of the eye with the
%                           optical axis of the camera
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
p.addParameter('kappaAngle',[],@(x)(isempty(x) || isnumeric(x)));
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
                
        %% Cornea front surface
        % The properties of the cornea are typically described by the
        % radius of curvature (R) at the vertex and its asphericity (Q).
        % These values are taken from Atchison 2006, Table 1. The radius of
        % curvature of the front surface at the apex varies by spherical
        % ametropia of the eye; Q does not vary.
        eye.corneaFrontSurfaceR = 7.77 + 0.022 * p.Results.sphericalAmetropia;
        eye.corneaFrontSurfaceQ = -0.15;
        
        % The cornea is modeled as a prolate ellipsoid that is radially
        % symmetric about the optical axis. We calculate here the radii of
        % the ellipsoid. The radii of an ellipse along the primary and
        % secondy axes (a, b) are related to R and Q by:
        %   R = b^2/a
        %	Q = (b^2 / a^2) - 1
        % when Q < 0. Therefore, given R and Q, we can obtain a and b,
        % which correspond to the radii of the ellipsoid model, with a
        % corresponding to the axial dimension, and b to the horizontal and
        % verical dimensions.
        % Checking my algebra here:
        %{
            syms a b R Q
            eqn1 = R == b^2/a;
            eqn2 = Q == (b^2 / a^2) - 1;
            solution = solve([eqn1, eqn2]);
            solution.a
            solution.b
        %}        
        a = eye.corneaFrontSurfaceR / ( eye.corneaFrontSurfaceQ + 1 );
        b = eye.corneaFrontSurfaceR * sqrt(1/(eye.corneaFrontSurfaceQ+1)) ;
        eye.corneaFrontSurfaceRadii(1) = a;
        eye.corneaFrontSurfaceRadii(2:3) = b;
        
        % We set the axial apex of the corneal front surface at position
        % [0, 0, 0]
        eye.corneaFrontSurfaceCenter = [-eye.corneaFrontSurfaceRadii(1) 0 0];
        
        
        %% Cornea back surface
        % The radius of curvature for the back corneal surface was not
        % found to vary by spherical ametropia. The asphericity Q for the
        % back corneal surface was set by Atchison to -0.275.
        eye.corneaBackSurfaceR = 6.4;
        eye.corneaBackSurfaceQ = -0.275;
        
        % Compute the radii of the ellipsoid
        a = eye.corneaBackSurfaceR / ( eye.corneaBackSurfaceQ + 1 );
        b = eye.corneaBackSurfaceR * sqrt(1/(eye.corneaBackSurfaceQ+1)) ;
        eye.corneaBackSurfaceRadii(1) = a;
        eye.corneaBackSurfaceRadii(2:3) = b;
        
        % The center of the cornea circle for the back surface is
        % positioned so that there is 0.55 mm of corneal thickness between
        % the front and back surface of the cornea at the apex, following
        % Atchison 2006.
        eye.corneaBackSurfaceCenter = [-0.55-eye.corneaBackSurfaceRadii(1) 0 0];
        
        
        %% Pupil
        % We position the pupil plane at the depth of the anterior point of
        % the lens. The coordinate space of the model eye is define with
        % respect to the center of the pupil, so the p2 and p3 values are
        % zero
        eye.pupilCenter = [-3.7 0 0];
        
        
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
        eye.irisRadius = 5.92;
        
        % We align the iris center withi the optical axis, although we note
        % that there are some reporst that the iris is shifted slightly
        % temporally and upward with respect to the pupil center:
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
                eye.irisCenter = [-3.7 0 0];
            case 'Left'
                eye.irisCenter = [-3.7 0 0];
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
        eye.lensFrontSurfaceR = 11.48;
        eye.lensFrontSurfaceQ = -5;
        a = eye.lensFrontSurfaceR * sqrt(abs( 1 / (eye.lensFrontSurfaceQ - 1 ) )) * sign(eye.lensFrontSurfaceQ);
        b = eye.lensFrontSurfaceR / (eye.lensFrontSurfaceQ - 1 );
        eye.lensFrontSurfaceRadii(1) = b;
        eye.lensFrontSurfaceRadii(2:3) = a;
        eye.lensFrontSurfaceCenter = [-3.7-eye.lensFrontSurfaceRadii(1) 0 0];
        
        eye.lensBackSurfaceR = -5.9;
        eye.lensBackSurfaceQ = -2;
        a = eye.lensBackSurfaceR * sqrt(abs( 1 / (eye.lensBackSurfaceQ - 1 ) )) * sign(eye.lensBackSurfaceQ);
        b = eye.lensBackSurfaceR / (eye.lensBackSurfaceQ - 1 );
        eye.lensBackSurfaceRadii(1) = b;
        eye.lensBackSurfaceRadii(2:3) = a;
        eye.lensBackSurfaceCenter = [-7.3-eye.lensBackSurfaceRadii(1) 0 0];
        
        
        %% Posterior chamber
        % The posterior chamber of the eye is modeled as an ellipsoid.
        % Atchison 2006 provides radii of curvature and asphericity for the
        % posterior chamber as they vary by spherical ametropia. We perform
        % the calculations here and save only the corresponding radii.
        % Calculated using the formula for a positive Q value. We can
        % compare the posterior chamber radii calculated here to those
        % reported in Atchison 2005, and we find they are very similar.
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
        eye.posteriorChamberRadii = [ -(Rzx/(Qzx+1)) -(Rzx*sqrt(1/(Qzx+1))) -(Rzy*sqrt(1/(Qzy+1)))];

        % The model holds the depth of the anterior chamber constant. To
        % position the posterior chamber, we need to know the distance
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
        %   10.165 * 2 = 20.33
        %
        % Therefore, the apex of the posterior ellipsoid is:
        %   23.5800 - 20.33 = 3.2964 mm
        % behind the corneal apex.
        posteriorChamberApexDepth = 3.25;
        
        % Compute and store axial length
        if isempty(p.Results.axialLength)
            eye.axialLength = posteriorChamberApexDepth + eye.posteriorChamberRadii(1)*2;
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
            scaleFactor = (p.Results.axialLength - posteriorChamberApexDepth) / (eye.posteriorChamberRadii(1)*2);
            eye.posteriorChamberRadii = eye.posteriorChamberRadii .* scaleFactor;
            eye.axialLength = p.Results.axialLength;
        end
        
        % Set the depth of the center of the posterior chamber
        eye.posteriorChamberCenter = ...
            [(-posteriorChamberApexDepth - eye.posteriorChamberRadii(1)) 0 0];
        

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
        % Note that the Fry & Hill measurements superseed the earlier, Park
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
            eye.rotationCenters.azi(dim) = eye.rotationCenters.azi(dim) .* (eye.posteriorChamberRadii(dim)/emmetropicPostChamberRadii(dim));
            eye.rotationCenters.ele(dim) = eye.rotationCenters.ele(dim) .* (eye.posteriorChamberRadii(dim)/emmetropicPostChamberRadii(dim));
            eye.rotationCenters.tor(dim) = eye.rotationCenters.tor(dim) .* (eye.posteriorChamberRadii(dim)/emmetropicPostChamberRadii(dim));
        end
        
        %% Kappa
        % We now calculate kappa, which is the angle (in degrees) between
        % the pupil and visual axes of the eye. The visual axis is
        % displaced nasally and superiorly within the visual field relative
        % to the pupil axis. Horizontal kappa is usually defined with
        % positive values being more nasal. We adopt a different convention
        % in which kappa is defined in head-fixed coordinates. Thus,
        % positive values for the right eye, and negative values for the
        % left eye, are more nasal. Positive values for vertical kappa are
        % upward.
        %
        % A source for an estimate of kappa comes from Mathur 2013:
        %
        %	Mathur, Ankit, Julia Gehrmann, and David A. Atchison. "Pupil
        %	shape as viewed along the horizontal visual field." Journal of
        %	vision 13.6 (2013): 3-3.
        %
        % They measured the shape of the entrance pupil as a function of
        % viewing angle relative to the fixation point of the eye. Their
        % data is well fit by a kappa of [5, -2] degrees (see
        % TEST_entrancePupilShape.m).
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
        % kappa of 5 degrees and a vertical (inferiorly directed) kappa of
        % -2 degrees.
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
        % Until better evidene is available, we adopt a vertical kappa of
        % -2 degrees for the emmetropic model eye.        
        if isempty(p.Results.kappaAngle)
            switch eyeLaterality
                case 'Right'
                    eye.kappaAngle(1) = atand((15.0924/(eye.axialLength-8.5000))*tand(5));
                case 'Left'
                    eye.kappaAngle(1) = -atand((15.0924/(eye.axialLength-8.5000))*tand(5));
            end
            eye.kappaAngle(2) = atand((15.0924/(eye.axialLength-8.5000))*tand(2));
        else
            eye.kappaAngle = p.Results.kappaAngle;
        end
        
        
        %% Refractive indices
        % Obtain refractive index values for this spectral domain.
        eye.corneaRefractiveIndex = returnRefractiveIndex( 'cornea', p.Results.spectralDomain );
        eye.aqueousRefractiveIndex = returnRefractiveIndex( 'aqueous', p.Results.spectralDomain );
        eye.lensRefractiveIndex = returnRefractiveIndex( 'lens', p.Results.spectralDomain );

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
        eye.corneaFrontSurfaceR = 8.375;
        eye.corneaFrontSurfaceQ = -0.15;        
        a = eye.corneaFrontSurfaceR / ( eye.corneaFrontSurfaceQ + 1 );
        b = eye.corneaFrontSurfaceR * sqrt(1/(eye.corneaFrontSurfaceQ+1)) ;
        eye.corneaFrontSurfaceRadii(1) = a;
        eye.corneaFrontSurfaceRadii(2:3) = b;
        
        % We set the axial apex of the corneal front surface at position
        % [0, 0, 0]
        eye.corneaFrontSurfaceCenter = [-eye.corneaFrontSurfaceRadii(1) 0 0];
        
        %% Cornea back surface
        % Asphericity is the human value.
        eye.corneaBackSurfaceR = 8;
        eye.corneaBackSurfaceQ = -0.275;
        
        % Compute the radii of the ellipsoid
        a = eye.corneaBackSurfaceR / ( eye.corneaBackSurfaceQ + 1 );
        b = eye.corneaBackSurfaceR * sqrt(1/(eye.corneaBackSurfaceQ+1)) ;
        eye.corneaBackSurfaceRadii(1) = a;
        eye.corneaBackSurfaceRadii(2:3) = b;
        
        % The thickness of the canine cornea is given as 0.587 mm by:
        %   Alario, Anthony F., and Christopher G. Pirie. "Central corneal
        %   thickness measurements in normal dogs: a comparison between
        %   ultrasound pachymetry and optical coherence tomography."
        %   Veterinary ophthalmology 17.3 (2014): 207-211.
        %
        % The center of the cornea circle for the back surface is
        % positioned to provide this thickness  between
        % the front and back surface of the cornea at the apex. 
        eye.corneaBackSurfaceCenter = [-0.587-eye.corneaBackSurfaceRadii(1) 0 0];
        

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
        eye.pupilCenter = [-4.877 0 0];
        
        
        %% Iris
        % Need values for this. Apparently the iris plane is tilted
        % substantially in the dog, so some estimate of this will be
        % needed.
        eye.irisRadius = 7;
       	eye.irisCenter = [-4.877 0 0];


        %% Posterior chamber
        eye.posteriorChamberRadii = [ 8.25 8.25 8.25];
        eye.axialLength = p.Results.axialLength;
        
        % Set the depth of the center of the posterior chamber
        eye.posteriorChamberCenter = ...
            [(-4.2 - eye.posteriorChamberRadii(1)) 0 0];
        
        eye.rotationCenters.azi = [-10 0 0];
        eye.rotationCenters.ele = [-10 0 0];
        eye.rotationCenters.tor = [0 0 0];
        eye.corneaRefractiveIndex = returnRefractiveIndex( 'cornea', p.Results.spectralDomain );
        eye.aqueousRefractiveIndex = returnRefractiveIndex( 'aqueous', p.Results.spectralDomain );
        eye.lensRefractiveIndex = returnRefractiveIndex( 'lens', p.Results.spectralDomain );

        
    otherwise
        error('Please specify a valid species for the eye model');
end

% Meta data regarding the units of the model
eye.meta = p.Results;
eye.meta.units = 'mm';
eye.meta.coordinates = 'eyeWorld';
eye.meta.dimensions = {'depth (axial)' 'horizontal' 'vertical'};
eye.meta.kappa = 'Degrees angle of visual axis w.r.t. pupil axis.';

end % function

