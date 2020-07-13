function [srRays] = stingray(srModel, srInitialize, srArc, srControl, srRays, N)

% EDITED BY ZOE KRAUSS TO WRITE ALL VARIABLES TO ASCII FILES, THEN READ THEM MANUALLY INTO A UNIX COMMAND
% INSTEAD OF USING MATLAB MEX FILES


%STINGRAY - Calculates travel times and ray paths using graph theory.
%             (Stingray toolbox)
%
%   srRays = stingray(srModel,srInitialize,srArc,srControl) returns the
%   result in the first plane of a new structure srRays.
%
%   srRays = stingray(srModel,srInitialize,srArc,srControl,srRays,N)
%   returns the result in the N'th plane of the input structure srRays.
%
%   INPUT:
%
%           srModel:            Stingray structure
%           srInitialize:       Stingray structure
%           srArc:              Stingray structure
%           srControl:          Stingray structure
%
%   OPTIONAL INPUT:
%
%           srRays:             Stingray structure
%                               Stingray(N).model specifies P or S model.
%                               If 4 input arguments, N=1 and P is assumed.
%           N:                  N'th plane where results should be put.
%
%   OUTPUT:
%
%           srRays:             Stingray structure
%
%
%   Stingray is a matlab pre-processor for the stingray mexmaci file.
%   Fortran code is not double precision.  This function fills the single
%   precision variables, using data from the Stingray structures:
%
%   The %val construct is used inside stingray_gateway.
%   stingray_gateway operates directly on t, iprec and inSetS. It is not
%   recommended to call stingray_gateway from outside this m-file, in order
%   to avoid memory problems.
%
%  After calling stingray the srRays structure is filled out

%  Copyright 2010 Blue Tech Seismics, Inc.

%% Check number of input arguments

if nargin < 4, help(mfilename), error('Not enough arguments');end
% if nargin ==5, help(mfilename),error('Specify number plane for srRays');end
if nargin ==4
    N = 1; 
    PorS = 'P'; 
else
    % PorS = char(srRays(N).model);
    PorS = srRays.phase;
    N = 1;
end




%% Define variables to call stingray, put them into single columns


u            = srModel.(PorS).u(:);
ghead        = srModel.ghead(:);


    %%%%% Changed to account for anisotropy symmetry system (GMA, 2018).
    %%%%% Negative represents slow symmetry system (i.e., cracks), positive
    %%%%% represents fast symmetry system (i.e., olivine).
    % a_r          = srModel.anis_sym(:) .* srModel.(PorS).anis_fraction(:); 
a_r          = srModel.anis_sym(:) .* srModel.(PorS).anis_fraction(:); 

a_r          = a_r(:);
a_t          = srModel.(PorS).anis_theta(:);
a_p          = srModel.(PorS).anis_phi(:);


knx          = srModel.nx;
kny          = srModel.ny;
knodes       = srModel.nodes;


inSetS       = srInitialize.inSetS(:);
t            = srInitialize.time(:);
iprec        = srInitialize.iprec(:);


zhang        = srModel.elevation(:);

%  Convert flags to single precision

% turning tf_anisotropy into an integer so it can be saved to ascii
tf_anisotropy1     = srControl.tf_anisotropy;
if tf_anisotropy1==true
    tf_anisotropy=1;
else
    tf_anisotropy=0;
end

tf_line_integrate = srControl.tf_line_integrate(:);
if tf_line_integrate ==true
    tf_line_integ=1;
else
    tf_line_integ=0;
end
tf_line_integ = tf_line_integ(:);


arcList      = srArc.arcList(:);
arcHead      = [srArc.mx srArc.my srArc.mz srArc.nfs];
kmx          = srArc.mx;
kmy          = srArc.my;
kmz          = srArc.mz;
kmaxfs       = srArc.nfs;

% Turning the rest of the variables into single columns

% combining the simple integer values
k7 = [knx kny kmx kmy kmz kmaxfs knodes];

k7 = k7(:);
arcHead = arcHead(:);
tf_anisotropy = tf_anisotropy(:);


% Create and change to new, unique directory before saving files
% This is to avoid redundancy issues when running in parallel

[~,iDir] = min(srInitialize.time);
filename = ['.stingray_' int2str(iDir) '-' num2str(now,'%17.10f')];
unix(['mkdir ' pwd '/' filename])
eval(['cd ' pwd '/' filename])

% save all variables to ascii files

save .stingray_u u -ascii -single;
save .stingray_t t -ascii -single;
save .stingray_iprec iprec -ascii -single;
save .stingray_inSetS inSetS -ascii -single;
save .stingray_ghead ghead -ascii -single;
save .stingray_k7 k7 -ascii -single;
save .stingray_arcList arcList -ascii -single;
save .stingray_arcHead arcHead -ascii -single;
save .stingray_zhang zhang -ascii -single;
save .stingray_tf_line_integ tf_line_integ -ascii -single;
save .stingray_tf_anisotropy tf_anisotropy -ascii -single;
save .stingray_a_r a_r -ascii -single;
save .stingray_a_t a_t -ascii -single;
save .stingray_a_p a_p -ascii -single;



%%  Call stingray_src_unix program

tic
!/Users/zoekrauss/Stingray_Updated/source/stingray_src_unix
toc

%%  Fill srRays structure
%
%  Required fields:
%
%       srRays.ghead               (1:8)
%       srRays.time                (nx, ny, nz)
%       srRays.iprec               (nx, ny, nz)
%
%  Derived fields:
%
%       srRays.nx                  nodes in x-direction
%       srRays.ny                  nodes in y-direction
%       srRays.nz                  nodes in z-direction
%       srRays.gx                  node-spacing in x
%       srRays.gy                  node-spacing in y
%       srRays.gz                  node-spacing in z
%       srRays.nodes               total number of nodes
%       srRays.xg                  x-location of nodes
%       srRays.yg                  y-location of nodes
%       srRays.zg                  z-location of nodes
%       srRays.elevation           mesh of elevation at nodes
%       srRays.srGeometry          srGeometry holds origin and rotaiton
%       srRays.modelname           srModel.filename (velocity model)
%
%  srRays should be able to describe itself completely.

%%  fill srRays

% catch to make sure dimensions specified in the include file, stingray_basedims, were made large enough
load .stingray_catch -ascii;
if X == 0
    disp('ERROR: The dimensions in the include file were not large enough.');
    return
end



srRays(N).ghead      = srModel.ghead;
load .stingray_t -ascii;
srRays(N).time       = reshape(X,srModel.nx,srModel.ny,srModel.nz);
load .stingray_iprec -ascii;
srRays(N).iprec      = reshape(X,srModel.nx,srModel.ny,srModel.nz);

srRays(N).nx         = srModel.nx;
srRays(N).ny         = srModel.ny;
srRays(N).nz         = srModel.nz;
srRays(N).gx         = srModel.gx;
srRays(N).gy         = srModel.gy;
srRays(N).gz         = srModel.gz;
srRays(N).nodes      = srModel.nodes;
srRays(N).xg         = srModel.xg;
srRays(N).yg         = srModel.yg;
srRays(N).zg         = srModel.zg;
srRays(N).elevation  = srModel.elevation;
srRays(N).srGeometry = srModel.srGeometry;
srRays(N).modelname  = srModel.filename;
srRays(N).ghead      = srModel.ghead;
srRays(N).srControl  = srControl;
srRays(N).LON  = srModel.LON;
srRays(N).LAT  = srModel.LAT;

% Delete temporary unique directory and change back to previous directory

cd ..
unix(['rm -r '  pwd '/' filename])