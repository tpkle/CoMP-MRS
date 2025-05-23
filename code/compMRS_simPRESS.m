% compMRS_simPRESS.m
% Dana Goerzen and Jamie Near, McGill University 2021.
% Revised 2025 by Jamie Near & Diana Rotaru
% 
% USAGE:
% out=compMRS_simPRESS(spinSys,params);
% 
% DESCRIPTION:
% This script simulates a PRESS experiment with fully shaped refocusing 
% pulses. Coherence order filtering is employed to only simulate desired signals
% This results in a 4x speed up compared to phase cycling (see deprecated run_simPressShaped_phCyc.m)
% Furthermore, simulations are run at various locations in space to account for the 
% within-voxel spatial variation of the metabolite signal.  Summation 
% across spatial positions is performed. The MATLAB parallel computing toolbox 
% (parfor loop) was used to accelerate the simulations.  Acceleration 
% is currently performed in the direction of the slice selective pulse along
% the x-direction, but this can be changed.  Up to a factor of 12 acceleration
% can be achieved using this approach. To achieve 
% faster perfomance compared to the original 'run_simPressShaped_coFilt.m' function,
% this code uses the method described by Yan Zhang et al. Med Phys 2017;44(8): 
% 4169-78.  Some additional acceleration is currently performed using parfor 
% loops in both x and y directions.  To enable the use of the MATLAB
% parallel computing toolbox, initialize the multiple worked nodes using
% "matlabpool size X" where "X" is the number of available processing
% nodes.  If the parallel processing toolbox is not available, then replace
% the "parfor" loop with a "for" loop.
% 
% INPUTS:
% To run this script, there is technically only one input argument:
% spinSys           = spin system to simulate 
%
% However, the user should also edit the following parameters as 
% desired before running the function:
% refocWaveform     = name of refocusing pulse waveform.
% refTp             = duration of refocusing pulses[ms]
% Bfield            = Magnetic field strength in [T]
% Npts              = number of spectral points
% sw                = spectral width [Hz]
% Bfield            = magnetic field strength [Tesla]
% lw                = linewidth of the output spectrum [Hz]
% thkX              = slice thickness of x refocusing pulse [cm]
% thkY              = slice thickness of y refocusing pulse [cm]
% fovX              = full simulation FOV in the x direction [cm]
% fovY              = full simulation FOV in the y direction [cm]
% nX                = number of spatial grid points to simulate in x-direction
% nY                = number of spatial grid points to simulate in y-direction
% taus              = vector of pulse sequence timings  [ms]
%
% OUTPUTS:
% out               = Simulation results, summed over all space.

function out=compMRS_simPRESS(spinSys,params)
tic
% ************INPUT PARAMETERS**********************************
refocWaveform=params.refocWaveform; %name of refocusing pulse waveform.
refTp=params.refTp; %duration of refocusing pulses[ms]
flipAngle=180;  %Flip Angle of the refocusing pulses [degrees] (e.g. Use 180 for Siemens PRESS.  Use 137 for GE PRESS).
Npts=4096; %number of spectral points
sw=4000; %spectral width [Hz]
Bfield=params.B0; %magnetic field strength [Tesla]
lw=2; %linewidth of the output spectrum [Hz]
thkX=2; %slice thickness of x refocusing pulse [cm]
thkY=2; %slice thickness of y refocusing pulse [cm]
fovX=4; %size of the full simulation Field of View in the x-direction [cm]
fovY=4; %size of the full simulation Field of View in the y-direction [cm]
nX=48; %Number of grid points to simulate in the x-direction
nY=48; %Number of grid points to simulate in the y-direction
tau1=params.te1; %TE1 for first spin echo [ms]
tau2=params.te2; %TE2 for second spin echo [ms]
centreFreq=2.3; %Centre frequency of simulation [ppm]
% ************END OF INPUT PARAMETERS**********************************

%set up spatial grid
x=linspace(-fovX/2,fovX/2,nX); %X positions to simulate [cm]
y=linspace(-fovY/2,fovY/2,nY); %y positions to simulate [cm]

%Load RF waveform
refRF=io_loadRFwaveform(refocWaveform,'ref',0);

gamma=42577000; %gyromagnetic ratio

%Load spin systems
%load spinSystems
%sys=eval(['sys' spinSys]);

%Resample refocusing RF pulse from 400 pts to 100 pts to reduce
%computational workload
refRF=rf_resample(refRF,100);

Gx=(refRF.tbw/(refTp/1000))/(gamma*thkX/10000); %[G/cm]
Gy=(refRF.tbw/(refTp/1000))/(gamma*thkY/10000); %[G/cm]

%Initialize structures:
d_temp=cell(1,1);
d=cell(1,1);

%loop through space: If you are using the parfor loops below, and you are 
%using an older version of MATLAB (e.g.R2012), don't forget to initialize 
%the parallel processing toolbox workers using 'matlabpool open N' (for N 
%workers, 12 max).  I don't think this is necessary for newer version of 
%MATLAB.  

%First loop through all x-positions, simulating only the first refocusing
%pulse.  
%First loop through x direction (first refoc pulse only);

%for X=1:length(x)  %Use this if you don't have the MATLAB parallel processing toolbox
parfor X=1:length(x)  %Use this if you have the MATLAB parallel processing toolbox
        disp(['Executing X-position ' num2str(X) ' of ' num2str(length(x)) '!!!']);
        d_temp{X}=sim_press_shaped_fastRef1(Bfield,sys,tau1,tau2,refRF,refTp,x(X),Gx,flipAngle,centreFreq);
end

%calculate the average density matrix (Doing this inside a separate for 
%loop because I couldn't figure out how to do this inside the parfor loop): 
for X=1:length(x)
        d{1}=sim_dAdd(d{1},d_temp{X});
end


% %Initialize structures:
out_temp=cell(length(y),1);
out=struct([]);

%Now loop through y direction (second refoc pulse only);
%for Y=1:length(y) %Use this if you don't have the MATLAB parallel processing toolbox
parfor Y=1:length(y) %Use this if you do have the MATLAB parallel processing toolbox
%            disp(['Executing Y-position ' num2str(Y) ' of ' num2str(length(y)) '!!!']);
            out_temp{Y}=sim_press_shaped_fastRef2(d{1},Npts,sw,Bfield,lw,sys,tau1,tau2,...
                refRF,refTp,y(Y),Gy,flipAngle,centreFreq);
end

%Now combine the outputs;  Again, doing this inside a separate for loop
%becuase I can't figure out how to do this inside the parfor loop:
for Y=1:length(y) 
            out=op_addScans(out,out_temp{Y});
end

%For consistent scaling across different shaped simulations, we need to :
%1.  Scale down by the total number of simulations run (since these were
%    all added together.
numSims=(nX*nY);
out=op_ampScale(out,1/numSims);

%2.  Scale by the total size of the simulated region, relative to the size
%    of the voxel.
voxRatio=(thkX*thkY)/(fovX*fovY);
out=op_ampScale(out,1/voxRatio);

toc
end








%Nested Function #1
function d = sim_press_shaped_fastRef1(Bfield,sys,tau1,tau2,RF,tp,dx,Gx,flipAngle,centreFreq)
% 
% USAGE:
% d = sim_press_shaped_fastRef1(n,sw,Bfield,linewidth,sys,tau1,tau2,RF,tp,dx,Gx,phCyc1,flipAngle)
% 
% DESCRIPTION:
% This function simulates only the first bit of the PRESS experiment, up to 
% the beginning of the second refocusing pulse.  The excitation is
% simulated as an instantaneous rotation, and the refocusing pulse is
% simulated as a shaped rotation.
%
% This code is designed to be used in highly-accelerated shaped simulations,
% using the method described by Yan Zhang et al. Med Phys 2017;44(8): 
% 4169-78.
% 
% This code employs coherence order filtering to select only desired
% signal. See run_simPressShaped_coFilt description for more details.
% Finally, this code simulates the spectrum at a given point in space (x),
% given the values of the slice selection gradient (Gx).  In order
% to fully simulate the MEGA-PRESS experiment, you have to run this
% simulation many times at various points in space (x), followed by 
% sim_press_shaped_fastRef2.m, at all points in space (y).  
% 
% INPUTS:
% n         = number of points in fid/spectrum
% sw        = desired spectral width in [Hz]
% Bfield    = main magnetic field strength in [T]
% linewidth = linewidth in [Hz]
% sys       = spin system definition structure
% tau1      = echo time 1 in [ms].
% tau2      = echo time 2 in [ms].
% RF        = RF pulse definition structure for refoc pulses (obtain using 'io_loadRFwaveform.m')
% tp        = RF pulse duration in [ms]
% dx        = position offset in x-direction (corresponding to first refocusing pulse) [cm]
% dy        = position offset in y-direction (corresponding to second refocusing pulse) [cm]
% Gx        = gradient strength for first selective refocusing pulse [G/cm]
% Gy        = gradient strength for second selective refocusing pulse [G/cm]
% flipAngle = flip angle of refocusing pulses [degrees] (Optional.  Default = 180 deg)
%
% OUTPUTS:
% out       = simulated spectrum, in FID-A structure format, using PRESS 
%             sequence.
if nargin<10
    centreFreq=2.3;
    if nargin<9
        flipAngle=180;
    end
end
    
if tau1<tp/1000
    error('ERROR:  Echo-time 1 cannot be less than duration of refocusing pulse! ABORTING!!');
end
if tau2<tp/1000
    error('ERROR:  Echo-time 2 cannot be less than duration of refocusing pulse! ABORTING!!');
end

%Set centre frequency
for k=1:length(sys)
    sys(k).shifts=sys(k).shifts-centreFreq;
end

%Calculate Hamiltonian matrices and starting density matrix.
[H,d]=sim_Hamiltonian(sys,Bfield);

%Calculate new delays by subtracting the pulse duration from tau1 and tau2;
delays=zeros(2);
delays(1)=tau1-tp;
delays(2)=tau2-tp;
if sum(delays<0)
    error(['ERROR! The following taus are too short: ' num2str(find(delays<0)) '.']);
end

%BEGIN PULSE SEQUENCE************
d=sim_excite(d,H,'x');                                    %EXCITE
d=sim_COF(H,d,-1);
d=sim_evolve(d,H,delays(1)/2000);                            %Evolve by delays(1)/2
d=sim_shapedRF(d,H,RF,tp,flipAngle,90,dx,Gx);          %1st shaped 180 degree refocusing pulse
d=sim_COF(H,d,1);
d=sim_evolve(d,H,(delays(1)+delays(2))/2000);                     %Evolve by (delays(1)+delays(2))/2
%END PULSE SEQUENCE**************

%After running this many times along x, the density matrices should be
%averaged, and then the average density matrix should be passed through
%'sim_press_shaped_fastRef2' at various different y-positions. 


end









%Nested Function #2
function out = sim_press_shaped_fastRef2(d,n,sw,Bfield,linewidth,sys,tau1,tau2,RF,tp,dy,Gy,flipAngle,centreFreq)
%
% USAGE:
% out = sim_press_shaped_fastRef2(d,n,sw,Bfield,linewidth,sys,tau2,RF,tp,dy,Gy,phCyc2,flipAngle)
% 
% DESCRIPTION:
% This function simulates only the last bit of the PRESS experiment, from the 
% the beginning of the second refocusing pulse, to the end.  The refocusing 
%pulse is simulated as a shaped rotation.
%
% This code is designed to be used in highly-accelerated shaped simulations,
% using the method described by Yan Zhang et al. Med Phys 2017;44(8): 
% 4169-78.
 
% This code employs coherence order filtering to select only desired
% signal. See run_simPressShaped_coFilt description for more details
% 
% Finally, this code simulates the spectrum at a given point in space (y),
% given the values of the slice selection gradient (Gy).  In order
% to fully simulate the MEGA-PRESS experiment, you have to first run
% sim_press_shaped_fastRef1.m at all points in space (x), followed by 
% this code, at all points in space (y).  
% 
% INPUTS:
% d         = starting density matrix (obtained using 'sim_press_shaped_fastRef1.m')
% n         = number of points in fid/spectrum
% sw        = desired spectral width in [Hz]
% Bfield    = main magnetic field strength in [T]
% linewidth = linewidth in [Hz]
% sys       = spin system definition structure
% tau1      = echo time 1 in [ms].
% tau2      = echo time 2 in [ms].
% RF        = RF pulse definition structure for refoc pulses (obtain using 'io_loadRFwaveform.m')
% tp        = RF pulse duration in [ms]
% dx        = position offset in x-direction (corresponding to first refocusing pulse) [cm]
% dy        = position offset in y-direction (corresponding to second refocusing pulse) [cm]
% Gx        = gradient strength for first selective refocusing pulse [G/cm]
% Gy        = gradient strength for second selective refocusing pulse [G/cm]
% flipAngle = flip angle of refocusing pulses [degrees] (Optional.  Default = 180 deg)
%
% OUTPUTS:
% out       = simulated spectrum, in FID-A structure format, using PRESS 
%             sequence.

if nargin<14
    centreFreq=2.3;
    if nargin<13
        flipAngle=180;
    end
end
   
if tau1<tp/1000
    error('ERROR:  Echo-time 1 cannot be less than duration of refocusing pulse! ABORTING!!');
end
if tau2<tp/1000
    error('ERROR:  Echo-time 2 cannot be less than duration of refocusing pulse! ABORTING!!');
end

%Set centre Frequency
for k=1:length(sys)
    sys(k).shifts=sys(k).shifts-centreFreq;
end

%Calculate Hamiltonian matrices and starting density matrix.
[H]=sim_Hamiltonian(sys,Bfield);

%Calculate new delays by subtracting the pulse duration from tau1 and tau2;
delays=zeros(2);
delays(1)=tau1-tp;
delays(2)=tau2-tp;
if sum(delays<0)
    error(['ERROR! The following taus are too short: ' num2str(find(delays<0)) '.']);
end

%BEGIN PULSE SEQUENCE************
d=sim_shapedRF(d,H,RF,tp,flipAngle,90,dy,Gy);          %2nd shaped 180 degree refocusing pulse
d=sim_COF(H,d,-1);
d=sim_evolve(d,H,delays(2)/2000);                            %Evolve by delays(2)/2
[out,~]=sim_readout(d,H,n,sw,linewidth,90);      %Readout along y (90 degree phase);
%END PULSE SEQUENCE**************

%Correct the ppm scale:
out.ppm=out.ppm-(4.65-centreFreq);

%Fill in structure header fields:
out.seq='press';
out.te=tau1+tau2;
out.sim='shaped';

%Additional fields for compatibility with FID-A processing tools.
out.sz=size(out.specs);
out.date=date;
out.dims.t=1;
out.dims.coils=0;
out.dims.averages=0;
out.dims.subSpecs=0;
out.dims.extras=0;
out.averages=1;
out.rawAverages=1;
out.subspecs=1;
out.rawSubspecs=1;
out.flags.writtentostruct=1;
out.flags.gotparams=1;
out.flags.leftshifted=0;
out.flags.filtered=0;
out.flags.zeropadded=0;
out.flags.freqcorrected=0;
out.flags.phasecorrected=0;
out.flags.averaged=1;
out.flags.addedrcvrs=1;
out.flags.subtracted=1;
out.flags.writtentotext=0;
out.flags.downsampled=0;
out.flags.isFourSteps=0;

end
