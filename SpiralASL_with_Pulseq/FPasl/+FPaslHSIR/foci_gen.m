function [sFoci, psRF1, psRF2,  psGrad_ARB] = foci_gen(MrProt, psRF1, psRF2,  psGrad_ARB, sFoci)

% This code was adopted from the original IDEA pulse code from Renzo Huber 
% see: https://layerfmri.com/2018/11/29/tr-foci-pulse-optimisations/#more-1225
%--------------------------------------------------------------
% macro defs for calculation of the phase of the RF pulse
% macro's are used instead of functions to increase calculation speed, though this might not be relevant at all
% these are analytically egrated expressions of frequency in paper by Payne and Leach
%--------------------------------------------------------------

%%
GRAD_RASTER_TIME = 10;
%%
PHI_HS = @(Mu, BetaT) (-(Mu)*log(cosh(BetaT)));
PHI_FOCI = @(Mu, BetaT) (-(Mu)*(cosh(BetaT)-1));

%% with frequency offset DeltaOmega [rad/s] for a slice shift
PHI_HS_FO = @(Mu,BetaT,DeltaOmega,T) (-(Mu) * log(cosh(BetaT)) + (DeltaOmega)*(T));
PHI_FOCI_FO = @(Mu,BetaT,DeltaOmegaDivBeta) (-(Mu) * (cosh(BetaT)-1) + (DeltaOmegaDivBeta) * sinh(BetaT));

%% parameters

dummy = 0.0;
dummy1 = 0.0;
debug = 0;

b1max =0. ;

Gs = 1.0;

Amax = 3.32;
w = 0.30;
r1 = 0.64;
r2 = 0.27;
r3 = 0.59;
r4 = 0.00;
r5 = 1.00;

beta = 3.90;
t1 = 0.25;
t2 = 0.40;
B1max =  0.327152;

c_ = Amax*(1.-r1);
Amin = r2*c_;
b1 = r3*(c_-Amin)/((w-1)*(w-1));
b2 = r4*(1-r3)*(c_-Amin)/((w-1)*(w-1)*(w-1)*(w-1));
b3 = r5*(1-r4)*(1-r3)*(c_-Amin)/((w-1)*(w-1)*(w-1)*(w-1)*(w-1)*(w-1));
b4 = (1-r5)*(1-r4)*(1-r3)*(c_-Amin)/((w-1)*(w-1)*(w-1)*(w-1)*(w-1)*(w-1)*(w-1)*(w-1));

A1 = @( t) ( Amax*(1.-r1*(t+1.)/w ));
A2 = @( x) ( Amin + b1*x*x + b2*x*x*x*x + b3*x*x*x*x*x*x + b4*x*x*x*x*x*x*x*x);
A3 = @( t) ( Amax*(1.-r1*(-t+1.)/w ));
T = @( x)  ( (t1*x*x*x*x*x+t2*x*x*x+x)/(t1+t2+1));
grad1 = @( x) ( A1(x)*Gs/A1(-1));
grad2 = @( x) ( A2(x)*Gs/A1(-1));
grad3 = @( x) ( A3(x)*Gs/A1(-1));
sech = @( x)  (2/(exp(x)+exp(-x)));
mytanh = @( x) ( (exp(x)-exp(-x))/(exp(x)+exp(-x)));

B1_1 = @(x)  (A1(x)*sech(beta* T(x))/B1max );
B1_2 = @(x)  (A2(x)*sech(beta* T(x))/B1max );
B1_3 = @(x)  (A3(x)*sech(beta* T(x))/B1max );
freq1 = @(x) (-A1(x)*beta*mytanh(beta*T(x))); 
freq2 = @(x) (-A2(x)*beta*mytanh(beta*T(x))); 
freq3 = @(x) (-A3(x)*beta*mytanh(beta*T(x))); 




%% local variables for defining the pulse
m_dB1Cutoff = 50; %RENZO: Es scheint ziemlich egal zu sein, war da steht, da m_dB1Cutoff eh in ...Bassi_RFPulses.cpp definiert wird
if isfield(sFoci,'m_dFOCI_Cutoff')
    if m_dB1Cutoff ~= sFoci.m_dB1Cutoff
        m_dB1Cutoff = sFoci.m_dB1Cutoff;
    end
end
m_dFOCI_Cutoff = 20.0;
if isfield(sFoci,'m_dFOCI_Cutoff')
    if m_dFOCI_Cutoff ~= sFoci.m_dFOCI_Cutoff
        m_dFOCI_Cutoff = sFoci.m_dFOCI_Cutoff;
    end
end
m_dMu = 7.71;
if isfield(sFoci,'m_dMu')
    if m_dMu ~= sFoci.m_dMu
        m_dMu = sFoci.m_dMu;
    end
end
m_lGradRampTime = 500;    % sets gradient ramp up/down time
if isfield(sFoci,'m_lGradRampTime')
    if m_lGradRampTime ~= sFoci.m_lGradRampTime
        m_lGradRampTime = sFoci.m_lGradRampTime;
    end
end
m_dBeta = 0.0;
m_dBW_kHz = 0.0;
m_dTs = 0.0;
m_bIsRTUpdate = false;    %% RTFOCI
m_dShift1 = -19522.0;
m_dShift2 = -19522.0;
m_bCalculated = false;
m_lRF1Samples = 0;
m_lRF2Samples = 0;
m_lGradSamples = 0;
m_lGradRampSamples = 0;
m_lGradFltSamples = 0;
m_lGradTotalTime = 0;
m_dGradMomentumTOT = 0.0;
m_dAmplIntegral = 0.0;
m_psRF1Sample = [];
m_psRF2Sample = [];
m_pfGradSample = [];

% read from MrProt
l_LongEntry_trfociAmpl = MrProt.spcl.g1.trfociAmpl;
l_LongEntry_trfociBWDTH = MrProt.spcl.g1.trfociBWDTH;
l_LongEntry_trfociD = MrProt.spcl.g1.trfociD;



%% prepare
dTwoPi = 2.0 * pi();

iN = 0;  %% counter

wip4  = 527. * l_LongEntry_trfociAmpl/ 100 ; %% AMPL
wip5  = l_LongEntry_trfociBWDTH/ 100 ; %% BWDTH %% in procent
wip6  = 0; % will be filled below
wip7  = 10. ; %l_LongEntry_trfociD ; % Thickness in procent
wip8  = l_LongEntry_trfociD ;
wip11 = 5120. ;
wip12 = 100. ;
% dimo - make the same as in SMS ASL

% Geom for offset calculation

e_position = MrProt.sliceGroupList(1).position;
wip6 = MrProt.sliceGroupList(1).shift ;
%wip6 = 0.; %/ Renzo doesn't wants to have it global only  pMrProt->sliceGroupList()[0].shift() ;
abs_offset = wip6 ;

m_dTs = 1.0E-6 * psRF1.Duration / psRF1.SampleSize;


%--------------------------------------------------------------
% HS pulse parameters
%--------------------------------------------------------------
% Steffen parameter Map

% special card -> bandwidth
mu = 7.71* wip5;


%RENZO:auskommentiert: m_dBeta = (2.0 / (1.0E-6 * psRF1->getDuration())) * acosh (1.0 / m_dB1Cutoff);

% bandwidth in kHz of the HS pulse %ist auf RENZOS Änderungen angepasst %brauch ich, um Gradient anzupassen
m_dBW_kHz = 1.0E-3 * mu * beta* Amax / pi()* 1/(1.0E-6 *psRF1.Duration);

fprintf('%s %d\n', '=== Tr-FOCI bandwidth in kHz (m_dBW_kHz): ', m_dBW_kHz);

%--------------------------------------------------------------
% FOCI pulse parameters
%--------------------------------------------------------------
% time of switch between constant and modulated frequency RENZO: in C-FOCI, in tr-FOCI the modulation funtions change

%RENZO:auskommentiert :double dT_Cutoff = acosh (m_dFOCI_Cutoff) / m_dBeta;

%--------------------------------------------------------------
% Gradient amplitude for FOCI pulse
%--------------------------------------------------------------
%Renzo hat ausgeklammer

% Steffen: direkt GS amplitude auf slab thickness einstellen -> Test mit alterm Ausdruck!
psRF1.GSAmplitude = psRF1.RequiredGSPolarity * 2* pi()* m_dBW_kHz  / (267522209. * wip8) * 1.0E9; % Deni's version

% Amplitude in mT/m 			%							       in HZ  	%267522209 ist Gamma  % in Meter

psRF2.GSAmplitude = psRF1.GSAmplitude;

% Gradient requirements: check with hardware performance


%--------------------------------------------------------------
% calculate FOCI pulse samples
%--------------------------------------------------------------
% pre-calculate variables used in loop
dN0 = (psRF1.SampleSize - 1) / 2.0;

%----------------------------------------------------------------------
% Calculate RF samples
%----------------------------------------------------------------------

%Offset Anfang

gradientenamplitude = psRF1.GSAmplitude;%*(pMrProt->wipMemBlock().adFree[2]-1.) /100.; %mT/m  soll mal duch parameter von special card ersetzt werden

%renzo_Shift = 0.;  %in mm   %RENZO wieder auskommentieren

dummy = 0.0 ;
dummy1 = 0.0 ;
renzo_dummyphase = 401.* 2.* 3.141596;
renzoN = 5120  ; % sollte die Anzahl der Samples sein psRF1->getSamples() Dass muss geändert werden, wenn man was an psRF1->getSamples() aendert, z.B. in RFexite.cpp
phase = zeros(renzoN,1) ;
phase1 = zeros(renzoN,1) ;
%berechne Phase aus numerischem Integral über Frequenz

renzo_Shift = -0.2* psRF1.LarmorConst * psRF1.GSAmplitude * abs_offset / dTwoPi * 2;

for  j = 0: renzoN-1
    
    if(double(j) < double(renzoN) * w/2)
        dummy  = dummy - (renzo_Shift * grad1(double(-1.+2.* double(j)/double(renzoN)))+mu * freq1(double(-1. +  2.*double(j)/double(renzoN))))/double(renzoN);
        dummy1 = dummy1- ( 0.* grad1(double(-1.+2.*double(j)/double(renzoN))) + mu *freq1(double(-1. +  2.*double(j)/double(renzoN))))/double(renzoN);
        
    elseif ((double(j) >= double(renzoN* w/2)) && (double(j) <= double(renzoN) * (1.-w/2.)))
        dummy  = dummy  - (renzo_Shift* grad2(double(-1.+2.*double(j)/double(renzoN))) + mu *freq2(double(-1. + 2.*double(j)/double(renzoN))))/double(renzoN);
        dummy1 = dummy1 - ( 0.* grad2(double(-1.+2.*double(j)/double(renzoN))) + mu *freq2(double(-1. + 2.*double(j)/double(renzoN))))/double(renzoN);
        if ( j == renzoN /2 )
            rad = 0.; % ((double)pMrProt->wipMemBlock().adFree[0]-1.)/180. * 3.1415927;
         
            %dummy = dummy + wip8 /180. * 3.1415927 ;
            %dummy1 = dummy1 + wip8 /180. * 3.1415927 ;
        end
    elseif(double(j) > double(renzoN)* (1.-w/2.))
        dummy  = dummy  +(-renzo_Shift * grad1(-double(-1.+2.*double(j)/double(renzoN)))+mu *freq1(double(1. -  2.*double(j)/double(renzoN) )))/double(renzoN);
        dummy1 = dummy1 +(-0. * grad1(-double(-1.+2.*double(j)/double(renzoN)))+mu *freq1(double(1. -  2.*double(j)/double(renzoN) )))/double(renzoN);
    end
    phase(j+1)  = dummy  - renzo_dummyphase; % dummy Phase ist eingefühgt, weil Phasenverlauf nicht positiv werden darf. (sonst Bug). Weil aber
    phase1(j+1) = dummy1 - renzo_dummyphase;
    % später eh alles auf Rest zu 2 Pi gekürtzt wird, dürfte das kein Problem darstellen
end

dAmplSumRe=0.0; dAmplSumIm=0.0; renzoInt =0.0; 
for iN = 0: psRF1.SampleSize-1
    
    dTn = 0;
    dPha = 0.;
    dummy = 0.0 ;
   
    % Steffen Parameter Map
    %dTn = (iN - dN0) * m_dTs * 2./(5.12E-3)* 5120./pMrProt->wipMemBlock().adFree[12];  % Zeit im Pulse normiert auf -1 bis 1
    dTn = (iN - dN0) * m_dTs * 2./(5.12E-3)* 5120./ 2./wip11;  % Zeit im Pulse normiert auf -1 bis 1
    %ddn(iN+1) = dTn;
    if (dTn < w-1.)
        
        m_psRF1Sample(iN+1).flAbs  = single(B1_1(single(dTn)));
        m_psRF2Sample(iN+1).flAbs  = single(B1_1(single(dTn)));
        
        m_psRF1Sample(iN+1).flPha = -single(fmod(phase(iN+1),2.*pi()));
        m_psRF2Sample(iN+1).flPha = -single(fmod(phase1(iN+1),2.*pi()));
        
        dPha = -single(fmod(phase(iN+1),2.*pi()));
    elseif (dTn <= 1.-w && dTn >= w-1. )
        
        m_psRF1Sample(iN+1).flAbs  = single(B1_2(single(dTn)));
        m_psRF2Sample(iN+1).flAbs  = single(B1_2(single(dTn)));
        
        
        m_psRF1Sample(iN+1).flPha = single(fmod(dummy,pi()));
        m_psRF1Sample(iN+1).flPha = -single(fmod( phase(iN+1),2.*pi()));
        m_psRF2Sample(iN+1).flPha = -single(fmod(phase1(iN+1),2.*pi()));
        dPha = -single(fmod(phase(iN+1),2.*pi()));
        % the above commands in cpp are fmod, which can result in negative
        % number, following the strict number of mod
        % but in matlab mod always returns a positive number
      
    elseif (single(dTn) > (1-w))
        
        m_psRF1Sample(iN+1).flAbs  = single(B1_3(single(dTn)));
        m_psRF2Sample(iN+1).flAbs  = single(B1_3(single(dTn)));
        
         m_psRF1Sample(iN+1).flPha = -single(fmod( phase(iN+1),2.*pi()));
         m_psRF2Sample(iN+1).flPha = -single(fmod(phase1(iN+1),2.*pi()));
         dPha = -single(fmod(phase(iN+1),2.*pi()));
         %disp(iN+1)
         %m_psRF1Sample(iN+1).flPha = phase(iN+1);
        % m_psRF2Sample(iN+1).flPha = phase1(iN+1);
        
    end
    % added for correction btw cpp and matlab
    m_psRF1Sample(iN+1).flPha = mod(m_psRF1Sample(iN+1).flPha,2*pi);
    m_psRF2Sample(iN+1).flPha = mod(m_psRF2Sample(iN+1).flPha,2*pi);
    % end of correction
    renzoInt = renzoInt + m_psRF1Sample(iN+1).flAbs;
    % for pulseq:
    %dPha = - dPha;
    % sum for amplitude integral
    dAmplSumRe = dAmplSumRe + m_psRF1Sample(iN+1).flAbs * cos(dPha);
    dAmplSumIm = dAmplSumIm + m_psRF1Sample(iN+1).flAbs * sin(dPha);
    
end

% plot the RF:
for iN = 1:length(m_psRF1Sample)
    rf1_sample(iN) = m_psRF1Sample(iN).flAbs*exp(1i*m_psRF1Sample(iN).flPha);
end

%figure,plot(abs(rf1_sample))
%----------------------------------------------------------------------
% RF amplitude integral
%----------------------------------------------------------------------
m_dAmplIntegral = sqrt (dAmplSumRe*dAmplSumRe + dAmplSumIm*dAmplSumIm) ;% * m_dTs/0.83295324;   %renzo
sFoci.FA_scale =  0.6427877/0.9874/m_dAmplIntegral*1067;

% Steffen test
% m_dAmplIntegral = m_dAmplIntegral * 90.721917/wip4;
% m_dAmplIntegral = m_dAmplIntegral * 5120./ wip11 ;


%--------------------------------------------------------------
% calculate gradient data
%--------------------------------------------------------------
% -dsb 11/4/05: psGrad changed to psGrad_ARB

m_lGradRampSamples = m_lGradRampTime / GRAD_RASTER_TIME;
m_lGradFltSamples = psRF1.Duration / GRAD_RASTER_TIME;
m_lGradSamples = 2 * m_lGradRampSamples + m_lGradFltSamples;
m_lGradTotalTime = m_lGradSamples * GRAD_RASTER_TIME;

% allocate memory for samples if required
m_pfGradSample = zeros(m_lGradSamples,1);

%----------------------------------------------------------------------
% FOCI-modulated "flat" top data
%----------------------------------------------------------------------
dN0 = (m_lGradFltSamples - 1) / 2.0;
m_dGradMomentumTOT = 0.0;

for iN = 0:m_lGradFltSamples-1
    dGTs = 1.0E-6 * GRAD_RASTER_TIME;  % gradient raster time [s]  (Zeit pro sample)
    dTn = single(2.* (iN - dN0) * dGTs/(5.12E-3)* 5120./wip11 / 2.); %  Zeitnormierung anpassen bei veränderlicher Pulsdauer nach special card
    % die 2 muss da hin weil ich zu viele Samples hab
    %RENZO GRAD ANFANG  %Gradienten verlauf geht von 0 bin 1
    if (dTn < w-1)
        m_pfGradSample(m_lGradRampSamples+iN+1) = single(grad1(dTn));
    elseif (dTn > (1-w))
        m_pfGradSample(m_lGradRampSamples+iN+1) = single(grad1(-dTn));
    end
    if (abs(dTn) < 1-w )
        m_pfGradSample(m_lGradRampSamples+iN+1) = single(grad2(dTn));
    end
    %RENZO GRAD ENDE
    m_dGradMomentumTOT = m_dGradMomentumTOT + m_pfGradSample(m_lGradRampSamples+iN+1);
end

%----------------------------------------------------------------------
% ramp up/down data
%----------------------------------------------------------------------
for iN = 0 : m_lGradRampSamples-1
    % ramp data
    fRampData = single(iN/m_lGradRampSamples);
    % ramp up
    m_pfGradSample(iN+1) = fRampData;
    % ramp down
    m_pfGradSample(m_lGradSamples-iN) = fRampData;
    % moment
    m_dGradMomentumTOT = m_dGradMomentumTOT + 2.0 * fRampData;
end
%----------------------------------------------------------------------
%  download shape to gradient controller
%----------------------------------------------------------------------

%psGrad_ARB.RampShape((float*) m_pfGradSample, (long) m_lGradSamples, 0, true);
psGrad_ARB.RampShape = m_pfGradSample;
psGrad_ARB.Amplitude = psRF1.GSAmplitude;

psGrad_ARB.Duration = m_lGradTotalTime;

%----------------------------------------------------------------------
% scale gradient moment
%----------------------------------------------------------------------
m_dGradMomentumTOT = m_dGradMomentumTOT* psGrad_ARB.Amplitude * GRAD_RASTER_TIME;


%--------------------------------------------------------------
% done!
%--------------------------------------------------------------
m_bCalculated = true; % don't calculate the same twice!

sFoci.m_bCalculated = m_bCalculated;
sFoci.m_dB1Cutoff = m_dB1Cutoff ;
sFoci.m_dFOCI_Cutoff = m_dFOCI_Cutoff ;
sFoci.m_dMu = m_dMu ;
sFoci.m_lGradRampTime = m_lGradRampTime ;
sFoci.m_dBeta = m_dBeta ;
sFoci.m_dBW_kHz = m_dBW_kHz;
sFoci.m_dTs = m_dTs ;
sFoci.m_bIsRTUpdate = m_bIsRTUpdate;
sFoci.m_bCalculated = m_bCalculated;
sFoci.m_lRF1Samples = m_lRF1Samples ;
sFoci.m_lRF2Samples = m_lRF2Samples ;
sFoci.m_lGradSamples = m_lGradSamples;
sFoci.m_lGradRampSamples = m_lGradRampSamples ;
sFoci.m_lGradFltSamples = m_lGradFltSamples ;
sFoci.m_lGradTotalTime = m_lGradTotalTime ;
sFoci.m_dGradMomentumTOT = m_dGradMomentumTOT ;
sFoci.m_dAmplIntegral = m_dAmplIntegral ;
sFoci.m_psRF1Sample = m_psRF1Sample;
sFoci.m_psRF2Sample = m_psRF2Sample;
sFoci.m_pfGradSample = m_pfGradSample ;
sFoci.B1max = B1max;
sFoci.phase = phase;
sFoci.phase1 = phase1;



%%
for iN = 1:length(m_psRF1Sample)
    rf1_sample(iN) = m_psRF1Sample(iN).flAbs*exp(1i*m_psRF1Sample(iN).flPha);
    rf2_sample(iN) = m_psRF2Sample(iN).flAbs*exp(1i*m_psRF2Sample(iN).flPha);
end
% figure,
% subplot(1,3,1),plot(real(rf1_sample));
% subplot(1,3,2),plot(angle(rf1_sample))
% subplot(1,3,3),plot(m_pfGradSample);

psRF1.Samples = rf1_sample;
psRF2.Samples = rf2_sample;
psRF1.Asymmetry = 0.5;
psRF2.Asymmetry = 0.5;
psRF1.aslstate = 'label';
psRF2.aslstate = 'control';

psGrad_ARB.Samples = m_pfGradSample * psRF1.GSAmplitude;
