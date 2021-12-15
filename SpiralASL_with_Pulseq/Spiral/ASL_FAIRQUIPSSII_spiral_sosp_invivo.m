buildpath = 'I:\pulseq-master\matlab\Spiral';
workpath = 'I:\pulseq-master\matlab';

addpath(fullfile(buildpath,'FPasl'));
addpath(buildpath);
addpath(fullfile(buildpath,'Parameter_define'));
addpath(fullfile(buildpath,'Spiral_grad'));
addpath(fullfile(buildpath,'RF_design'));
addpath(fullfile(buildpath,'Preparation_Blocks'));
addpath(fullfile(buildpath,'Pulseq_BuildingBlocks'));

cd(workpath)
%% sequence limit define

seq=mr.Sequence();          % Create a new sequence object

% Set system limits
lims = mr.opts('MaxGrad',200,'GradUnit','mT/m',...
    'MaxSlew',200,'SlewUnit','T/m/s',...
    'adcDeadTime', 10e-6,'rfRasterTime',1e-6);  %%'rfRingdownTime', 30e-6,'rfDeadtime', 100e-6

% (in Siemens interpreter from January 2019 duration is limited to 8.192 ms, and although product EPI uses 10.24 ms, 8 ms seems to be sufficient)
lims_sp = mr.opts('MaxGrad',200,'GradUnit','mT/m',...
    'MaxSlew',300,'SlewUnit','T/m/s',...
    'adcDeadTime', 10e-6,'rfRasterTime',1e-6);  %%'rfRingdownTime', 30e-6,'rfDeadtime', 100e-6
% spiral gradient seems for follow a differnt limit for slew rate, no idea why
B0=6.98; % 1.5 2.89 3.0
larmor_freq = 42.5756*1e3; %(MHz/T)
sat_ppm=-3.45;
sat_freq=sat_ppm*1e-6*B0*lims.gamma;

rfras_us = lims.rfRasterTime*1e6;
%%
%warning('Please check MR parameters');
%edit Exam_card_param_sosp;

%% ##############   load parameters           #############################
Exam_card_param_sosp_invivo;
MrProt.perfusion=1;
System_param;
SeqLim_param;
if MrProt.perfusion == 1
    PASL_param;
    PreSat_param;
    SupSat_param;
    InfSat_param
    HSIR_param;
end
CSatFat_param;
SpoilGrad_param;

spiral_param_sosp;


 RFExcite_sosp_param;
%% ##############    prepare the preparation pulses as well as the timing: FAIR QUIPSII + fatsat + spoilers
if MrProt.perfusion == 1
    prepare_trFOCI;
    prepare_PreSat;
    prepare_InfSat;
    prepare_SupSat;
end
prepare_CSatFat;
prepare_Spoil;

prepare_RF_timing;

%% ##############   prepare the delay between spoil gradient and Excitation unit #
delay_spoil_exc = mr.makeDelay(RFExcite.predelay * 1e-6);
%% ##############   prepare the excitation unit 
if strcmp(MrProt.private.rftype,'sinc')
    prepare_RFexcite;
elseif strcmp(MrProt.private.rftype,'rec')
    prepare_RFexcite_rec;
end

%% ##############   make delay before spiral readout but after skope trigger #####################
delay_pre_spiral = mr.makeDelay((RFExcite.postdelay/2-10)*1e-6);
%% ##############   prepare SOSP RF and gradients
prepare_grad_sosp;
prepare_RF_sosp;
prepare_spoil_sosp;

prepare_sosp_timing;

%% ##############   SKOPE: for skope sync scan
pre_sync_delay = mr.makeDelay(1e-6*(PASL.lSBBPaslDurationPerRequest_us+ ...
    RFExcite.predelay+(RFExcite.exc_GPSSel.duration+RFExcite.refoc.duration)*1e6));
inter_sync_delay = mr.makeDelay( (RFExcite.exc_GPSSel.duration+RFExcite.refoc.duration));
%% ##############   SKOPE: prepare skope trigger (makeDigitalOutputPulse  #############################
skope_trig = mr.makeDigitalOutputPulse('osc0', 'delay',RFExcite.postdelay/2*1e-6,'duration',10e-6,'system',lims);
 % possible channels: 'osc0','osc1','ext1'
skope_trig_delayreplace = mr.makeDelay(mr.calcDuration(skope_trig));
