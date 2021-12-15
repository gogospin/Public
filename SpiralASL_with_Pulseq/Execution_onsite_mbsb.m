%clear all
MrProt.sliceGroupList(1).shift = 0; % 5.7 L, 28.6 A, 48.0 H
MrProt.sliceGroupList(1).thickness = 1.2; %mm
MrProt.sliceGroupList(1).sliceperslab = 40;
MrProt.dynamics = 10;%;
MrProt.flipAngle = 60;
MrProt.sliceSeries.size = MrProt.sliceGroupList(1).sliceperslab;
MrProt.sliceGroupList(1).distanceFactor = 50;
MrProt.spcl.g3.Spiral_Type = 0;
MrProt.spcl.single_shot_mode = 1;
MrProt.spcl.accr_kz = 1;
MrProt.spcl.private.rotate_mode_ax = 'none';% or 'golden_angle'
MrProt.spcl.private.increm_mode_sl = 'linear';% or 'binomial' or 'arbit'
MrProt.spcl.g3.RO_samples = 10240*1.2;

MrProt.sliceGroupList(1).shiftflag=1;
run('I:\pulseq-master\matlab\Spiral\ASL_FAIRQUIPSSII_spiral_mbsb_invivo.m')
InfSat.ifs_RF.rf.freqOffset = -InfSat.ifs_RF.rf.freqOffset;
PreSat.pre_RF.rf.freqOffset = -PreSat.pre_RF.rf.freqOffset;
SupSat.sus_RF.rf.freqOffset = -SupSat.sus_RF.rf.freqOffset;
CSatFat.csat_RF.rf.freqOffset = -CSatFat.csat_RF.rf.freqOffset;
run('I:\pulseq-master\matlab\Spiral\ASL_FAIRQUIPSSII_spiral_mbsb_addBlock_invivo.m')
%seq.write('C:\Users\P70071217\WMWare_share\pulseq\sequences\test_0411\asl_mbsb.seq');
%% save data
save_dir = 'Z:\test\to_test_07_30\'; % change it accordingly
RZ = MrProt.spcl.private.increm_mode_sl;
RXY = MrProt.spcl.private.rotate_mode_ax;
iPAT = num2str(MrProt.spcl.accr_kz);
if strcmp(RXY, 'golden_angle')
    RXY =  'goldang';
elseif strcmp(RXY, 'none')
    RXY = 'norot';
end

namestr = strcat(save_dir,'external_h45_','iPAT',iPAT,'_RZ',RZ,'_RXY',RXY,'.seq');
seq.write(namestr);


%% optional: if you want to see the trajectory
counti=1;
for i = 1:10240*1.2*2.5/10-2
    if counti == 1
        kx_evol(counti) = spiral.GradBuf(1,i);
        ky_evol(counti) = spiral.GradBuf(2,i);
        kz_evol(counti) = spiral.caipiblip.Grad(i)*spiral.caipiblip.GradAmp;
    else
        kx_evol(counti) = kx_evol(counti-1)+spiral.GradBuf(1,i);
        ky_evol(counti) = ky_evol(counti-1)+spiral.GradBuf(2,i);
        kz_evol(counti) = kz_evol(counti-1)+spiral.caipiblip.Grad(i)*spiral.caipiblip.GradAmp;
    end
    counti = counti+1;
end
kx_evol = single(kx_evol)/100*2*pi*43.576;
ky_evol = single(ky_evol)/100*2*pi*43.576;
kz_evol = single(kz_evol)/100*2*pi*43.576;
