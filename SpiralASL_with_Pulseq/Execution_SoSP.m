clear all
MrProt.sliceGroupList(1).shift = 0; % 5.7 L, 28.6 A, 48.0 H -> offcenter shift
MrProt.sliceGroupList(1).thickness = 2; %mm 
MrProt.sliceGroupList(1).sliceperslab = 24;
MrProt.dynamics = 2;
MrProt.flipAngle = 10;

MrProt.spcl.g3.Spiral_Type = 0;
MrProt.spcl.single_shot_mode = 1;
MrProt.spcl.accr_kz = 1;
MrProt.spcl.private.rotate_mode_ax = 'none';%'none';%'caipi';% or 'golden_angle'
MrProt.spcl.private.increm_mode_sl = 'linear';% or 'binomial' or 'arbit'
MrProt.spcl.g3.RO_samples = 10240;

MrProt.sliceGroupList(1).shiftflag=1;
run('I:\pulseq-master\matlab\Spiral\ASL_FAIRQUIPSSII_spiral_sosp_invivo.m')
InfSat.ifs_RF.rf.freqOffset = -InfSat.ifs_RF.rf.freqOffset;
PreSat.pre_RF.rf.freqOffset = -PreSat.pre_RF.rf.freqOffset;
SupSat.sus_RF.rf.freqOffset = -SupSat.sus_RF.rf.freqOffset;
CSatFat.csat_RF.rf.freqOffset = -CSatFat.csat_RF.rf.freqOffset;
run('I:\pulseq-master\matlab\Spiral\ASL_FAIRQUIPSSII_spiral_sosp_addBlock_invivo.m')
seq.write('test.seq');
