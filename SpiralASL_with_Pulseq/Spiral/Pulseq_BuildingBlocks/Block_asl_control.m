%% ASL Control
for i = 1:PASL.m_lNPreSats
seq.addBlock(PreSat.pre_RF.rf,PreSat.pre_GPSSel.grad);
seq.addBlock(PreSat.pre_GPSpoil.grad);
end

seq.addBlock(SupSat.sus_RF.rf,SupSat.sus_GPSSel.grad);
seq.addBlock(SupSat.sus_GPSpoil.grad);

if SUPSAT_DELAY1 > 0
    seq.addBlock(mr.makeDelay(SUPSAT_DELAY1*1e-6));
end

seq.addBlock(rf_hsir_label,gz_hsir_control)
 seq.addBlock(HSIR.ss_GPSpoil.grad);

seq.addBlock(mr.makeDelay(PASL.m_lTI1FillTime_us*1e-6));

seq.addBlock(InfSat.ifs_RF.rf,InfSat.ifs_GPSSel.grad);
seq.addBlock(InfSat.ifs_GPSpoil.grad);

if INFSAT_DELAY1 > 0
seq.addBlock(mr.makeDelay(INFSAT_DELAY1*1e-6));
end

seq.addBlock(InfSat.ifs_RF.rf,InfSat.ifs_GPSSel.grad);
seq.addBlock(InfSat.ifs_GPSpoil.grad);

if INFSAT_DELAY2 > 0
seq.addBlock(mr.makeDelay(INFSAT_DELAY2*1e-6));
end
      
seq.addBlock(InfSat.ifs_RF.rf,InfSat.ifs_GPSSel.grad);
seq.addBlock(InfSat.ifs_GPSpoil.grad);

if INFSAT_DELAY3 > 0
seq.addBlock(mr.makeDelay(INFSAT_DELAY3*1e-6));
end

seq.addBlock(InfSat.ifs_RF.rf,InfSat.ifs_GPSSel.grad);
seq.addBlock(InfSat.ifs_GPSpoil.grad);

if INFSAT_DELAY4 > 0
seq.addBlock(mr.makeDelay(INFSAT_DELAY4*1e-6));
end
   
if PASL.m_lTI2FillTime_us > 0
seq.addBlock(mr.makeDelay(PASL.m_lTI2FillTime_us*1e-6));
end