%% ##############   design spiral grad        #############################
lims0 = mr.opts('MaxGrad',200,'GradUnit','mT/m',...
    'MaxSlew',300,'SlewUnit','T/m/s',...
    'adcDeadTime', 10e-6,'rfRasterTime',1e-6);

[spiral,MrProt,SeqLim] = sprial_prep(spiral,MrProt,SeqLim);

spiral.GradBuf(1:2,:) = spiral.GradBuf(1:2,:)*spiral.dMaxGradSpirAct*1000; %spiral.dMaxGradSpirAct unit is T/m)

spiral.adcDwell = (size(spiral.GradBuf,2)-spiral.RAMP_DOWN_POINTS)*spiral.GradDwellTime/spiral.paramLongROSamples;
spiral.adc = mr.makeAdc(spiral.paramLongROSamples ,'Dwell',spiral.adcDwell*1e-6);
spiral.grad_x = mr.makeArbitraryGrad('x',spiral.GradBuf(1,:)*larmor_freq,'system',lims0);
spiral.grad_y = mr.makeArbitraryGrad('y',spiral.GradBuf(2,:)*larmor_freq,'system',lims0);
if spiral.UseCAIPI
    spiral.GradBuf(3,:) = spiral.GradBuf(3,:)*spiral.caipiblip.GradAmp; %spiral.caipiblip.GradAmp unit is mT/m
    spiral.grad_z = mr.makeArbitraryGrad('z',spiral.GradBuf(3,:)*larmor_freq,'system',lims0);
end