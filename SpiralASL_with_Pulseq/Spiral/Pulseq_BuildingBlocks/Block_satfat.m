%% SatFat and spoiler blocks

seq.addBlock(CSatFat.csat_grad.grad_x1,CSatFat.csat_grad.grad_y1,CSatFat.csat_grad.grad_z1);
seq.addBlock(CSatFat.csat_RF.rf);
seq.addBlock(CSatFat.csat_grad.grad_x2,CSatFat.csat_grad.grad_y2,CSatFat.csat_grad.grad_z2);

seq.addBlock(SpoilGrad.grad_y,SpoilGrad.grad_z);

seq.addBlock(delay_spoil_exc);