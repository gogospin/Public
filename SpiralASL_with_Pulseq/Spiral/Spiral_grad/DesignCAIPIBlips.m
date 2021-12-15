function spiral = DesignCAIPIBlips( paramLongSpGradAmp, paramLongSpSlewRate, ...
    Separation, Delay_us, Period_us, TotalTime_us,spiral )

dGamma  = 42576000.0; 
GRAD_RASTER_TIME = 10;

MaxSlewRate=paramLongSpSlewRate/(pi()*1000.0);

DesiredMoment=1000000*0.5/(Separation*dGamma/1000000.0);
	
MaxBlockMoment=paramLongSpSlewRate*GRAD_RASTER_TIME*GRAD_RASTER_TIME/1000.0;

dBlockNeeded=DesiredMoment/MaxBlockMoment;
   
ClosestOdd=ceil(sqrt(dBlockNeeded)*2+1);
    
Closest=ClosestOdd;
    
ClosestEvenHalf=(ClosestOdd-1)/2;
if (ClosestEvenHalf+1)*ClosestEvenHalf>dBlockNeeded
    Closest=ClosestOdd-1;
end

nSteps=Closest;
    
FHalfnSteps=floor(nSteps/2);
    
if(mod(nSteps,2)==0)
    nBlocks=FHalfnSteps*(FHalfnSteps+1);
else
    nBlocks=FHalfnSteps*(FHalfnSteps+1)+FHalfnSteps+1;
end
    
BlockMoment=DesiredMoment/nBlocks;
BlockH=BlockMoment/GRAD_RASTER_TIME;
    
nTotalSteps=TotalTime_us/GRAD_RASTER_TIME;
VecLen=nSteps;
    
GradAmp=BlockH*FHalfnSteps;
if(mod(nSteps,2)==1)
    GradAmp=BlockH*(FHalfnSteps+1);
end
GradAmp=max(GradAmp,0.0001);

for i=1:FHalfnSteps
    TmpGrad(i)=BlockH*(i)/GradAmp;
    TmpGrad(nSteps-i+1)=BlockH*(i)/GradAmp;
end
if mod(nSteps,2)==1
    TmpGrad(FHalfnSteps+2)=BlockH*(FHalfnSteps+1)/GradAmp;
end

DelaySteps=Delay_us/GRAD_RASTER_TIME;
Grad = zeros(nTotalSteps,1);
        
PeriodSteps=Period_us/GRAD_RASTER_TIME;
    if(PeriodSteps<VecLen)  
        PeriodSteps=VecLen; 
    end
    
    j=DelaySteps;
    DurDir=1;
    
    while(j<nTotalSteps) 
        if(j<nTotalSteps-VecLen) 
            for i=1:VecLen 
                Grad(i+j)=DurDir*TmpGrad(i);
            end
        end
        j = j+PeriodSteps;
        DurDir = DurDir*(-1);
    end
    
spiral.caipiblip.Grad = Grad;
spiral.caipiblip.GradAmp = GradAmp;

end

