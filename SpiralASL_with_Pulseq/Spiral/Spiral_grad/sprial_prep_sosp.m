function [spiral,MrProt,SeqLim] = sprial_prep_sosp(spiral,MrProt,SeqLim)
slInnerSpSlewRate=spiral.paramLongSpSlewRate;
slInnerSpGradAmp=spiral.paramLongSpGradAmp;
lSetsToMeasure=spiral.paramLongShiftsToMeasure;
lSlicesToMeasure = MrProt.sliceSeries.size;
lSpiralType = spiral.paramLongSpiralType;
paramLongROSamples = spiral.paramLongROSamples;
paramLongSpGradAmp = spiral.paramLongSpGradAmp;
paramLongSpSlewRate = spiral.paramLongSpSlewRate;
paramLongInterleaves = spiral.paramLongInterleaves;

SeparationInmm = spiral.paramDoubleEchoShift;
CAIPIWantedInterval = spiral.paramDoubleCAIPIPeriod;
CAIPIDelay = spiral.paramDoubleCAIPIDelay;

RAMP_DOWN_POINTS = 50;   % Ramp down of the gradients at the end of a spiral gradient
spiral.RAMP_DOWN_POINTS = RAMP_DOWN_POINTS;

VD = spiral.VD;
spBW = spiral.spBW;
AccR = spiral.AccR;
nInterleaves = spiral.paramLongInterleaves;

SeqLim.ReadoutOSFactor = 1.0;
lADCsToMeasure = spiral.paramLongROSamples/1024;
spiral.lADCsToMeasure = lADCsToMeasure;
spiral.adcperSegment = 1024;
MrProt.fastImaging.segments = lADCsToMeasure;		
%dSamplingIntBase = MrProt.rxSpec.effDwellTime/SeqLim.ReadoutOSFactor;
dSamplingIntBase=1e9/spiral.spBW;

BASE_RES = spiral.BASE_RES;

%Dwell = MrProt.rxSpec.effDwellTime/SeqLim.ReadoutOSFactor;

f=SeqLim.ReadoutOSFactor;
LL=floor(dSamplingIntBase + 0.5);
dFOV = MrProt.sliceSeries.aFront.readoutFOV / 1000.0;

PreparedGrad = (lSpiralType > 9);
if lSpiralType > 9
    error('the code has not been implemented yet, please go to line 277 of gSpiral.cpp');
else
    if VD < 0.1
        error('the code has not been implemented yet, please go to line 301 of gSpiral.cpp');
    else
        nPointsADC=paramLongROSamples;
        spiral.OutInOut = (lSpiralType==-3);
        spiral.OutIn = (lSpiralType==-2);
        
        if(lSpiralType>1)
            spiral.OutInx = lSpiralType;
            nPointsADC = paramLongROSamples;
        end
        if(spiral.OutInOut)
            nPointsADC = nPointsADC/3;
        end
        if(spiral.OutIn)
            nPointsADC=nPointsADC/2;
            MrProt.contrasts = 2;
        end
        if(lSpiralType>1)
            nPointsADC = nPointsADC/lSpiralType;
        end
        BW = spBW;
        FOVm = dFOV;
        FOV = FOVm;
        Acc = AccR;
        nInterleaves = spiral.paramLongInterleaves;
        alpha=VD;
        
        NTrajAll = nPointsADC*Acc*nInterleaves;
        TendWanted = nPointsADC/BW;
        BaseRes=sqrt( NTrajAll*4/pi());
        
        trueRes=BaseRes;
        
        if (mod(trueRes, BASE_RES) >= (BASE_RES/2.0))
            xyRes = floor(ceil(trueRes/BASE_RES)) * floor(BASE_RES);
        else
            xyRes = floor(floor(trueRes/BASE_RES)) * floor(BASE_RES);
        end
        
        MaxGradAmp=paramLongSpGradAmp/(2*pi()*1000.0);
        MaxSlewRate=paramLongSpSlewRate/(2*pi());
        
        res=FOVm/BaseRes;
        N1=FOV/res;
        
        NLow=0.01;
        NHigh=N1/(2*nInterleaves)-0.001;
        NCur=1;
        
        TEndOut=0;
        
        nIter=0;
        while(abs(TEndOut-TendWanted)>1e-6)
            nIter = nIter+1 ;
            if(nIter > 20)
                break;
            end
            TEndOut = VDgetTend(FOV, res, NCur, nInterleaves, alpha, MaxGradAmp, MaxSlewRate, spiral);
            if(TEndOut<TendWanted)
                NHigh=NCur;
            else
                NLow=NCur;
            end
            NCur=(NHigh+NLow)/2;
            fprintf("nIter = %d,NCur = %.7f, TEndOut = %.5f\n",nIter,NCur,TEndOut);
        end
        
        ExtraDelay=spiral.paramLongZCToFlat;
        spiral = VDSpiral_sosp_V1(dFOV, dFOV/BaseRes, NCur, nInterleaves, alpha, MaxGradAmp,MaxSlewRate,1,ExtraDelay,spiral.nsosp_segments,spiral.nz,spiral.FOVz,spiral.nsosp_rz,spiral); 
        % test
%         for ll = 1:4
%             for i = 1:3
%                 counti = 1;
%                 for j = 1:4
%                     for k = 1:length(spiral.gx_seg(ll).shot(i).kslice(j).grad)
%                         if counti == 1
%                             kx_evol(ll,i,counti) = spiral.gx_seg(ll).shot(i).kslice(j).grad(k);
%                             ky_evol(ll,i,counti) = spiral.gy_seg(ll).shot(i).kslice(j).grad(k);
%                             kz_evol(ll,i,counti) = spiral.gz_seg(ll).shot(i).kslice(j).grad(k);
%                         else
%                             kx_evol(ll,i,counti) = kx_evol(ll,i,counti-1)+spiral.gx_seg(ll).shot(i).kslice(j).grad(k);
%                             ky_evol(ll,i,counti) = ky_evol(ll,i,counti-1)+spiral.gy_seg(ll).shot(i).kslice(j).grad(k);
%                             kz_evol(ll,i,counti) = kz_evol(ll,i,counti-1)+spiral.gz_seg(ll).shot(i).kslice(j).grad(k);
%                         end
%                         counti = counti+1;
%                     end
%                 end
%             end
%         end
        
        
        % STILL NEED TO IMPLEMENT RAMPDOWN
        
        
        
        
        
%         ulSpirSamples = spiral.ulSpirSamples;
%         
%         MaxLoc=ulSpirSamples-1;
%         if(spiral.OutIn)
%             MaxLoc=ulSpirSamples/2-1;
%         end
%         
%         if(lSpiralType>1)
%             MaxLoc=ulSpirSamples/lSpiralType-1;
%         end
%         
%         EffMaxRes=0;
%         CurR=0;
%         for ii=1:MaxLoc
%             CurR=norm([spiral.TrajBuf(1,ii),spiral.TrajBuf(2,ii)]);
%             EffMaxRes=max(EffMaxRes,CurR);
%         end
%         EffMaxRes = EffMaxRes*2*FOVm/(2*pi());
%         
%         EffRes_mm=1000.0*FOVm/EffMaxRes;
        
        % CAIPI prep goes here, not implemented yet
        
        % need to define ADC
        %             if(NewAcq())
        %                 lADCsToMeasure=4;
        %                 for lI = 1:lADCsToMeasure
        %                     sADC(lI).prep(6400, 1000); % prep(how many samples should be done, dwell time in nanoseconds)
        %                 end
        %             else
        %                 for lI = 1:lADCsToMeasure
        %                     sADC(lI).prep(1024, LL);
        %                 end
        %             end
        
        % ##### CAIPI #### %
%         if(spiral.UseCAIPI)
%             GradBuf (3,1:ulSpirSamples) = 0;
%             
%             GRAD_RASTER_TIME = 10;
%             spiral = DesignCAIPIBlips( paramLongSpGradAmp,paramLongSpSlewRate, ...
%             SeparationInmm,CAIPIDelay,CAIPIWantedInterval,ulSpirSamples*GRAD_RASTER_TIME,spiral);
%             
%             CAIPIGrad = spiral.caipiblip.Grad;
%             spiral.GradBuf(3,:) = CAIPIGrad(:);
%             spiral.CAIPIAmp = spiral.caipiblip.GradAmp;
%             
%         end
        % ##### END OF CAIPI #####%
        
%         spiral = CheckGradZerosCrossings(spiral.GradBuf(1,:),ulSpirSamples,spiral);
%         
%         FirstGradWaveLength = spiral.FirstGradWaveLength;
%         LastGradWaveLength = spiral.LastGradWaveLength;
%         
%         
%         GradDwellTime =spiral.GradDwellTime;
%         dMaxGradSpirAct = spiral.dMaxGradSpirAct;
%         
%         FirstGradWaveLength = FirstGradWaveLength *GradDwellTime;
%         LastGradWaveLength = LastGradWaveLength *GradDwellTime;
%         
%         %float MaxSlewRateXAct=GetMaxVecDiff(GradBuf[0],ulSpirSamples)*1e6*dMaxGradSpirAct/GradDwellTime;
%         %float MaxSlewRateYAct=GetMaxVecDiff(GradBuf[1],ulSpirSamples)*1e6*dMaxGradSpirAct/GradDwellTime;
%         
%         %float MaxSlewRateXActAfter20=GetMaxVecDiff(&(GradBuf[0][20]),ulSpirSamples-20)*1e6*dMaxGradSpirAct/GradDwellTime;
%         %float MaxSlewRateYActAfter20=GetMaxVecDiff(&(GradBuf[1][20]),ulSpirSamples-20)*1e6*dMaxGradSpirAct/GradDwellTime;
%         
%         
%         MaxSlewRateXAct=GetMaxVecDiff(spiral.GradBuf(1,:),ulSpirSamples)*1e6*dMaxGradSpirAct/GradDwellTime;
%         MaxSlewRateYAct=GetMaxVecDiff(spiral.GradBuf(2,:),ulSpirSamples)*1e6*dMaxGradSpirAct/GradDwellTime;
%         
%         MaxSlewRateXActAfter20=GetMaxVecDiff(spiral.GradBuf(1,21:end),ulSpirSamples-20)*1e6*dMaxGradSpirAct/GradDwellTime;
%         MaxSlewRateYActAfter20=GetMaxVecDiff(spiral.GradBuf(2,21:end),ulSpirSamples-20)*1e6*dMaxGradSpirAct/GradDwellTime;
%         
%         
%         MaxSlewRateZAct=0.0;
%         MaxSlewRateZActAfter20=0.0;
        
        
    end
    
end
%%%%%
if(PreparedGrad == 0)
    %error('the case of PreparedGrad equaling to zero is not impletemted yet, please refer to the commented part below to implement the code');
    %// Ramp down
    %for lI = 1: RAMP_DOWN_POINTS
    %    spiral.GradBuf(1,ulSpirSamples+lI) = spiral.GradBuf(1,ulSpirSamples) * (RAMP_DOWN_POINTS-lI) / (RAMP_DOWN_POINTS-1.0);
    %    spiral.GradBuf(2,ulSpirSamples+lI) = spiral.GradBuf(2,ulSpirSamples) * (RAMP_DOWN_POINTS-lI) / (RAMP_DOWN_POINTS-1.0);
    %    if size(spiral.GradBuf,1) == 3
    %        spiral.GradBuf(3,ulSpirSamples+lI) = spiral.GradBuf(3,ulSpirSamples) * (RAMP_DOWN_POINTS-lI) / (RAMP_DOWN_POINTS-1.0);
    %    end
    %end
    %
    % 		  % Find gradient sums for winders
    % 		  dGSumX=VecSum(GradBuf(1,:),ulSpirSamples + RAMP_DOWN_POINTS)*dGradientInt*dMaxGradSpirAct*1e9;
    % 		  dGSumY=VecSum(GradBuf(2,:),ulSpirSamples + RAMP_DOWN_POINTS)*dGradientInt*dMaxGradSpirAct*1e9;
    % 		  % For spiral in reflect the gradient shape and reverse sign
    % 		  if(lSpiralType==1)
    % 			  ReflectVec(GradBuf(1,:), ulSpirSamples+RAMP_DOWN_POINTS);
    % 			  MinusVec(GradBuf(1,:), ulSpirSamples+RAMP_DOWN_POINTS);
    % 			  ReflectVec(GradBuf(2,:), ulSpirSamples+RAMP_DOWN_POINTS);
    % 			  MinusVec(GradBuf(2,:), ulSpirSamples+RAMP_DOWN_POINTS);
    %           end
end

	 

