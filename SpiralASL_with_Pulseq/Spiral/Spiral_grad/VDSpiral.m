function spiral = VDSpiral(FOV, res, AccR, nInterleaves, alpha, MaxGradAmp, MaxSlewRate, GradConnRatio, ExtraDelay,spiral)

%	field of view in meters
%   resolution (Nyquist distance) in meters
%   undersampling factor along frequency encoding direction (normalized units, positive, may be lower than one)
%   undersampling factor (normalized units, positive, may be lower than one)
%   number of interleaves
%   variable density factor (alpha in [1])
%   maximum amplitude (in T/m)
%   maximum slew rate (in T/m/s)*/

%     //double gamma = 2.678e8; // in rad/T/s
%     /*double FOV=7;//y[0];
%     double res=8;//y[1];
%     double f_sampling=3;//y[2];
%     double R=1;//y[3];
%     double nInterleaves=1;//y[4];
%     double alpha=2;//y[5];
%     double gm=2;//y[6];
%     double sm=1;//y[7];*/
dGammaRad = spiral.dGammaRad;
dGamma = spiral.dGamma;
GRAD_RASTER_TIME = 10;
dGradientInt      = GRAD_RASTER_TIME * 1e-6;


FacA=((GRAD_RASTER_TIME*1e-6)*dGammaRad)/1000;
dkToGFac=FOV/2/pi();

lambda = .5/res; % in m^(-1)
n = 1/(1-pow((1-nInterleaves*AccR/(FOV*lambda)),(1/alpha)));

w = 2*pi()*n;
Tea = lambda*w/(dGammaRad*MaxGradAmp*(alpha+1)); % in s
Tes = sqrt(lambda*(w*w)/(MaxSlewRate*dGammaRad))/(alpha/2+1); % in s
Ts2a = pow(pow(Tes,((alpha+1)/(alpha/2+1)))*(alpha/2+1)/(Tea*(alpha+1)),(1+2/alpha)); %in s

T2saSmallerThanTes = Ts2a<Tes;

tautrans=0;
dt=Tea*(1e-4); % in s
if(Ts2a<Tes)
    tautrans = pow(Ts2a/Tes,(1/(alpha/2+1)));
    Tend = Tea;
else
    Tend = Tes;
    
end

Dt=dGradientInt;

%VecLen=Tend/Dt +1; % cpp version
VecLen=ceil(Tend/Dt); % matlab version
if mod(VecLen,10)~=0
    VecLen=round(VecLen/10)*10;
end
VecLenm1=VecLen-1;

if(Ts2a<Tes)
    for i=0:VecLen-1
        Curt=i*Dt;
        
        tmp=kFunc1(lambda,w,Curt,Tes,alpha,Ts2a,Tea,tautrans);
        
        curw=cmuld(tmp,2*pi());
        TrajBuf(1,i+1)=curw.real;
        TrajBuf(2,i+1)=curw.imag;
        
    end
else
    for i=0:VecLen-1
        Curt=i*Dt;
        
        tmp=kFunc2(lambda,w,Curt,Tes,alpha,Ts2a,Tea,tautrans);
        
        curw=cmuld(tmp,2*pi());
        TrajBuf(1,i+1)=curw.real;
        TrajBuf(2,i+1)=curw.imag;
    end
end
    
for i=VecLen-1:-1:ExtraDelay
    TrajBuf(1,1+i)=TrajBuf(1,i-ExtraDelay+1);
    TrajBuf(2,1+i)=TrajBuf(2,i-ExtraDelay+1);
end
for i=1:ExtraDelay
    TrajBuf(1,i)=0;
    TrajBuf(2,i)=0;
end

   % the case of spiral.OutInx>1 is not included here
   % pleaes see the cpp code for the implementation
   
   if(spiral.OutInOut || spiral.OutIn)
       
       if(spiral.OutInOut)
           for i=1:VecLen
               TrajBuf(1,i+2*VecLen)=TrajBuf(1,i);
               TrajBuf(2,i+2*VecLen)=TrajBuf(2,i);
           end
       end
       
       
       RefI=VecLen-1;
       RefI = RefI+1; % from C to matlab
       CurPhi=atan2(TrajBuf(2,RefI),TrajBuf(1,RefI));
       CurR=sqrt(TrajBuf(1,RefI)*TrajBuf(1,RefI)+TrajBuf(2,RefI)*TrajBuf(2,RefI));
       RefPhi=atan2(TrajBuf(2,RefI),TrajBuf(1,RefI));
       RefPhi2=atan2(TrajBuf(2,RefI-1),TrajBuf(1,RefI-1));
       dPhi=RefPhi-RefPhi2;
       CurPhi = CurPhi+dPhi;
       
       
       
       TrajBuf(1,VecLen+1)=CurR*cos(CurPhi);
       TrajBuf(2,VecLen+1)=CurR*sin(CurPhi);
       
       for i=1:VecLen-1
            
           RefI = VecLen-1-i;
           RefI = RefI+1; % from C to matlab
           CurR=sqrt(TrajBuf(1,RefI)*TrajBuf(1,RefI)+TrajBuf(2,RefI)*TrajBuf(2,RefI));
           RefPhi=atan2(TrajBuf(2,RefI+1),TrajBuf(1,RefI+1));
           RefPhi2=atan2(TrajBuf(2,RefI),TrajBuf(1,RefI));
           dPhi=RefPhi-RefPhi2;
           CurPhi = CurPhi+dPhi;
           
           TrajBuf(1,1+i+VecLen)=CurR*cos(CurPhi);
           TrajBuf(2,1+i+VecLen)=CurR*sin(CurPhi);
       end
       
       if(spiral.OutInOut)
           VecLen=VecLen*3;
       end
       if(spiral.OutIn)
           VecLen=VecLen*2;
       end
   end
   
   GradBuf(1,1)=0;
   GradBuf(2,1)=0;
   % now diff them
   for i=0:VecLen-2
       GradBuf(1,i+1+1)=TrajBuf(1,i+1+1)-TrajBuf(1,i+1);
       GradBuf(2,i+1+1)=TrajBuf(2,i+1+1)-TrajBuf(2,i+1);
   end
   
   GradBuf(1,:) = VecMul(GradBuf(1,:),VecLen,1/(Dt*dGammaRad));
   GradBuf(2,:) = VecMul(GradBuf(2,:),VecLen,1/(Dt*dGammaRad));
   
   ulSpirSamples = VecLen;
   
   % Find the real maximum
   dMaxGradSpirAct=max(VecAbsMax(GradBuf(1,:),ulSpirSamples),VecAbsMax(GradBuf(2,:),ulSpirSamples));
   % Scale to one
   GradBuf(1,:) = VecMul(GradBuf(1,:),ulSpirSamples,1/dMaxGradSpirAct);
   GradBuf(2,:) = VecMul(GradBuf(2,:),ulSpirSamples,1/dMaxGradSpirAct);
   
   
   spiral.GradBuf = GradBuf;
   spiral.TrajBuf = TrajBuf;
   spiral.ulSpirSamples = ulSpirSamples;
   spiral.dMaxGradSpirAct = dMaxGradSpirAct;
   

end

