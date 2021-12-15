function spiral = VDSpiral_sosp(FOV, res, AccR, nInterleaves, alpha, MaxGradAmp, MaxSlewRate, GradConnRatio, ExtraDelay,nsosp_segments,nz,FOVz,nsosp_rz,spiral)

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
%     nsosp_segments: number of segments to form one full spiral
%     nz: number of slices
%     FOVz: slab thickness
dGammaRad = spiral.dGammaRad;
dGamma = spiral.dGamma;
GRAD_RASTER_TIME = 10;
dGradientInt      = GRAD_RASTER_TIME * 1e-6;

% SOSP 
FOV_ori = FOV;
FOV = FOV_ori/nsosp_segments;
% SOSP

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

% SOSP rotate 
delta_angle_rot = 2*pi/(nsosp_segments);
TrajBuf0 = TrajBuf; clear TrajBuf;
for i = 1:nsosp_segments
    temp = complex(TrajBuf0(1,:),TrajBuf0(2,:));
    angle_rot = (i-1)*delta_angle_rot;
    temp = temp*exp(1i*angle_rot);
    TrajBuf(1,i,:) = real(temp);
    TrajBuf(2,i,:) = imag(temp);
end

clear GradBuf;
GradBuf = zeros(size(TrajBuf));
GradBuf(1,:,1)=0;
GradBuf(2,:,1)=0;
GradBuf(1,:,2:end)=TrajBuf(1,:,2:end)-TrajBuf(1,:,1:end-1);
GradBuf(2,:,2:end)=TrajBuf(2,:,2:end)-TrajBuf(2,:,1:end-1);
GradBuf = GradBuf/Dt/dGammaRad*1000; %in mT/m
   
   
% SOSP traveling across z
%CAIPIDelay = 200;
%CAIPIWantedInterval = 200;
%ulSpirSamples = 2560;
%grad_init = DesignSOSPBlips(FOVz/nz*2*1000,CAIPIDelay,CAIPIWantedInterval,ulSpirSamples*GRAD_RASTER_TIME,spiral);
%[gz_init,gz_init_Amp] = DesignSOSPBlips(FOVz/nz*2*1000, spiral ); % 
%[gz_delta,gz_delta_Amp] = DesignSOSPBlips(FOVz*1000, spiral );
%time_ent_slices = length(gz_delta);

MaxBlockMoment=dGammaRad*spiral.paramLongSpSlewRate* ...
    GRAD_RASTER_TIME*1e-6*GRAD_RASTER_TIME*1e-6; % in rad/m
MaxBlockGrad = spiral.paramLongSpSlewRate*GRAD_RASTER_TIME*1e-6*1000; % in mT/m
%dBlockNeeded=DesiredMoment/MaxBlockMoment;
   
%ClosestOdd=ceil(sqrt(dBlockNeeded)*2+1);
clear gx_refoc gy_refoc rampdown_grad_x rampdown_grad_y
for i = 1:nsosp_segments
    % part 1: design rampdown gradients after spiral gradients
    
    rampdown_sample_x(i) = ceil(abs(GradBuf(1,i,end)/MaxBlockGrad)); % positive if gradiend is changing from positive side
    rampdown_sample_y(i) = ceil(abs(GradBuf(2,i,end)/MaxBlockGrad)); %
    rampdown_sr_x(i) = GradBuf(1,i,end)/abs(rampdown_sample_x(i))/GRAD_RASTER_TIME; %mT/m/us
    rampdown_sr_y(i) = GradBuf(2,i,end)/abs(rampdown_sample_y(i))/GRAD_RASTER_TIME; %mT/m/us
    
    for j = 1:abs(rampdown_sample_x(i))
        rampdown_grad_x(i).grad(j) = GradBuf(1,i,end)-j*rampdown_sr_x(i)*GRAD_RASTER_TIME;
    end
    for j = 1:abs(rampdown_sample_y(i))
        rampdown_grad_y(i).grad(j) = GradBuf(2,i,end)-j*rampdown_sr_y(i)*GRAD_RASTER_TIME;
    end
    
    % part 2 design blips to bring kx-ky traj to zero
    
    Traj_refoc_x(i) = TrajBuf(1,i,end)+sum(rampdown_grad_x(i).grad(:))*GRAD_RASTER_TIME*dGammaRad*1e-9; % in rad/m
    Traj_refoc_y(i) = TrajBuf(2,i,end)+sum(rampdown_grad_y(i).grad(:))*GRAD_RASTER_TIME*dGammaRad*1e-9; % in rad/m
    
    [gx_refoc(i).grad,gx_refoc_amp(i)] = DesignSOSPBlips(abs(1000/Traj_refoc_x(i)*2*pi), spiral );
    [gy_refoc(i).grad,gy_refoc_amp(i)] = DesignSOSPBlips(abs(1000/Traj_refoc_y(i)*2*pi), spiral );
    if Traj_refoc_x(i) > 0
        gx_refoc(i).grad = -gx_refoc(i).grad;
        %gx_refoc_amp(i) = -gx_refoc_amp(i);
    end
    gx_refoc(i).grad = gx_refoc(i).grad*gx_refoc_amp(i);
    if Traj_refoc_y(i) > 0
        gy_refoc(i).grad = -gy_refoc(i).grad;
    end
    gy_refoc(i).grad = gy_refoc(i).grad*gy_refoc_amp(i);
    gx_refoc_time(i) = length(gx_refoc(i).grad);
    gy_refoc_time(i) = length(gy_refoc(i).grad);
    refoc_time(i) = max(gx_refoc_time(i),gy_refoc_time(i));
    rampdown_time(i) = max(rampdown_sample_x(i),rampdown_sample_y(i));
end
%refoc_time_max = max(refoc_time(:));
%rampdown_time_max = max(rampdown_time(:));

% organize kz direction sampling order:
kzfov_half = 1/(FOVz/nz)/2/1000; %rad/2pi/mm
delta_kz = 1/FOVz/1000; %rad/2pi/mmm
%[gz_init,gz_init_Amp] = DesignSOSPBlips(FOVz/nz*2*1000, spiral ); % 
%[gz_delta,gz_delta_Amp] = DesignSOSPBlips(FOVz*1000, spiral );
for i = 1:nsosp_rz
     kz_targ(i,1) = kzfov_half - (i-1)*delta_kz;
    for j = 2:nz/nsosp_rz
        kz_targ(i,j) = kzfov_half - (i-1)*delta_kz - nsosp_rz*delta_kz*(j-1);
    end
end
kz_switch_time = 0;
kz_switch_time_postrf = 0;
for i = 1:nsosp_rz
    [gz0_shot(i).kslice(1).grad(:),gz0_shot(i).kslice(1).amp] = DesignSOSPBlips(abs(1/kz_targ(i,1)), spiral );
    kz_switch_time_postrf = max(kz_switch_time_postrf,length(gz0_shot(i).kslice(1).grad(:)));
    if kz_targ(i,1) < 0
        gz0_shot(i).kslice(1).grad(:) = -gz0_shot(i).kslice(1).grad(:);
    end
    gz0_shot(i).kslice(1).grad(:) = gz0_shot(i).kslice(1).grad(:)*gz0_shot(i).kslice(1).amp;
    for j = 2:nz/nsosp_rz
        [gz0_shot(i).kslice(j).grad(:),gz0_shot(i).kslice(j).amp] = DesignSOSPBlips(abs(1/(kz_targ(i,j)-kz_targ(i,j-1))), spiral );
        if kz_targ(i,j)-kz_targ(i,j-1)<0
            gz0_shot(i).kslice(j).grad(:) = -gz0_shot(i).kslice(j).grad(:);
        end
        kz_switch_time = max(kz_switch_time,length(gz0_shot(i).kslice(j).grad(:)));
        gz0_shot(i).kslice(j).grad(:) = gz0_shot(i).kslice(j).grad(:)*gz0_shot(i).kslice(j).amp;
    end
    
end
        
% build the entire gradient readout
refoc_time_xy = max(refoc_time(:));
rampdown_time_xy = max(rampdown_time(:));

length_interv(1) = kz_switch_time_postrf;
length_interv(2) = size(GradBuf,3);
length_interv(3) = max(kz_switch_time,rampdown_time_xy+refoc_time_xy);
end_interv(1) = length_interv(1);
end_interv(2) = length_interv(1)+length_interv(2);
end_interv(3) = length_interv(1)+length_interv(2)+length_interv(3);

RAMP_DOWN_POINTS = 50;   % Ramp down of the gradients at the end of a spiral gradient
spiral.RAMP_DOWN_POINTS = RAMP_DOWN_POINTS;

for k = 1:nsosp_segments
    for i = 1:nsosp_rz
        gz_seg(k).shot(i).kslice(1).grad(1:length(gz0_shot(i).kslice(1).grad)) = -gz0_shot(i).kslice(1).grad;
        gz_seg(k).shot(i).kslice(1).grad(end_interv(1)+1:end_interv(2)) = 0;
        gx_seg(k).shot(i).kslice(1).grad(1:end_interv(1)) = 0;
        gx_seg(k).shot(i).kslice(1).grad(1+end_interv(1):end_interv(2)) = GradBuf(1,k,:);
        gy_seg(k).shot(i).kslice(1).grad(1:end_interv(1)) = 0;
        gy_seg(k).shot(i).kslice(1).grad(1+end_interv(1):end_interv(2)) = GradBuf(2,k,:);
        adc_shot(i).start(1) = 1+end_interv(1);
        adc_shot(i).end(1) = end_interv(2);
        adc_shot(i).end(1) = round(adc_shot(i).end(1)/10)*10;
        for j = 2:nz/nsosp_rz
            gz_seg(k).shot(i).kslice(j).grad(1:length(gz0_shot(i).kslice(j).grad)) = -gz0_shot(i).kslice(j).grad;
            gz_seg(k).shot(i).kslice(j).grad(length_interv(3)+1:length_interv(2)+length_interv(3)) = 0;
            gx_seg(k).shot(i).kslice(j).grad = zeros(length_interv(2)+length_interv(3),1);
            gx_seg(k).shot(i).kslice(j).grad(1:length(rampdown_grad_x(k).grad)) = rampdown_grad_x(k).grad;
            gx_seg(k).shot(i).kslice(j).grad(rampdown_time_xy+1:rampdown_time_xy+length(gx_refoc(k).grad)) = gx_refoc(k).grad;
            gx_seg(k).shot(i).kslice(j).grad(length_interv(3)+1:length_interv(2)+length_interv(3)) = GradBuf(1,k,:);
            gy_seg(k).shot(i).kslice(j).grad = zeros(length_interv(2)+length_interv(3),1);
            gy_seg(k).shot(i).kslice(j).grad(1:length(rampdown_grad_y(k).grad)) = rampdown_grad_y(k).grad;
            gy_seg(k).shot(i).kslice(j).grad(rampdown_time_xy+1:rampdown_time_xy+length(gy_refoc(k).grad)) = gy_refoc(k).grad;
            gy_seg(k).shot(i).kslice(j).grad(length_interv(3)+1:length_interv(2)+length_interv(3)) = GradBuf(2,k,:);
            adc_shot(i).start(j) = 1+end_interv(2)+(j-2)*(length_interv(3)+length_interv(2))+length_interv(3);
            adc_shot(i).end(j) = end_interv(2)+(j-1)*(length_interv(3)+length_interv(2));            
        end
        for lI = 1: RAMP_DOWN_POINTS
            gx_seg(k).shot(i).rampdown(lI) = gx_seg(k).shot(i).kslice(end).grad(end) * (RAMP_DOWN_POINTS-lI) / (RAMP_DOWN_POINTS-1.0);
            gy_seg(k).shot(i).rampdown(lI) = gy_seg(k).shot(i).kslice(end).grad(end) * (RAMP_DOWN_POINTS-lI) / (RAMP_DOWN_POINTS-1.0);
            gz_seg(k).shot(i).rampdown(lI) = gz_seg(k).shot(i).kslice(end).grad(end) * (RAMP_DOWN_POINTS-lI) / (RAMP_DOWN_POINTS-1.0);
        end
    end
end

%test
% for i = 1:3
%     counti = 1;
%     for j = 1:4
%         for k = 1:length(gx_shot(i).seg(j).grad)
%             if counti == 1
%                 kx_evol(i,counti) = gx_shot(i).seg(j).grad(k);
%             else
%                 kx_evol(i,counti) = kx_evol(i,counti-1)+gx_shot(i).seg(j).grad(k);
%             end
%             counti = counti+1;
%         end
%     end
% end
        



   spiral.gx_seg = gx_seg;
   spiral.gy_seg = gy_seg;
   spiral.gz_seg = gz_seg;
   spiral.adc_shot = adc_shot;
   spiral.kz_switch_time_postrf = kz_switch_time_postrf;
   spiral.kz_switch_time = kz_switch_time;
   spiral.refoc_time_xy = refoc_time_xy;
   spiral.rampdown_time_xy = rampdown_time_xy;
   spiral.TrajBuf = TrajBuf;
   spiral.ulSpirSamples = size(GradBuf,3);
   %spiral.dMaxGradSpirAct = dMaxGradSpirAct;
   

end

