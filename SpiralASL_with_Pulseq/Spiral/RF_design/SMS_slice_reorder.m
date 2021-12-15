function lNewSliceIndex = SMS_slice_reorder(MrProt)

l_LongEntry_SmsFactor = MrProt.spcl.g4.SMS_factor ;

newSliceSize = MrProt.sliceSeries.size/l_LongEntry_SmsFactor;
%MrProt.sliceSeries.Slices = 'Interleaved';
for lSlice = 1:newSliceSize
    %lOrigSlcIndex = lSlice; %% needed to reset after EPI readout
	%%//recalculate the slice-index 
	lNewSliceIndex(lSlice)= getSmsRecalculatedSliceIndex( MrProt,  lSlice,  l_LongEntry_SmsFactor);
end
end

% benpos caipi --> helper function for slice index calculation
function lNewSliceIndex = getSmsRecalculatedSliceIndex(MrProt, lOrigSlcIndex, lSmsFactor)


	switch MrProt.sliceSeries.Slices

        case 'Ascending'
		
			lNewSliceIndex = lOrigSlcIndex; 
		
        case 'Descending'
		
			lNewSliceIndex = lOrigSlcIndex + (lSmsFactor-1) * MrProt.sliceSeries.size / lSmsFactor ;
		    lNewSliceIndex = MrProt.sliceSeries.size - lNewSliceIndex + 1;
            
        case 'Interleaved' % un peu bizzare parce que IDEA donne la simuluation tres bizzare
            lNewSliceIndex = lOrigSlcIndex; 
			if (mod(lNewSliceIndex,2) == 0)
                lNewSliceIndex = lNewSliceIndex/2+MrProt.sliceSeries.size / lSmsFactor / 2;
            else
                lNewSliceIndex = (lNewSliceIndex+1)/2;
            end
			%//interleaved always really means *interleaved ascending* order
			%//for an odd  number of slices it goes 0-2-4-6... 1-3-5     --> for an odd number slices the first "interleave" will have one more slice than the seccond
			%//for an even number of slices it goes 1-3-5-7... 0-2-4-6   --> for an even number of slices, the first and second "interleave" of slices will have an equal number 
			%// --> we can take the ceil value to work out the number of slices for the first "interleave" (true for both the of the full set of slices and slice-accelerated set of slices)
			%lNSlcFirstInterleaveFull = ceil( MrProt.sliceSeries.size/2  );
			%lNSlcFirstInterleaveSMS = ceil (MrProt.sliceSeries.size / lSmsFactor / 2 );
			%lAnatomicalSliceNo = pSLC->getSliceIndex();

			%// rejig the slice index if  the "anatomical range" and "slice index range" are exceeded
			%if ( ( lAnatomicalSliceNo >= pMrProt->sliceSeries().size()/lSmsFactor )  && (lOrigSlcIndex >= lNSlcFirstInterleaveSMS ) ) 
			%
			%	lNewSliceIndex = lOrigSlcIndex - lNSlcFirstInterleaveSMS + lNSlcFirstInterleaveFull  ;
			
			%else
			
			%	lNewSliceIndex = lOrigSlcIndex;
			
		
            %end
    end
		
end
 

