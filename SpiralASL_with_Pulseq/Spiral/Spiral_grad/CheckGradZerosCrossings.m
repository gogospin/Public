function spiral = CheckGradZerosCrossings(Vec, N,spiral) 
	
%// find first zero crossing
for i=1:N-1
    if(Vec(i)*Vec(i+1)<0)
        break;
    end
end
FirstZC=i;
for i=FirstZC+1:N-1
    if (Vec(i)*Vec(i+1)<0)
        break;
    end
end
SecondZC=i;

%/ find last zero crossing
for i=N-1:-1:1
    if(Vec(i)*Vec(i+1)<0)
        break;
    end
end
LastZC=i;
for i=LastZC-1:-1:1
    if(Vec(i)*Vec(i+1)<0)
        break;
    end
end
SecondLastZC=i;

	spiral.FirstGradWaveLength=SecondZC-FirstZC;
	spiral.LastGradWaveLength=LastZC-SecondLastZC;
end