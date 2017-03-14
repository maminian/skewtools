function progressbar(idx,last,steps)

checkpoints = linspace(1,last,steps);

if any(idx==checkpoints)
     pct=idx/last*100;
     
     disp( strcat(num2str(pct), '%') )
end

end
