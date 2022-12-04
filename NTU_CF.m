function [Cmin,Cmax,Cr,NTU] = NTU_CF(Eff,Ch,Cc,Ccmixed) %#ok<*STOUT> 
%NTU_CF Summary of this function goes here
%   Detailed explanation goes here

    if Cc<Ch 
        Cmin=Cc;
        Cmax=Ch;
    end

    if Ch<Cc
        Cmin=Ch;
        Cmax=Cc;
    end
    Cr=Cmin/Cmax;

    if Cc>Ch & Ccmixed==true | Cc<Ch & Ccmixed==false %#ok<*OR2,*AND2> 
        
        NTU=-log(1+log(1-Eff*Cr)/Cr);
    end
    if Cc<Ch & Ccmixed==true | Cc>Ch & Ccmixed==false
        NTU=-log(Cr*log(1-Eff)+1)/Cr;
    end
end

