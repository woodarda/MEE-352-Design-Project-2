function [Cmin,Cmax] = FindCminCmax(Ch,Cc)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
if Ch<Cc
    Cmin=Ch;
    Cmax=Cc;
end

if Cc<Ch
    Cmin=Cc;
    Cmax=Ch;
end

end

