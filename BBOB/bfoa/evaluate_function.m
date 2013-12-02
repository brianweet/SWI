function [ fcurrent, fbest, xbest, sBestJ, sBestX] = evaluate_function( FUN, p, fbest, xbest, sBestJ, sBestX )
%EVALUATE_FUNCTION Summary of this function goes here
%   Detailed explanation goes here
    fcurrent = feval(FUN, p);
    if fcurrent < fbest
        fbest = fcurrent;
        xbest = p;
    end
    
    if fcurrent < sBestJ
        sBestJ = fcurrent;
        sBestX = p; 
    end
end

