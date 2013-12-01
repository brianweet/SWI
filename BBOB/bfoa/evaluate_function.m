function [ fcurrent, fbest, xbest ] = evaluate_function( FUN, p, fbest, xbest )
%EVALUATE_FUNCTION Summary of this function goes here
%   Detailed explanation goes here
    fcurrent = feval(FUN, p);
    if fcurrent < fbest
        fbest = fcurrent;
        xbest = p;
    end
end

