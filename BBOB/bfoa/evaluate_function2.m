function [ fcurrent, fbest, xbest ] = evaluate_function2( FUN, p, fbest, xbest )
%EVALUATE_FUNCTION Summary of this function goes here
%   Detailed explanation goes here
    fcurrent = feval(FUN, p);
    if fcurrent < fbest
        fbest = fcurrent;
        xbest = p;
    end
end
