function [estimates, model,sse] = Hill_fit(concentrations,Frequency_response)


% Call fminsearch with a random starting point.;
start_point = [max(Frequency_response)-min(Frequency_response) concentrations(floor(length(concentrations)/2)) min(Frequency_response)];
model = @Hill_curve;
estimates = fminsearch(model, start_point);
[sse, FittedCurve] = Hill_curve(estimates);

% Hill_curve accepts curve parameters as inputs, and outputs sse,
% the sum of squares error for A * exp(-lambda * xdata) - ydata, 
% and the FittedCurve. FMINSEARCH only needs sse, but we want to 
% plot the FittedCurve at the end.
    function [sse, FittedCurve] = Hill_curve(params)
        
        
        length (concentrations);
        A = params(1);
        c0 = params(2);
        b=params(3);
        FittedCurve = (concentrations)./(concentrations+c0);
        ErrorVector = FittedCurve - Frequency_response;
        sse = sum(ErrorVector .^ 2);
    end
end
