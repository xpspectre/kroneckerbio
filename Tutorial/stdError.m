function handle = stdError(times, stdData)

handle = @sdHandle;

    function [sigma dsigmady d2sigmady2] = sdHandle(t, yInd, yVal)
        sigma = interp1(times, stdData, t);
        
        if nargout >= 2
            % First derivative
            % TODO: calculate derivative - spread increases over time
            dsigmady = zeros(size(yVal));
            
            if nargout >=3
                % Second derivative
                d2sigmady2 = zeros(size(yVal));
            end
        end
    end
end