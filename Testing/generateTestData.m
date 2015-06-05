function [outputsList, timesList, measurementsList] = generateTestData(m, expt, times, outputs, sd)
% Generate test data for sample fitting.

% Make times a col vector for convenience later
times = vec(times);

% Simulate system
sim = SimulateSystem(m, expt, times(end));

% Get output indices
outputNames = {m.Outputs.Name};
outputIdxs = find(ismember(outputNames, outputs));

% Generate data with noise
nOutputs = length(outputs);
nTimes = length(times);
measurements = zeros(nOutputs, nTimes);
for i = 1:nOutputs
    outputIdx = outputIdxs(i);
    for j = 1:nTimes
        time = times(j);
        val = sim.y(time, outputIdx);
        noise = normrnd(0, sd(time, outputIdx, val));
        measurements(i,j) = max(val + noise, 0); % ensure nonnegative value
    end
end

% Index and reshape outputs to form objective fun expects
outputsList = reshape(repmat(outputIdxs, nTimes, 1), nTimes*nOutputs, 1);
timesList = repmat(times, nOutputs, 1);
measurementsList = reshape(measurements',1, nTimes*nOutputs)';