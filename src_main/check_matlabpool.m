%% check matlabpool for backward compatibility of MATLAB versions older than 2014!

if verLessThan('matlab', '8.3.0')
    isOpen = (matlabpool('size') > 0); %#ok<*DPOOL>
    if (~isOpen)
        mycluster = parcluster;
        matlabpool('local',mycluster.NumWorkers);
    end
end