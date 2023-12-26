function [departure, flyby, arrival] = loadWindows()

departure{1}.lb = 1.15e4;
departure{1}.ub = 1.42e4;
flyby{1}.lb = 1.3e4;
flyby{1}.ub = 1.52e4;
arrival{1}.lb = 1.53e4;
arrival{1}.ub = 1.57e4;

departure{2}.lb = 1.15e4;
departure{2}.ub = 1.42e4;
flyby{2}.lb = 1.27e4;
flyby{2}.ub = 1.46e4;
arrival{2}.lb = 1.47e4;
arrival{2}.ub = 1.51e4;

departure{3}.lb = 1.15e4;
departure{3}.ub = 1.42e4;
flyby{3}.lb = 1.29e4;
flyby{3}.ub = 1.51e4;
arrival{3}.lb = 1.53e4;
arrival{3}.ub = 1.58e4;

departure{4}.lb = 1.15e4;
departure{4}.ub = 1.42e4;
flyby{4}.lb = 1.29e4;
flyby{4}.ub = 1.55e4;
arrival{4}.lb = 1.59e4;
arrival{4}.ub = 1.64e4;

departure{5}.lb = 1.15e4;
departure{5}.ub = 1.42e4;
flyby{5}.lb = 1.29e4;
flyby{5}.ub = 1.56e4;
arrival{5}.lb = 1.65e4;
arrival{5}.ub = 1.71e4;

departure{6}.lb = 1.15e4;
departure{6}.ub = 1.42e4;
flyby{6}.lb = 1.3e4;
flyby{6}.ub = 1.55e4;
arrival{6}.lb = 1.71e4;
arrival{6}.ub = 1.76e4;

departure{7}.lb = 1.15e4;
departure{7}.ub = 1.42e4;
flyby{7}.lb = 1.3e4;
flyby{7}.ub = 1.55e4;
arrival{7}.lb = 1.77e4;
arrival{7}.ub = 1.82e4;

departure{8}.lb = 1.15e4;
departure{8}.ub = 1.42e4;
flyby{8}.lb = 1.3e4;
flyby{8}.ub = 1.55e4;
arrival{8}.lb = 1.83e4;
arrival{8}.ub = 1.88e4;

departure{9}.lb = 1.15e4;
departure{9}.ub = 1.42e4;
flyby{9}.lb = 1.32e4;
flyby{9}.ub = 1.55e4;
arrival{9}.lb = 1.89e4;
arrival{9}.ub = 1.94e4;

departure{10}.lb = 1.15e4;
departure{10}.ub = 1.42e4;
flyby{10}.lb = 1.33e4;
flyby{10}.ub = 1.55e4;
arrival{10}.lb = 1.95e4;
arrival{10}.ub = 2e4;

departure{11}.lb = 1.15e4;
departure{11}.ub = 1.42e4;
flyby{11}.lb = 1.32e4;
flyby{11}.ub = 1.54e4;
arrival{11}.lb = 2.01e4;
arrival{11}.ub = 2.06e4;

departure{12}.lb = 1.15e4;
departure{12}.ub = 1.42e4;
flyby{12}.lb = 1.31e4;
flyby{12}.ub = 1.56e4;
arrival{12}.lb = 2.07e4;
arrival{12}.ub = 2.12e4;

% new column
departure{13}.lb = 1.15e4;
departure{13}.ub = 1.42e4;
flyby{13}.lb = 1.54e4;
flyby{13}.ub = 1.66e4;
arrival{13}.lb = 1.66e4;
arrival{13}.ub = 1.71e4;

departure{14}.lb = 1.15e4;
departure{14}.ub = 1.42e4;
flyby{14}.lb = 1.55e4;
flyby{14}.ub = 1.73e4;
arrival{14}.lb = 1.72e4;
arrival{14}.ub = 1.77e4;

departure{15}.lb = 1.15e4;
departure{15}.ub = 1.42e4;
flyby{15}.lb = 1.55e4;
flyby{15}.ub = 1.8e4;
arrival{15}.lb = 1.79e4;
arrival{15}.ub = 1.83e4;

departure{16}.lb = 1.15e4;
departure{16}.ub = 1.42e4;
flyby{16}.lb = 1.56e4;
flyby{16}.ub = 1.84e4;
arrival{16}.lb = 1.86e4;
arrival{16}.ub = 1.89e4;

departure{17}.lb = 1.15e4;
departure{17}.ub = 1.42e4;
flyby{17}.lb = 1.57e4;
flyby{17}.ub = 1.79e4;
arrival{17}.lb = 1.92e4;
arrival{17}.ub = 1.96e4;

departure{18}.lb = 1.15e4;
departure{18}.ub = 1.42e4;
flyby{18}.lb = 1.57e4;
flyby{18}.ub = 1.79e4;
arrival{18}.lb = 1.98e4;
arrival{18}.ub = 2.01e4;

departure{19}.lb = 1.15e4;
departure{19}.ub = 1.42e4;
flyby{19}.lb = 1.58e4;
flyby{19}.ub = 1.78e4;
arrival{19}.lb = 2.04e4;
arrival{19}.ub = 2.07e4;

% new column
departure{20}.lb = 1.15e4;
departure{20}.ub = 1.42e4;
flyby{20}.lb = 1.7e4;
flyby{20}.ub = 1.82e4;
arrival{20}.lb = 1.83e4;
arrival{20}.ub = 1.86e4;

departure{21}.lb = 1.15e4;
departure{21}.ub = 1.42e4;
flyby{21}.lb = 1.71e4;
flyby{21}.ub = 1.88e4;
arrival{21}.lb = 1.89e4;
arrival{21}.ub = 1.93e4;

departure{22}.lb = 1.15e4;
departure{22}.ub = 1.42e4;
flyby{22}.lb = 1.72e4;
flyby{22}.ub = 1.93e4;
arrival{22}.lb = 1.95e4;
arrival{22}.ub = 1.99e4;

departure{23}.lb = 1.15e4;
departure{23}.ub = 1.42e4;
flyby{23}.lb = 1.75e4;
flyby{23}.ub = 1.96e4;
arrival{23}.lb = 2.01e4;
arrival{23}.ub = 2.05e4;

departure{24}.lb = 1.15e4;
departure{24}.ub = 1.42e4;
flyby{24}.lb = 1.73e4;
flyby{24}.ub = 1.99e4;
arrival{24}.lb = 2.08e4;
arrival{24}.ub = 2.12e4;

end

