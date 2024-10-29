function perform_optimization()
global OPT

switch OPT.options.optimizer
    case 'fmincon-active-set'
        OPT.history = runfmincon(OPT.dv,@(x)obj(x),@(x)nonlcon(x));
    case 'mma'
        OPT.history = runmma(OPT.dv,@(x)obj(x),@(x)nonlcon(x));
end
