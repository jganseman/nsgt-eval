function definput=arg_firwin(definput)
  
  definput.flags.wintype={'hanning','hann','sqrthan','sqrthann','hamming', ...
                      'sqrtham','square','halfsquare','rect', ...
                      'sqrtsquare','sqrtrect', 'tria','triangular', ...
                      'sqrttria','blackman','nuttall', 'ogg','itersine', ...
                      'sine'};
