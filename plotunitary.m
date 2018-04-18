%-*-octave-*--
%
% plotunitary: plot traces of the solution of the schroedinger equaiton
%
% USAGE:
% 
% plotunitary()
%
% INPUT:
%
function usave = plotunitary(us, verbose)
  if nargin < 1
    printf("ERROR: no solution provided\n");
    return;
  end

  if nargin < 2
    verbose=0;
  end

  T=20;
  nsteps = length(us(1,1,:));
  
  N1 = length(us(:,1,1));
  N2 = length(us(1,:,1));
  if verbose==1
    printf("Data has dimensions %d x %d x %d\n", N1, N2, nsteps);
  end

  if (N1 != N2)
    printf("ERROR: N1=%d and N2=%d must be equal!\n");
    return;
  end

  t = linspace(0,T,nsteps);

# one figure for the response of each basis vector
  for q=1:N2
    figure(q);
    h=plot(t, abs(us(1,q,:)), t, abs(us(2,q,:)), t, abs(us(3,q,:)), t, abs(us(4,q,:)));
    set(h,"linewidth",2);
    legend("u0", "u1", "u2", "u3", "location", "north");
    tstr = sprintf("Resonse to initial data e%1d", q-1);
    title(tstr);
    xlabel("Time");
  end
    
end
