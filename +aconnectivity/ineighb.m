function n = ineighb(x,flag)
    % adjacent neighbours on a square matrix
    n = [ x - [1 0];
          x - [0 1];
          x + [1 0];
          x + [0 1] ];


    if nargin > 1 && flag
        n = [ x - [1 0];
              x - [1 1];
              x - [0 1];
              x + [1 0];
              x + [1 1];
              x + [0 1];
              x + [-1 1];
              x + [1 -1] ];
    end

end

