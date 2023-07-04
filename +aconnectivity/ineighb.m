function n = ineighb(x)
    % adjacent neighbours on a square matrix
    n = [ x - [1 0];
          x - [1 1];
          x - [0 1];
          x + [1 0];
          x + [1 1];
          x + [0 1];
          x + [-1 1];
          x + [1 -1] ];
    

end