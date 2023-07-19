function aimagesc(x)

% wrapper on imagesc that produces a symmetric colorbar
s = max(abs(x(:)))*1.1;

imagesc(x);caxis([-s s]);