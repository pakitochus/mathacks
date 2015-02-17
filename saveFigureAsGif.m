function saveFigureAsGif(gifname, seg)
% Saves an animated 3D figure as a GIF 
fig=gcf;
h=gca;
gif_fps = 24;
nframes=gif_fps*seg;
paso=round(360/nframes);
initPos=get(h,'View');
pos=initPos;
for k=1:nframes
    frame = getframe(fig);
    if k == 1
        [animated, cmap] = rgb2ind(frame.cdata, 256, 'nodither');
    else
        animated(:,:,1,k) = rgb2ind(frame.cdata, cmap, 'nodither');
    end
    pos(1)=pos(1)+paso;
    set(h,'View',pos);
end
imwrite(animated, cmap, gifname, 'DelayTime', 1/gif_fps, ...
   'LoopCount', inf);
web(gifname)

end
