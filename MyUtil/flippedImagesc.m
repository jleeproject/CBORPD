function img = flippedImagesc(domain1, domain2, data)
    img = imagesc((domain1), domain2, (data));
    set(gca,'YDir','normal');
end