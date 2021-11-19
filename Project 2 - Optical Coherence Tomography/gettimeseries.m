function [subPixelMotion, time] = gettimeseries(timeSeries, deltaT, pixel, lambda0, n)
    subPixelMotion = lambda0 / (4*pi*n) * unwrap(angle(timeSeries(pixel, :)));
    time = 0:deltaT:(size(timeSeries, 2)-1)*deltaT;
end