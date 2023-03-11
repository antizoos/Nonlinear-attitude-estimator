function dy=dr(y)
global sigmaR ko
dy=ko*CrossMatrix(sigmaR)*y;
end