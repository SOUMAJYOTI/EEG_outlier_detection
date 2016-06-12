function M=MAD(y)

M = median(abs(y-median(y)));