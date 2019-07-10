function diagMat=UpdateDiagonal(diagMat, num, alpha, step);

diagMat(1:(num+1):end)=nan;
diagstd=nanstd(diagMat,1);
diagMat(1:(num+1):end)=diagstd*num*exp(2*alpha*step);
