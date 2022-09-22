function recon = CG_Recon(AhA,Ahd,Niter,tol,Display)

x = zeros(size(Ahd));
r = Ahd - AhA(x);
p = r;
rsold = r(:)'*r(:);

for k = 1:Niter
    z = AhA(p);
    alpha = rsold / (p(:)'*z(:));
    x = x + alpha * p;
    r = r - alpha *z;
    rsnew = r(:)'*r(:);
    if sqrt(rsnew) < tol
        break;
    end
    p = r + (rsnew / rsold) * p;
    rsold = rsnew;
 
    if Display
    drawnow,
    imagesc(abs(x(:,:))),colormap gray,axis image off,title(['iter-' num2str(k)])    
    end
end

recon = x;