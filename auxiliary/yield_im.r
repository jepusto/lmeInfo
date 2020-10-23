yim = function(params)
{
  iset = params$iset; H = length(iset);
  dat = params$dat[iset]
  nmu = params$nmu; ivc = params$ivc[iset];
  vc = params$vc[iset]; mode = params$mode[iset];

  nvc = numeric(H)
  for(ii in 1:H){ nvc[ii] = length(vc[[ii]]) }

  I11 = list(H)
  I12 = list(H)
  I22 = list(H)
  I33 = list(H)
  I2g1 = list(H)
  for(ii in 1:H){
    I11[[ii]] = matrix(0,nmu,nmu)
    I12[[ii]] = matrix(0,nmu,1)
    I22[[ii]] = 0
    I33[[ii]] = matrix(0,nvc[ii],nvc[ii])
    if(ivc[ii] != "Path"){
      Source = factor(dat[[ii]]$Source)
      SourceLevels = levels(Source)
      nslev = length(SourceLevels)
      isvc = vc[[ii]]
      te = ts1 = ts12 = tp1 = tp12 = tp13 = 0
      if(mode[ii] == "REML"){
        U = matrix(0,nmu+1,nmu+1) 
        Te1 = Ts1 = Ts12 = Tp1 = Tp12 = Tp13 = matrix(0,nmu+1,nmu+1)
        Te2 = Ts2 = Tp2 = matrix(0,nmu+1,nmu+1)
      }
      for(jj in 1:nslev){
        tsdat = dat[[ii]][dat[[ii]]$Source == SourceLevels[jj],]
        ntsdat = nrow(tsdat)
        X = cbind(rep(1,ntsdat), tsdat$Log10.Delta.)
        if(ivc[ii] == "Full"){
          Path = factor(tsdat$Path)
          n = as.numeric(table(Path))
          ln = length(n)
          n = c(1,n)
          cn = cumsum(n)
          Z2 = matrix(0,ntsdat,ln)
          for(kk in 1:ln){
            Z2[cn[kk]:(cn[kk+1]-1),kk] = 1
          }
        }
        if(ivc[ii] != "Residual"){
          Omega = matrix(isvc[1],ntsdat,ntsdat)
        } else {
          Omega = matrix(0,ntsdat,ntsdat)
        }
        if (ivc[ii] == "Full") {
          Omega = Omega + isvc[2] * Z2 %*% t(Z2)
        }
        for(kk in 1:ntsdat){
          Omega[kk,kk] = Omega[kk,kk] + isvc[nvc[ii]]
        }
        cOmega = chol(Omega)
        Xt = backsolve(cOmega,X,transpose=TRUE)
        yt = backsolve(cOmega,rep(1,ntsdat),transpose=TRUE)
        yt2 = t(yt) %*% yt
        
        #
        # Information matrix fixed effects
        #
        I11[[ii]] = I11[[ii]] + t(Xt) %*% Xt
        I12[[ii]] = I12[[ii]] + tsdat$Log10.Yield.[1] * (t(Xt) %*% yt)
        I22[[ii]] = I22[[ii]] + (tsdat$Log10.Yield.[1])^2 * yt2

        dOmegaR = diag(ntsdat)
        dOmegaRt = backsolve(cOmega,dOmegaR,transpose=TRUE)
        dOmegaRtt = backsolve(cOmega,dOmegaRt)
        te = te + sum(diag(dOmegaRtt %*% dOmegaRtt))
        if(ivc[ii] != "Residual"){
          dOmegaS = matrix(1,ntsdat,ntsdat)
          dOmegaSt = backsolve(cOmega,dOmegaS,transpose=TRUE)
          dOmegaStt = backsolve(cOmega,dOmegaSt)
          ts1 = ts1 + sum(diag(dOmegaStt %*% dOmegaStt))
          ts12 = ts12 + sum(diag(dOmegaStt %*% dOmegaRtt))
        }
        if(ivc[ii] == "Full"){
          dOmegaP = Z2 %*% t(Z2)
          dOmegaPt = backsolve(cOmega,dOmegaP,transpose=TRUE)
          dOmegaPtt = backsolve(cOmega,dOmegaPt)
          tp1 = tp1 + sum(diag(dOmegaPtt %*% dOmegaPtt))
          tp12 = tp12 + sum(diag(dOmegaStt %*% dOmegaPtt))
          tp13 = tp13 + sum(diag(dOmegaPtt %*% dOmegaRtt))
        }
        if(mode[ii] == "REML"){
          C = cbind(X,rep(1,ntsdat)*tsdat$Log10.Yield.[1])
          Ct = cbind(Xt,yt*tsdat$Log10.Yield.[1])
          Ctt = backsolve(cOmega,Ct)
          U = U + t(Ct) %*% Ct
          Te1 = Te1 + t(C) %*% dOmegaRtt %*% dOmegaRtt %*% Ctt
          Te2 = Te2 + t(C) %*% dOmegaRtt %*% Ctt
          if(ivc[ii] != "Residual"){
            Ts1 = Ts1 + t(C) %*% dOmegaStt %*% dOmegaStt %*% Ctt
            Ts12 = Ts12 + t(C) %*% dOmegaStt %*% dOmegaRtt %*% Ctt
            Ts2 = Ts2 + t(C) %*% dOmegaStt %*% Ctt
          }
          if(ivc[ii] == "Full"){
            Tp1 = Tp1 + t(C) %*% dOmegaPtt %*% dOmegaPtt %*% Ctt
            Tp12 = Tp12 + t(C) %*% dOmegaStt %*% dOmegaPtt %*% Ctt
            Tp13 = Tp13 + t(C) %*% dOmegaPtt %*% dOmegaRtt %*% Ctt
            Tp2 = Tp2 + t(C) %*% dOmegaPtt %*% Ctt
          }
        }
      }
      if(mode[ii] == "REML"){
        cU = chol(U)
        Te1t = backsolve(cU,Te1,transpose=TRUE)
        Te1tt = sum(diag(backsolve(cU,Te1t)))
        Te2t = backsolve(cU,Te2,transpose=TRUE)
        Te2tt = backsolve(cU,Te2t)
        if(ivc[ii] != "Residual"){
          Ts1t = backsolve(cU,Ts1,transpose=TRUE)
          Ts1tt = sum(diag(backsolve(cU,Ts1t)))
          Ts12t = backsolve(cU,Ts12,transpose=TRUE)
          Ts12tt = sum(diag(backsolve(cU,Ts12t)))
          Ts2t = backsolve(cU,Ts2,transpose=TRUE)
          Ts2tt = backsolve(cU,Ts2t)
        }
        if(ivc[ii] == "Full"){
          Tp1t = backsolve(cU,Tp1,transpose=TRUE)
          Tp1tt = sum(diag(backsolve(cU,Tp1t)))
          Tp12t = backsolve(cU,Tp12,transpose=TRUE)
          Tp12tt = sum(diag(backsolve(cU,Tp12t)))
          Tp13t = backsolve(cU,Tp13,transpose=TRUE)
          Tp13tt = sum(diag(backsolve(cU,Tp13t)))
          Tp2t = backsolve(cU,Tp2,transpose=TRUE)
          Tp2tt = backsolve(cU,Tp2t)
        }
      }

      #
      # Information matrix variance components
      #
      if(ivc[ii] == "Residual"){
        I33[[ii]][1,1] = te
      } else if(ivc[ii] == "Source"){
        I33[[ii]][1,1] = ts1
        I33[[ii]][1,2] = ts12
        I33[[ii]][2,1] = I33[[ii]][1,2]
        I33[[ii]][2,2] = te
      } else {
        I33[[ii]][1,1] = ts1
        I33[[ii]][1,2] = tp12
        I33[[ii]][2,1] = I33[[ii]][1,2]
        I33[[ii]][1,3] = ts12
        I33[[ii]][3,1] = I33[[ii]][1,3]
        I33[[ii]][2,2] = tp1
        I33[[ii]][2,3] = tp13
        I33[[ii]][3,2] = I33[[ii]][2,3]
        I33[[ii]][3,3] = te
      }
      if(mode[ii] == "REML"){
        if(ivc[ii] == "Residual"){
          I33[[ii]][1,1] = I33[[ii]][1,1] + sum(diag(Te2tt %*% Te2tt)) -
                           2*Te1tt
        } else if(ivc[ii] == "Source"){
          I33[[ii]][1,1] = I33[[ii]][1,1] + sum(diag(Ts2tt %*% Ts2tt)) -
                           2*Ts1tt
          I33[[ii]][1,2] = I33[[ii]][1,2] + sum(diag(Ts2tt %*% Te2tt)) -
                           2*Ts12tt
          I33[[ii]][2,1] = I33[[ii]][1,2]
          I33[[ii]][2,2] = I33[[ii]][2,2] + sum(diag(Te2tt %*% Te2tt)) -
                           2*Te1tt
        } else {
          I33[[ii]][1,1] = I33[[ii]][1,1] + sum(diag(Ts2tt %*% Ts2tt)) -
                           2*Ts1tt
          I33[[ii]][1,2] = I33[[ii]][1,2] + sum(diag(Ts2tt %*% Tp2tt)) -
                           2*Tp12tt
          I33[[ii]][2,1] = I33[[ii]][1,2]
          I33[[ii]][1,3] = I33[[ii]][1,3] + sum(diag(Ts2tt %*% Te2tt)) -
                           2*Ts12tt
          I33[[ii]][3,1] = I33[[ii]][1,3]
          I33[[ii]][2,2] = I33[[ii]][2,2] + sum(diag(Tp2tt %*% Tp2tt)) -
                           2*Tp1tt
          I33[[ii]][2,3] = I33[[ii]][2,3] + sum(diag(Tp2tt %*% Te2tt)) -
                           2*Tp13tt
          I33[[ii]][3,2] = I33[[ii]][2,3]
          I33[[ii]][3,3] = I33[[ii]][3,3] + sum(diag(Te2tt %*% Te2tt)) -
                           2*Te1tt
        }
      }
    } else {
      Path = factor(dat[[ii]]$Path)
      PathLevels = levels(Path)
      nplev = length(PathLevels)
      ipvc = vc[[ii]]
      te = tp1 = tp12 = 0
      if(mode[ii] == "REML"){
        U = matrix(0,nmu+1,nmu+1) 
        Te1 = Tp1 = Tp12 = matrix(0,nmu+1,nmu+1)
        Te2 = Tp2 = matrix(0,nmu+1,nmu+1)
      }
      for(jj in 1:nplev){
        tpdat = dat[[ii]][dat[[ii]]$Path == PathLevels[jj],]
        ntpdat = nrow(tpdat)
        X = cbind(rep(1,ntpdat), tpdat$Log10.Delta.)
        W = tpdat$Log10.Yield.
        Omega = matrix(ipvc[1],ntpdat,ntpdat)
        for(kk in 1:ntpdat){
          Omega[kk,kk] = Omega[kk,kk] + ipvc[nvc[ii]]
        }
        cOmega = chol(Omega)
        Xt = backsolve(cOmega,X,transpose=TRUE)
        yt = backsolve(cOmega,W,transpose=TRUE)
        yt2 = t(yt) %*% yt

        #
        # Information matrix fixed effects
        #
        I11[[ii]] = I11[[ii]] + t(Xt) %*% Xt
        I12[[ii]] = I12[[ii]] + t(Xt) %*% yt
        I22[[ii]] = I22[[ii]] + yt2

        dOmegaR = diag(ntpdat)
        dOmegaRt = backsolve(cOmega,dOmegaR,transpose=TRUE)
        dOmegaRtt = backsolve(cOmega,dOmegaRt)
        te = te + sum(diag(dOmegaRtt %*% dOmegaRtt))
        dOmegaP = matrix(1,ntpdat,ntpdat)
        dOmegaPt = backsolve(cOmega,dOmegaP,transpose=TRUE)
        dOmegaPtt = backsolve(cOmega,dOmegaPt)
        tp1 = tp1 + sum(diag(dOmegaPtt %*% dOmegaPtt))
        tp12 = tp12 + sum(diag(dOmegaPtt %*% dOmegaRtt))
        if(mode[ii] == "REML"){
          C = cbind(X,W)
          Ct = cbind(Xt,yt)
          Ctt = backsolve(cOmega,Ct)
          U = U + t(Ct) %*% Ct
          Te1 = Te1 + t(C) %*% dOmegaRtt %*% dOmegaRtt %*% Ctt
          Te2 = Te2 + t(C) %*% dOmegaRtt %*% Ctt
          Tp1 = Tp1 + t(C) %*% dOmegaPtt %*% dOmegaPtt %*% Ctt
          Tp12 = Tp12 + t(C) %*% dOmegaPtt %*% dOmegaRtt %*% Ctt
          Tp2 = Tp2 + t(C) %*% dOmegaPtt %*% Ctt
        }
      }
      if(mode[ii] == "REML"){
        cU = chol(U)
        Te1t = backsolve(cU,Te1,transpose=TRUE)
        Te1tt = sum(diag(backsolve(cU,Te1t)))
        Te2t = backsolve(cU,Te2,transpose=TRUE)
        Te2tt = backsolve(cU,Te2t)
        Tp1t = backsolve(cU,Tp1,transpose=TRUE)
        Tp1tt = sum(diag(backsolve(cU,Tp1t)))
        Tp12t = backsolve(cU,Tp12,transpose=TRUE)
        Tp12tt = sum(diag(backsolve(cU,Tp12t)))
        Tp2t = backsolve(cU,Tp2,transpose=TRUE)
        Tp2tt = backsolve(cU,Tp2t)
      }

      #
      # Information matrix variance components
      #
      I33[[ii]][1,1] = tp1
      I33[[ii]][1,2] = tp12
      I33[[ii]][2,1] = I33[[ii]][1,2]
      I33[[ii]][2,2] = te
      if(mode[ii] == "REML"){
        I33[[ii]][1,1] = I33[[ii]][1,1] + sum(diag(Tp2tt %*% Tp2tt)) -
                         2*Tp1tt
        I33[[ii]][1,2] = I33[[ii]][1,2] + sum(diag(Tp2tt %*% Te2tt)) -
                         2*Tp12tt
        I33[[ii]][2,1] = I33[[ii]][1,2]
        I33[[ii]][2,2] = I33[[ii]][2,2] + sum(diag(Te2tt %*% Te2tt)) -
                         2*Te1tt
      }
    }
    I33[[ii]] = I33[[ii]]/2
    cI11 = chol(I11[[ii]])
    at = backsolve(cI11,I12[[ii]],transpose=TRUE)
    I2g1[[ii]] = I22[[ii]] - t(at) %*% at
  }
  return(list(I11 = I11, I12 = I12, I22 = I22, I33 = I33, I2g1 = I2g1))
}
