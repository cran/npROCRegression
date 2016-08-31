!    ******************************************************************
!    ******************************************************************
     subroutine init_random_seed(seed)
     implicit none
     integer, parameter :: ssize = selected_int_kind (8)
     integer seed, i
     integer :: nseed = 12
     integer (kind = ssize) :: t
     integer (kind = ssize), allocatable :: iseed(:)
     
     call random_seed(size = nseed)
     allocate(iseed(nseed))
     do i = 1, nseed
          t = seed*i
          iseed(i) = t
     end do
     call random_seed(put = iseed)    
     deallocate(iseed)
     end subroutine init_random_seed
!    ******************************************************************
!    ******************************************************************
     subroutine sampleROC(Z0, X0, W0, n0, Z1, X1, W1, n1, nvar, &
     coutcome, Z0b, X0b, W0b, n0b, Z1b, X1b, W1b, n1b)
     implicit none
     integer n0, n1, n0b, n1b, nvar, ir0(n0), ir1(n1), ir(n0+n1), i, j       
     logical coutcome
     double precision X0(n0), X1(n1), Z0(n0,nvar), Z1(n1,nvar), &
          W0(n0), W1(n1), &
          X0b(n0+n1), X1b(n0+n1), &  
          Z0b(n0+n1,nvar), & 
          Z1b(n0+n1,nvar), &
          W0b(n0+n1), W1b(n0+n1)
     Z0b = 0.0
     X0b = 0.0
     W0b = 0.0
     Z1b = 0.0
     X1b = 0.0
     W1b = 0.0     
     if(coutcome) then
          n0b=n0
          n1b=n1
          call sample_int(n1,n1,ir1)
          do i=1,n1 
               X1b(i)=X1(ir1(i))
               W1b(i)=W1(ir1(i))
               do j=1,nvar
                    Z1b(i,j)=Z1(ir1(i),j)
               end do
          end do
          call sample_int(n0,n0,ir0)
          do i=1,n0 
               X0b(i)=X0(ir0(i))
               W0b(i)=W0(ir0(i))
               do j=1,nvar
                    Z0b(i,j)=Z0(ir0(i),j)
               end do
          end do
     else
          n0b=0
          n1b=0
          call sample_int(n0+n1,n0+n1,ir)
          do i=1,n0+n1
               if(ir(i).le.n0) then
                    n0b = n0b+1
                    X0b(n0b)=X0(ir(i)) 
                    W0b(n0b)=W0(ir(i)) 
                    do j=1,nvar
                         Z0b(n0b,j)=Z0(ir(i),j)
                    end do
               else
                    n1b = n1b+1
                    X1b(n1b)=X1(ir(i)-n0) 
                    W1b(n1b)=W1(ir(i)-n0) 
                    do j=1,nvar
                         Z1b(n1b,j)=Z1(ir(i)-n0,j)
                    end do
               end if
          end do
     end if
     end
!    ******************************************************************
!    ******************************************************************
     double precision function SD(x,Ed,W,n)
     implicit none
     integer i,n
     double precision x,Ed(n),W(n),sumw
     Sd=0
     SumW=0
     do i=1,n
          if (W(i).gt.0) then
               sumw=sumw+W(i)
               if (Ed(i).ge.x) Sd=Sd+W(i)
          end if
     end do
     if (sumw.gt.0) Sd=Sd/sumw
     end
!    ******************************************************************
!    ******************************************************************
     subroutine SH_(t,nt,EH,W,n,SH)
     implicit none
     integer nt,n,i,it
     double precision t(nt),EH(n),W(n),SH(nt)
     double precision SD2(n)
     double precision,external::SD, maximum, minimum
     
     do i=1,n
          Sd2(i)= SD(EH(i),EH,W,n)
     end do
     SH = maximum(EH,n)
     do it=1,nt
          do i=1,n
               if (Sd2(i).le.t(it).and.W(i).gt.0) then
                    if (Eh(i).le.SH(it)) SH(it)=Eh(i)
               end if
          end do
     end do
     end
!    ******************************************************************
!    ******************************************************************
     subroutine cROC(m0,m1,v0,v1,Err1,w1,n1,nt, &
          SH, ROC)
     implicit none
     integer n1, nt, it
     double precision m0, m1, v0, v1, w1(n1), &
          Err1(n1), auxM, auxV, auxV1, SH(nt), ROC(nt)
     double precision,external::SD
          auxM  = (m0-m1)/sqrt(v1)
          auxV1 = (sqrt(v0)/sqrt(v1))
          do it = 1, nt
               auxV =  SH(it)*auxV1
               ROC(it) = SD(auxM + auxV, Err1, w1, n1)
          end do
     end
!    ******************************************************************
!    ******************************************************************
     double precision function cAUC(ROC,t,nt)
     implicit none
     integer nt, it
     double precision ROC(nt), t(nt)
     cAUC=0.0
     do it=1,(nt-3)/2
          cAUC = cAUC + 2.0*ROC(2*it+1)
     end do
     do it =1,(nt-1)/2
          cAUC = cAUC + 4.0*ROC(2*it)
     end do
     cAUC = (((t(nt)-t(1))/(nt-1))/3.0)*(cAUC + ROC(1) + ROC(nt))
     end
!    ******************************************************************     
!    ******************************************************************
     subroutine Rfast(X,Y,n,W,h,p,Xb,Pb,kbin)
     implicit none
     integer, parameter :: nh = 21
     integer n,i,j,kbin,p,icont,ih
     double precision x(n),y(n),W(n),Xb(kbin),Yb(kbin),Wb(kbin), &
          Area(2), dis1, dis2, Pb(kbin), h, &
          ErrCV,erropt,sumw,med,h0
     double precision, allocatable:: X2(:), W2(:), B(:), pred(:), &
          haux(:)
     allocate (Pred(n),X2(kbin),W2(kbin),B(10),haux(nh))
     
     !******************************************************************
     ! Binning sample
     !******************************************************************
     Wb=0
     Yb=0
     do i=1,n     
          if (W(i).gt.0) then
               if (X(i).lt.Xb(1)) then
                    Wb(1)=wb(1)+W(i)
                    yb(1)=yb(1)+W(i)*Y(i)
               elseif (X(i).gt.Xb(kbin)) then
                    Wb(kbin)=wb(kbin)+W(i)
                    yb(kbin)=yb(kbin)+W(i)*Y(i)
               else
                    do j=1,kbin-1
                         if (Xb(j).le.X(i).and.X(i).le.Xb(j+1)) then
                              dis1=X(i)-Xb(j)
                              dis2=Xb(j+1)-X(i)
                              Area(1)=dis2/(dis1+dis2)
                              Area(2)=dis1/(dis1+dis2)
                              Wb(j)=Wb(j)+W(i)*Area(1)
                              Yb(j)=Yb(j)+Y(i)*W(i)*Area(1)
                              Wb(j+1)=Wb(j+1)+W(i)*Area(2)
                              Yb(j+1)=Yb(j+1)+Y(i)*W(i)*Area(2)
                         end if
                    end do
               end if
          end if
     end do
     do i=1,kbin
          if (Wb(i).gt.0) Yb(i)=Yb(i)/Wb(i)
     end do
     
     if (h.eq.-3) then ! no effect
          Pb = 0.0
          goto 1
     elseif (h.eq.-2) then ! constant effect
          Pb = 0
          sumw = 0
          med = 0
          do i = 1,n
               med = med + W(i)*Y(i)     
               sumw = sumw + W(i)
          end do     
          if (sumw.gt.0) pb = med/sumw
          goto 1
     elseif (h.eq.0) then ! linear effect
          icont = n
          call RegLinealPred(X,Y,W,icont,p, pred, Xb,pb,kbin)
          goto 1
     elseif (h.lt.0) then ! smooth effect by CV
          do ih = 1, nh
               haux(ih) = (ih-1)*1.0/(nh-1)
          end do          
          haux(1) = 0.05
          haux(nh) = 0.95
          icont = 0
          h0 = haux(1)
          erropt = 9e9
          do ih = 1, nh
               h = haux(ih)
               call Rfast_(h,p,Xb,Yb,Wb,Pb,kbin,1)
               ErrCV = 0
               do i = 1,kbin
                    ErrCV = ErrCV + Wb(i)*(Yb(i)-Pb(i))**2
               end do
               icont = icont + 1
               if (icont.eq.1) then
                    erropt = ErrCV
                    h0 = h
               else
                    if (ErrCV.lt.erropt) then
                         erropt = ErrCV
                         h0 = h
                    end if
               end if
          end do
          h = 1.0*h0
     end if
     call Rfast_(h,p,Xb,Yb,Wb,Pb,kbin,0)
1    continue
     deallocate(Pred,X2,W2, B, haux)
     end subroutine
!    ******************************************************************
!    ******************************************************************
     subroutine Rfast_(h,p,Xb,Yb,Wb,Pb, &
          kbin,ifcv)
     implicit none
     integer i, j, kbin, p, ifcv
     double precision, parameter :: deg2rad = 3.14159265/180.0
     double precision Xb(kbin),&
          Yb(kbin), Wb(kbin), Pb(kbin),&
          aux, h, sumw, h2, min, max, pdf
     double precision, allocatable:: X2(:), W2(:),&
          Kernel(:), B(:), V(:), pred(:), Xb2(:)
     allocate (Pred(kbin),X2(kbin),W2(kbin),&
          kernel(0:kbin),&
          B(10),V(kbin), Xb2(kbin))
     h2 = h
     W2 = 1
     call min_y_max(Xb,kbin,min,max,W2)
     do i = 1, kbin
          Xb2(i) = (Xb(i) - min)/(max - min)
     end do
     do i = 1,kbin
          sumw = 0
          W2 = 0
          V = 0
          X2 = 0
          do j=1,kbin
               aux = (Xb2(j)- Xb2(i))/h2
               pdf = exp(-0.5*(aux**2))/sqrt(2*3.141593)
               W2(j) = Wb(j)*pdf/h2
               V(j) = Yb(j)
               X2(j) = Xb(j)-Xb(i)
               sumw = sumw + W2(j)
          end do
          if (ifcv.gt.0) W2(i) = 0
          W2 = W2/sumw
          call Reglineal(X2,V,W2,kbin,p,B,pred)
          pb(i) = B(1)
     end do
     deallocate(Pred,X2,W2,kernel,B,V,Xb2)
     end
!    ******************************************************************
!    ******************************************************************
     double precision function QQ(X,n,alpha)
     implicit none
     integer n
     double precision alpha, newalpha(1),X(n), Q(1)
          newalpha(1) = alpha
          call quantile(X,n,newalpha,1,Q)
          QQ = Q(1)
     end
!    ******************************************************************
!    ******************************************************************
     subroutine quantile(X,n,alpha,nalpha,Q)
     implicit none
     integer n,nalpha,ip,j,ind(n)
     double precision X(n),alpha(nalpha),Q(nalpha),R
     call qsortd(x,ind,n)
     
          do j=1,nalpha
               IP = floor(alpha(j)*(n+1.))
               IF(IP .lt. 1) then
                    Q(j)=X(ind(1))
               elseif (IP.ge.n) then
                    Q(j)=X(ind(n))
               else
                    R=alpha(j)*(n+1.)-IP
                    Q(j)=(1.-R)*X(ind(IP)) + R*X(ind(IP+1))
               end if
          end do
     end
!    ******************************************************************
!    ******************************************************************
     double precision function Meanf(X,n)
     implicit none
     integer i,n
     double precision X(n)
          Meanf = 0.0
          do i=1,n
               Meanf = Meanf + X(i)
          end do
          Meanf = Meanf/(1.0*n)
     end
!    ******************************************************************
!    ******************************************************************
     subroutine sample_int(n,size,II)
     implicit none
     integer n,size,II(n),i
     double precision ru
     do i=1,size
          call random_number(ru)
          II(i)= int(ru*n + 1)
          if (II(i).le.1) II(i) = 1
          if (II(i).ge.n) II(i) = n
     end do
     end
!    ******************************************************************
!    ******************************************************************
     Module Data
     integer q
     integer,allocatable::nf(:)
     double precision,allocatable::mode(:,:),X(:,:,:), Xp(:,:,:), &
          Xpar(:,:), Xppar(:,:),B(:),Fact(:,:)
     end Module
!    ******************************************************************
!    ******************************************************************
     module Mod0
     integer minit, maxit
     double precision eps
     end Module
!    ******************************************************************
!    ******************************************************************
     subroutine FinGam()
     use Mod0
     use Data
     implicit none
     deallocate (X,Xp,Mode,Fact,nf,Xpar,Xppar,B)
     end
!    ***********************************************************************
!    ***********************************************************************
     subroutine IniGAM(n, npred,nvar,npar, Xdata, Xpred, Mode0, II, h,&
          vpar, n_vpar, vnpar, n_vnpar)
     use Data
     use Mod0
     implicit none
     integer n,nvar,npar,npred,r,s,i,j,l,Mode0(nvar-1), II(2,npar), & 
          vpar(npar), vnpar(npar), n_vpar, n_vnpar,k,aux
     double precision, allocatable:: Xaux(:,:), Xpaux(:,:)
     double precision Xdata(n,nvar-1), Xpred (npred, nvar-1),h(n,npar)
     allocate (X(n,2,npar), Xp(npred,2,npar), mode(2,npar), &
          Fact(n,npar),nf(npar))
     ! Type of effects
     do r=1,npar
          do s=1,2
               mode(s,r)=0
               if (ii(s,r).gt.0.and.ii(s,r).lt.nvar) & 
                    Mode(s,r) = mode0(ii(s,r))     
          end do
     end do     
     ! Data matrix
     do r=1,npar
          do s=1,2
               if (II(s,r).le.0) then
                    do i=1,n
                         X(i,s,r)=0
                    end do
                    do i=1,npred
                         Xp(i,s,r)=0
                    end do
               elseif (II(s,r).lt.nvar) then
                    do i=1,n
                         X(i,s,r)=Xdata(i,II(s,r))
                    end do
                    do i=1,npred
                         Xp(i,s,r)=Xpred(i,II(s,r))
                    end do
               end if
         end do
     end do     
     vpar=0
     vnpar=0
     n_vpar=0
     n_vnpar=0
     do j=1,npar
          if((Mode(1,j).eq.6.and.Mode(2,j).eq.0).or. &
               (Mode(1,j).eq.0.and.Mode(2,j).eq.6).or. &
               (((Mode(1,j).eq.5.and.Mode(2,j).eq.0).or. &
               (Mode(1,j).eq.0.and.Mode(2,j).eq.5)).and. &
                    h(1,j).eq.0)) then
               n_vpar=n_vpar+1
               vpar(n_vpar)=j
          else
               n_vnpar=n_vnpar+1
               vnpar(n_vnpar)=j
          end if
     end do
     q=0
     do i=1,n_vpar
          if(Mode(1,vpar(i)).eq.5.or.Mode(2,vpar(i)).eq.5) then ! Continuous
               q=q+1
               nf(vpar(i))=1
          else if(Mode(1,vpar(i)).eq.6.or.Mode(2,vpar(i)).eq.6) then ! Categorical
               aux=1
               if(Mode(1,vpar(i)).eq.0) aux=2
               call getLevels(X(1,aux,vpar(i)),n,Fact(1,vpar(i)),nf(vpar(i)))
               nf(vpar(i)) = nf(vpar(i))-1  ! Number of levels - 1
               q=q+nf(vpar(i))
          end if
     end do
     allocate (Xpar(n,q),Xppar(npred,q),B(q+1))
     k=0
     do j=1,n_vpar
          if(Mode(1,vpar(j)).eq.5.or.Mode(2,vpar(j)).eq.5) then
               aux=1
               if(Mode(1,vpar(j)).eq.0) aux=2
               k=k+1
               do i=1,n
                    Xpar(i,k)=X(i,aux,vpar(j))
               end do
               do i=1,npred
                    Xppar(i,k)=Xp(i,aux,vpar(j))
               end do
          else if(Mode(1,vpar(j)).eq.6.or.Mode(2,vpar(j)).eq.6) then
               aux=1
               if(Mode(1,vpar(j)).eq.0) aux=2
               allocate(Xaux(n,nf(vpar(j))),Xpaux(npred,nf(vpar(j))))
               call getModelMatrixFact(X(1,aux,vpar(j)),n,Xaux, &
                     Xp(1,aux,vpar(j)),npred,Xpaux,nf(vpar(j))+1)
               do i=1,n
                    do l=1,nf(vpar(j))
                         Xpar(i,k+l)=Xaux(i,l)
                    end do
               end do
               do i=1,npred
                    do l=1,nf(vpar(j))
                         Xppar(i,k+l)=Xpaux(i,l)
                    end do
               end do
               deallocate(Xaux,Xpaux)
               k=k+nf(vpar(j))
          end if
     end do
     end
!    ***********************************************************************
!    ***********************************************************************
     subroutine getLevels(X,n,Fact,nf)
     implicit none
     integer n,nf,i,j
     double precision X(n),Fact(n)
     logical ifnew
     nf = 0
     ifnew = .TRUE.
     do i = 1,n
          if (i.eq.1) then
                nf = nf + 1
                Fact(1) = X(i)
          else
               ifnew = .TRUE.
               do j = 1,nf
                    If(X(i).eq.Fact(j)) ifnew = .FALSE. 
               end do
               if (ifnew) then
                    nf = nf+1
                    Fact(nf) = X(i)
               end if
          end if
     end do          
     end subroutine
!    ***********************************************************************
!    ***********************************************************************
     subroutine getModelMatrixFact(X,n,X2,Xp,np,X2p,nf)
     implicit none
     integer n,nf,np,i,j
     double precision X(n),X2(n,nf-1),Xp(np),X2p(np,nf-1),Fact(n)
     call getLevels(X,n,Fact,nf)
     do i=1,n          
          do j=1,nf-1
               X2(i,j) = 0
               if (X(i).eq.Fact(j)) X2(i,j) = 1.0
               if (X(i).eq.Fact(nf)) X2(i,j) = -1.0
          end do
     end do
     do i=1,np          
          do j=1,nf-1
               X2p(i,j) = 0
               if (Xp(i).eq.Fact(j)) X2p(i,j) = 1.0
               if (Xp(i).eq.Fact(nf)) X2p(i,j) = -1.0
          end do
     end do
     end subroutine
!    ********************************************************************
!    ********************************************************************
!                              Local scoring algorithm
!    ********************************************************************
!    ********************************************************************
     subroutine GAM(n, nvar, npar, Mode0, II, Xdata, Y, W, h2, kbin, p,&
          family, F, coeff, muhat, Xpred, Fpred, muhatpred, npred)
     use data
     use mod0
     implicit none
     integer n,nvar,npar,II(2,npar),mode0(nvar-1),family, npred,&
          p(npar), kbin, nitl, maxitl
     double precision Xdata(n,nvar-1),Y(n),h2(n,npar),F(n,npar), &
          muhat(n),W(n), Xpred(npred, nvar-1),Fpred(npred,npar), &
          muhatpred(npred), coeff(20), Devian, thloc, linc
     double precision,allocatable::Eta(:),Etap(:)
     integer i, vpar(npar), vnpar(npar), n_vpar, n_vnpar
     double precision devnew,p0,var,Devold,aux,eta0,der,xmiss
     double precision,external::Slinc,Dev,diriv,weight
     double precision,allocatable::h0(:,:),Z(:),WZ(:),WX(:)
     allocate (h0(n,npar),Z(n),Wz(n),Wx(n),Eta(n), Etap(npred))
       
     minit=1
     maxit=10 
     thloc=.01
     Eps=.01
     maxitl=10
     xmiss=99999
     linc=family
     if(family.eq.2.or.family.eq.6) maxitl = 1  ! 6 censored data

     Wx=1
     do i=1,n
          if (Y(i).eq.xmiss) then
               w(i)=0
               Wx(i)=0
          end if
     end do
     
     call IniGAM(n,npred,nvar,npar,Xdata,Xpred,Mode0,II,h2,&
          vpar, n_vpar, vnpar, n_vnpar) ! Load the data
               
     call Mean_and_Var(Y,W,n,p0,var)
     Muhat=p0
     Eta0=Slinc(p0,linc)
     devnew=dev(n,Muhat,y,w,family)
     
     F = 0
     do i=1,n
         Eta(i)=Eta0
         Muhat(i)=p0
     end do
     Fpred = 0
     ! Loop
     do nitl = 1,maxitl
          do i=1,n
               der = diriv(muhat(i),linc)
               Z(i) = Eta(i)+(Y(i)-muhat(i))*der
               Wz(i) = weight(w(i),muhat(i),family,linc)
          end do
          h0=h2          
          call BackFitInter(n,npar,Z,Wz,Wx,h0,p,kbin,&  ! Input
               vpar,n_vpar,vnpar,n_vnpar,&              ! Input
               Eta,F,Etap,Fpred,npred)                  ! Output
          call Linv(n,Eta,Muhat,linc)
          devold = devnew
          devnew = dev(n,muhat,y,w,family)
          aux=abs((devold-devnew)/devold)
          if(aux.lt.thloc) then
               devian=devnew
               goto 23009
          end if
     end do
23009 continue
     h2=h0
     call Linv(npred,Etap,Muhatpred,linc)
     deallocate (h0,Z,Wz,Wx,Eta,Etap)
     
     do i=1,q+1
          coeff(i) = B(i)
     end do          
     call fingam()     
     end
!    ********************************************************************
!    ********************************************************************
!                                     Backfitting
!    ********************************************************************
!    ********************************************************************
     subroutine BackFitInter(n,npar,Y2,W2,Wx,h2,p,kbin,&
          vpar,n_vpar,vnpar,n_vnpar, &  
          Eta,F, Etap, Fp, np)
     Use Data
     Use Mod0
     implicit none
     integer i,j,n,npar, np, p(npar), n_vpar, n_vnpar, &
          vpar(n_vpar), vnpar(n_vnpar), kbin, nit
     double precision deltaf,ratio,normf,med,Eta(n),Y2(n),W2(n), &
          F(n,npar),aux, Wx(n), Etap(np), Fp(np,npar)
     double precision,allocatable::Z2(:),Old(:),etal(:),etalp(:)
     double precision,allocatable::h0(:,:),h3(:,:)
      
     double precision h2(n,npar)
     
     double precision,external::dnorm2
     allocate (Old(n),Z2(n),h0(n,npar),h3(n,npar),etal(n),&
          etalp(np))

     ! Parametric part
     do i=1,n    
          eta(i)=0
          do j=1,n_vnpar
               Eta(i)=Eta(i)+F(i,vnpar(j))
          end do
     end do
     ! Only one covariate or all parametric
     if (npar.le.1.or.n_vnpar.eq.0) then 
          minit=1
          maxit=1
     end if
     do nit=1,maxit
     !    Parametric part
          do i=1,n          
               Z2(i)=Y2(i)-eta(i)
          end do
          call Param(n,npar,Z2,W2,etal,F,etalp,Fp,np,vpar,n_vpar)
     !    Smooth part
          deltaf=0     
          do j=1,n_vnpar
               do i=1,n
                    old(i)=F(i,vnpar(j))
                    Z2(i)=Y2(i)-etal(i)-eta(i)+old(i)
               end do                    
               do i=1,n
                    h0(i,vnpar(j))=h2(i,vnpar(j))
               end do

               call Finter(X(1,1,vnpar(j)),Z2,n,W2,Wx,h0(1,vnpar(j)), &
               p(vnpar(j)),Fact(1,vnpar(j)),nf(vnpar(j)), &
               mode(1,vnpar(j)),kbin, F(1,vnpar(j)), &
               Xp(1,1,vnpar(j)), Fp(1,vnpar(j)), np)
               
               do i=1,n
                    h3(i,vnpar(j))=h0(i,vnpar(j))
               end do
     !         *************************
     !               Center the effect
     !         *************************
     !         Data               
               call Mean_and_Var(F(1,vnpar(j)),Wx,n,Med,aux)
               do i=1,n
                    F(i,vnpar(j))=F(i,vnpar(j))-Med
                    eta(i)=eta(i)+F(i,vnpar(j))-old(i)
               end do          
     !         Prediction
               do i=1,np
                    Fp(i,vnpar(j))=Fp(i,vnpar(j))-Med
               end do
               deltaf=deltaf+dnorm2(n,old,F(1,vnpar(j)),W2)          
          end do
          normf = 0e0
          ratio = 1e9   
          do i=1,n
               normf = normf + W2(i)*eta(i)**2
           end do
          if(normf.gt.0) ratio = sqrt(deltaf/normf)
          if (ratio.le.eps.and.nit.ge.minit) goto 1
     end do
1     continue
     h2 = h3
     !Add the parametric part
     do i=1, n          
          eta(i)= eta(i) + etal(i)
     end do
     !Additive predictor
     etap=0
     do i=1,np          
          etap(i) = etalp(i)
          do j=1,n_vnpar
               etap(i)= etap(i) + Fp(i,vnpar(j))
          end do
     end do
     deallocate (Old,Z2,h0,h3,etal,etalp)
     end
!    ********************************************************************
!    ********************************************************************
     subroutine Finter(X,Y2,n,W,Wx,h2,p,Fact,nf,mode,kbin,F,Xp,Fp, np)
     implicit none
     integer n, nf, np, p, kbin
     double precision X(n,2),Y2(n),W(n),Wx(n),mode(2),F(n), &
          Fact(n), Xp(np,2), Fp(np), h2(n)

     if (mode(1).eq.0.and.mode(2).eq.5) then        ! %%%%%%   ooo -- Continuous
          call RNP1DFast(X(1,2),Y2,n,W,WX,h2(2),p,kbin,F, &
          Xp(1,2), Fp, np)
     elseif (mode(1).eq.5.and.mode(2).eq.0) then    ! %%%%%    Continuous -- ooo
          call RNP1DFast(X(1,1),Y2,n,W,WX,h2(1),p,kbin,F, &
          Xp(1,2), Fp, np)
     elseif(mode(1).eq.0.and.mode(2).eq.6) then     ! %%%%%%   ooo - Factor
          call RegFact(X(1,2),Y2,W,n,F, Xp(1,2), Fp, np)
     elseif(mode(1).eq.6.and.mode(2).eq.0) then     ! %%%%%    Factor-ooo
          call RegFact(X(1,1),Y2,W,n,F, Xp (1,2), Fp, np)
     elseif (mode(1).eq.6.and.mode(2).eq.5) then    ! %%%%%    Factor-Continuous
          nf=0
          if (nf.le.0) call getLevels(X(1,1),n,Fact,nf)
          call FRNP1DFast(X(1,1),X(1,2),Y2,W,Wx,Fact,n,nf,h2,p,kbin, &
               F, Xp(1,1), Xp(1,2), Fp, np)
     elseif (mode(1).eq.5.and.mode(2).eq.6) then    ! %%%%%    Continuous-Factor
          nf=0
          if (nf.le.0) call getLevels(X(1,2),n,Fact,nf)
          call FRNP1DFast(X(1,2),X(1,1),Y2,W,Wx,Fact,n,nf,h2,p,kbin, &
               F, Xp(1,2), Xp(1,1), Fp, np)
     end if
     end
!    ***********************************************************************
!    ***********************************************************************
     subroutine linv(n,etahat,muhat,linc)
     implicit none
     integer n
     double precision linc,etahat(n),muhat(n)
     if(linc.eq.2) then !IDENTITY
          call linvid(n,etahat,muhat)
     elseif(linc.eq.1) then !LOGIT
          call linvlt(n,etahat,muhat)
     elseif(linc.eq.5) then !LOGAR
          call linvlo(n,etahat,muhat)
     elseif(linc.eq.4) then !INVER
          call linvin(n,etahat,muhat)
     elseif(linc.eq.7) then !PROBIT
          call linvpr(n,etahat,muhat)
     elseif(linc.eq.8) then !CLOGLOG
          call linvcll(n,etahat,muhat)
     end if
     end
!     *********************************************************
!     *********************************************************
     subroutine linvid(n,etahat,muhat)
     implicit none
     integer i,n
     double precision muhat(n),etahat(n)
     do i = 1,n
         muhat(i)=etahat(i)
     end do
     end
!     *********************************************************
!     *********************************************************
     subroutine linvlt(n,etahat,muhat)
     implicit none
     integer n,i
     double precision muhat(n),etahat(n)
     double precision pr
     do i = 1,n 
          if (etahat(i).gt.10) then
               pr=exp(10.0)
          elseif (etahat(i).lt.-10) then
               pr=exp(10.0)
          else
               pr=exp(etahat(i))
          end if
          pr=pr/(1+pr)
          muhat(i)=pr
     end do
     end
!     *********************************************************
!     *********************************************************
     subroutine linvcll(n,etahat,muhat)
     implicit none
     integer n,i
     double precision muhat(n),etahat(n), vmax, vmin
     vmax = 0.9999
     vmin = 0.0001
     do i=1,n
          muhat(i) = 1.0 - exp(-exp(etahat(i)))
          muhat(i) = min(muhat(i), vmax)
          muhat(i) = max(muhat(i), vmin)
     end do
     end
!     *********************************************************
!     *********************************************************
     subroutine linvpr(n,etahat,muhat)
     implicit none
     integer n,i
     double precision muhat(n),etahat(n)
     double precision, external :: normal   
     do i=1,n                     
          muhat(i)=normal(etahat(i))
     end do
     end
!     *********************************************************
!     *********************************************************
     subroutine linvlo(n,etahat,muhat)
     implicit none
     integer i,n
     double precision muhat(n),etahat(n)
     do i = 1,n 
          if (etahat(i).le.88) then
               muhat(i)=exp(etahat(i))
          else
               muhat(i)=exp(88.0)
          end if
     end do
     end
!     *********************************************************
!     *********************************************************
     subroutine linvin(n,etahat,muhat)
     implicit none
     integer n,i
     double precision muhat(n),etahat(n)
     do 23082 i=1,n 
          if(.not.(etahat(i).lt.0.0001))goto 23084
          muhat(i)=1.0/0.0001
          goto 23085
23084          continue
          muhat(i)=1/etahat(i)
23085          continue
23082 continue
     return
     end
!    ******************************************************************
!    ******************************************************************
     double precision function SLINC(muhat,linc)
     implicit none
     double precision muhat,linc
     double precision,external::lincid,linclt,linclo,lincin,lincpr,linccll

     if(linc.eq.2) then     !Identity
          SLINC=LINCID(muhat)
     elseif(linc.eq.1) then !Logit
          SLINC=LINCLT(muhat)
     elseif(linc.eq.5) then !Logar
          SLINC=LINCLO(muhat)
     elseif(linc.eq.4) then !Inver
          SLINC=LINCIN(muhat)
     elseif(linc.eq.7) then !Probit
          SLINC=LINCPR(muhat)
     elseif(linc.eq.8) then !cloglog
          SLINC=LINCCLL(muhat)
     else
          SLINC=LINCID(muhat)
     end if
     end
!     *********************************************************
!     ********************************************************* 
     double precision function LINCID(muhat)     !Link identity
     implicit none
     double precision muhat
     LINCID=muhat
     end
!     *********************************************************
!     *********************************************************
     double precision function LINCLT(muhat)     !Link Logit
     implicit none
     double precision muhat, logit, d
     logit=muhat
     d=1.0-logit
     if(d.lt.0.0001) d=0.0001
     logit=logit/d
     if(logit.lt.0.0001) then
          LINCLT=log(0.0001)
     elseif(logit.gt.9999.0) then
          LINCLT=log(9999.0)
     else
          LINCLT=log(muhat)
     end if
     end
!     *********************************************************
!     *********************************************************
     double precision function LINCPR(muhat)       !Link Probit
     implicit none
     double precision muhat
     double precision, external :: normdev
     LINCPR = normdev(muhat)      
     end
!     *********************************************************
!     *********************************************************
     double precision function LINCCLL(muhat)       !Link cloglog
     implicit none
     double precision muhat
     LINCCLL = log(-log(1.0 - muhat))      
     end     
!     *********************************************************
!     *********************************************************
     double precision function LINCLO(muhat)     !Link Logarithm
     implicit none
     double precision muhat
     if(muhat.le.0.0001) then
          LINCLO = log(.0001)
     else
          LINCLO = log(muhat)
     end if
     end
!     *********************************************************
!     *********************************************************
     double precision function LINCIN(muhat) !Link inverse
     implicit none
     double precision muhat
     LINCIN=1.0/muhat
     end
!     *********************************************************
!     *********************************************************
!     *********************************************************
!     *********************************************************
!                    Local linear smoothers
!
!                    RNP1DFAST  ---> Univariate case.
!                    FRNP1DFast ---> With factors.
!     *********************************************************
!     *********************************************************
!               RNP1DFAST   (Univariate local linear smoother) 
!     *********************************************************
!     *********************************************************
     subroutine RNP1DFast(X,Y,n,WY,WX,h,p,kbin,M0, &
          Xp,Fp,np)
     implicit none
     integer n, kbin, icont, i, p, np
     double precision X(n),Y(n),WY(n),WX(n),M0(n),&
          xmin,xmax, Xp(np), Fp(np), h
     double precision,allocatable::Xb(:),M0grid(:),yspl(:)
     double precision,External:: DEV
     
     allocate (Xb(kbin), M0grid(kbin),yspl(kbin))

     xmin=9e9
     xmax=-xmin
     do i=1,n
          xmin=min(xmin, X(i))
          xmax=max(xmax, X(i)) 
     end do
     do i=1,kbin
          Xb(i)=xmin+(i-1)*(xmax-xmin)/(kbin-1)
     end do
     Wx = 1.0
     if (h.eq.0) then
          icont=n
          call RegLinealPred(X,Y,WY,icont,p,M0,Xp,Fp,np)          
     else
          call Rfast(X,Y,n,Wy,h,p,Xb,M0grid,kbin)
          call spline(Xb, M0grid, kbin, yspl) 
          do i=1,n
               call splint(Xb, M0grid, yspl, kbin, x(i), m0(i))
          end do
          do i=1,np
          call splint(Xb, M0grid, yspl, kbin, xp(i), Fp(i))
          end do
     end if
     deallocate (Xb, M0grid, yspl)
     end
!    ******************************************************************
!    ******************************************************************
     subroutine FRNP1DFast(XF,X,Y,W,Wx,Fact,n,nf,h,p,kbin,F, &
          XFp, Xp, Fp, np)
     implicit none
     integer nf,n,kbin,i,j, np,p
     double precision X(n),XF(n),Y(n),W(n),Wx(n),Fact(nf),Med,F(n), &  
          Xp(np), Fp(np), XFp(np)
     double precision,allocatable::W2(:),Res(:),GrF(:,:), &
          GrFp(:,:)
     double precision h(nf)
     allocate(W2(n),GrF(n,nf),Res(n), GrFp(np, nf))

     do j=1,nf            
          W2=0
          do i=1,n
               if (XF(i).eq.Fact(j)) W2(i)=W(i)
          end do
          if (h(j).eq.-2) then
               do i=1,n
                    GrF(i,j)=0
               end do
          elseif (h(j).eq.0) then
               call RegLinealPred(X,Y,W2,n,p,GrF(1,j),Xp,GrFp(1,j),np)
          else
               call RNP1DFast(X,Y,n,W2,WX,h(j),p,kbin,GrF(1,j), &
                    Xp,GrFp(1,j), np)
          end if
     end do
     !****************************************************
     ! Center the effects
     !****************************************************
     ! Data
     do i=1,n              
          Med=0
          do j=1,nf
               Med=Med+GrF(i,j)/nf
          end do
          do j=1,nf
               GrF(i,j)=GrF(i,j)-Med
          end do
     end do
     ! Prediction
     do i=1,np
          Med=0
          do j=1,nf
               Med=Med+GrFp(i,j)/nf
          end do
          do j=1,nf
               GrFp(i,j)=GrFp(i,j)-Med
          end do
     end do
     !******************************************************
     ! Center the partial effects
     !******************************************************
     do i=1,nf
          do j=1,n
               W2(j)=0
               if (XF(j).eq.Fact(i)) W2(j)=Wx(j)
          end do
          call Mean(GrF(1,i),W2,n,Med)
          do j=1,n
               GrF(j,i) = GrF(j,i) - Med 
          end do
          do j=1,np
               GrFp(j,i) = GrFp(j,i) - Med 
          end do
     end do
     !*************************************************************
     ! All in one vector
     !*************************************************************
     ! Data
     do i=1,n
          do j=1,nf
               if (XF(i).eq.Fact(j)) F(i) = GrF(i,j)
          end do
     end do     
     ! Prediction          
     do i=1,np
          do j=1,nf
               if (XFp(i).eq.Fact(j)) Fp(i) = GrFp(i,j)
          end do
     end do
      deallocate(W2,GrF,Res, GrFp)
     end
!     *********************************************************
!     *********************************************************
     subroutine Mean(vector,w,n,med)
     implicit none
     integer n,i 
     double precision vector(n),med,w(n),wtot
     med=0
     wtot=0
     do i=1,n
          med=med+w(i)*vector(i)
          wtot=wtot+w(i)
     end do
     if (wtot.gt.0) then
          med=med/wtot
     end if
     end
!    ********************************************************************
!    ********************************************************************
     subroutine RegFact(F,Y,W,n,M0, Fp, Mp, np)
     implicit none
     integer nf,n,i,j,np
     double precision F(n),Y(n),W(n),m0(n),Fact(n), Fp(np), Mp(np)
     double precision,allocatable::X2(:,:),beta(:), X2p(:,:)
     
     call getLevels(F,n,Fact,nf)

     allocate (X2(n,nf-1),Beta(nf), X2p(np,nf-1))
     ! Sample
     do i=1,n          
          do j=1,nf-1
               X2(i,j)=0
               if (F(i).eq.Fact(j)) X2(i,j)=1.0
          end do
     end do
     ! Prediction
     do i=1,np          
          do j=1,nf-1
               X2p(i,j)=0
               if (Fp(i).eq.Fact(j)) X2p(i,j)=1.0
          end do
     end do
     call Regl(X2,Y,W,n,nf-1,beta,M0)
     do i=1,np
          Mp(i) = beta(1)          ! Intercept
          do j=1,nf-1
               Mp(i) = Mp(i) + X2p(i,j)*beta(j)
          end do
     end do
     deallocate(X2,Beta, X2p)
      end
!    ************************************************************     
!    ************************************************************
     subroutine Regl(X,Y,W,n,p,Beta,Pred)
     implicit none
     integer n, p, iopt
     double precision X(n,p),Y(n),W(n),Pred(n),beta(p+1),&
          sterr(p+1),se,r2
     iopt=0
          call WRegresion(X,Y,W,n,p,beta,sterr,se,r2,iopt)
          call PredLineal (X,n,p,Beta,Pred)
     end
!    ************************************************************
!    ************************************************************
     subroutine PredLineal (X,n,p,B,Pred)
     implicit none
     integer i,n,j,p
     double precision X(n,p),B(p+1),Pred(n)
          
     Pred=0
     do i=1,n
          Pred(i)=B(1)
          do j=1,p
               Pred(i)=Pred(i)+B(j+1)*X(i,j)          
          end do
     end do
     end
!    ************************************************************
!    ************************************************************
     subroutine RegL1d(X,Y,n,W,m0,m1)
     implicit none
     integer n,i
     double precision X(n),Y(n),W(n),m0(n),m1(n),Beta(2)
     call RegL(X,Y,W,n,1,Beta,m0)
     do i=1,n
          m1(i)=Beta(2)
     end do
     end
!    ******************************************************************
!    ******************************************************************     
     subroutine Reglineal(X,Y,W,n,p,B,pred)
     implicit none
     integer i,n,j,p
     double precision X(n),Y(n),W(n),Pred(n),B(p+1)
     double precision, allocatable:: X2(:,:)
     allocate (X2(n,p))
     do i=1,n
          do j=1,p
               X2(i,j)=X(i)**j
          end do
     end do
     call Regl(X2,Y,W,n,p,B,Pred)
     deallocate (X2)
     end
!    ******************************************************************
!    ******************************************************************
     subroutine RegLinealPred(X,Y,W,n,p,F,Xp,Yp,np)
     implicit none
     integer n,np,i,p,j
     double precision X(n), Y(n), W(n), Xp(np), Yp(np), F(n)
     double precision, allocatable:: Xp2(:,:), B(:), pred(:)
     allocate (B(p+1),Xp2(np,p),pred(n))

     call Reglineal (X,Y,W,n,p,B,F)

!     Prediccion
     do i=1,np
          Yp(i)=B(1)   
          do j=1,p
               Xp2(i,j)=Xp(i)**j
               Yp(i)=Yp(i)+B(j+1)*Xp2(i,j)
          end do
     end do
     deallocate(B,Xp2,pred)
     end subroutine
!    ************************************************************
!    ************************************************************
     subroutine min_y_max(x,n,xmin,xmax,W)
     implicit none
     integer i,n
     double precision x(n),W(n),xmin,xmax
     do i=1,n
          if (W(i).gt.0) then
               xmin=x(1)
               xmax=x(1)
               goto 1
          end if
     end do
1     continue
     do i=1,n
          if (W(i).gt.0) then
               xmin=min(x(i),xmin) 
               xmax=max(x(i),xmax)
          end if
     end do 
     end
!    ************************************************************
!    ************************************************************
     subroutine Mean_and_Var(vector,w,n,Mean,Variance)
     implicit none
     integer n,i 
     double precision vector(n),Mean,Variance,w(n),wtot
     Mean=0
     Variance=0
     wtot=0
     do i=1,n
          Mean = Mean+w(i)*vector(i)
          wtot = wtot+w(i)
     end do
     Mean=Mean/wtot
     do i=1,n
          Variance = Variance+w(i)*((vector(i)-Mean)**2)/wtot
     end do
     end     
!    ******************************************************************
!    ******************************************************************
     double precision function minimum(x,n)
     implicit none
     integer n,i
     double precision x(n),temp
     temp=x(1)
     do i=1,n
          if(x(i).le.temp) temp=x(i)
     end do
     minimum=temp
     end
!    ******************************************************************
!    ******************************************************************
     double precision function maximum(x,n)
     implicit none
     integer n,i
     double precision x(n),temp
     temp=x(1)
     do i=1,n
          if(x(i).ge.temp) temp=x(i)
     end do
     maximum=temp
     end
!    ******************************************************************
!    ******************************************************************
     double precision function Reg_0(X,Y,n)
     implicit none
     integer n,i
     double precision X(n),Y(n),sumx2,sumxy
     sumx2=0
     sumxy=0
     do i=1,n
          sumx2=sumx2+X(i)**2
          sumxy=sumxy+X(i)*Y(i)
     end do
          Reg_0 = sumxy/sumx2
     end
!    ******************************************************************
!    ******************************************************************
!                    Deviance functions
!    ********************************************************************
!    ********************************************************************
     double precision function DEV(n,fits,y,w,family)
     implicit none 
     integer n,family
     double precision fits(n),y(n),w(n)
     double precision,external::DEVGAM,DEVPOI,DEVB,DEVG
     if(family.eq.2) then                                     !Gaussian
          dev = devg(n,fits,y,w)
     elseif (family.eq.1.or.family.eq.7.or.family.eq.8) then  !Binomial
          dev = devb(n,fits,y,w)
     elseif(family.eq.4) then                                 !Gamma
          dev = devgam(n,fits,y,w)
     elseif(family.eq.5) then                                 !Poisson
          dev = devpoi(n,fits,y,w)
     else
          dev = devg(n,fits,y,w)     
     end if
     end
!    ************************************************************
!    ************************************************************
     double precision function DEVG(n,fits,y,w)
     implicit none
     integer i,n
     double precision fits(n),y(n),w(n),rss
     rss=0.0
     do i=1,n
          rss=rss+w(i)*(y(i)-fits(i))*(y(i)-fits(i))
     end do
     devg=rss
     end
!    ************************************************************
!    ************************************************************
     double precision function DEVB(n,fits,y,w)
     implicit none
     integer i,n
     double precision fits(n),y(n),w(n)
     double precision pr,entrop,entadd,dev
     dev=0.0
     do i=1,n 
          pr=fits(i)
          if(pr.lt.0.01) pr=0.01
          if(pr.gt.0.99) pr=0.99
          if((1.0-y(i))*y(i).le.0.0) then
               entrop=0.0
          else
               entrop=2.0*w(i)*(y(i)*log(y(i)) + &
               (1.0-y(i))*log(1.0-y(i)))
          end if
          entadd=2.0*w(i)*(y(i)*log(pr)+(1-y(i))*log(1.0-pr))
          dev=dev+entrop-entadd
     end do
     devb = dev
     end
!    ************************************************************
!    ************************************************************
     double precision function DEVGAM(n,fits,y,w)
     implicit none
     integer i,n
     double precision fits(n),y(n),w(n),tempy,tempf
     devgam=0.0
     do i=1,n
          tempy=y(i)
          if(tempy.lt.0.0001) tempy=.0001
          tempf=fits(i)
          if(tempf.lt.0.0001) tempf=.0001
          devgam=devgam+2*w(i)*(-log(tempy/tempf)+(y(i)-fits(i))/tempf)
     end do
     return
     end
!    ************************************************************
!    ************************************************************  
     double precision function DEVPOI(n,fits,y,w)
     implicit none
     integer i,n
     double precision fits(n),y(n),w(n),tempf
     devpoi=0.0
     do i=1,n
          tempf=fits(i)
          if(tempf.lt.0.0001) tempf=.0001
          devpoi=devpoi+2*w(i)*(-y(i)*log(tempf)- (y(i)-fits(i)))
          if(y(i).gt.0.0) devpoi=devpoi+2*w(i)*y(i)*log(y(i))
     end do
     return
     end
!    ************************************************************
!    ************************************************************
     double precision function weight(w,muhat,family,linc)
     implicit none
     integer family
     double precision w,muhat,linc,temp1,temp,aux,vmax, vmin
     double precision,external::diriv
     
     temp = DIRIV(muhat,linc)
     vmax = 0.9999
     vmin = 0.0001
     
     if(family.eq.2) then !Gaussian
          weight = w/(temp*temp)
     elseif(family.eq.1.or.family.eq.7.or.family.eq.8) then   ! Binomial 
          if (temp.eq.0) then                                 ! family=1: logit
               weight = 0                                     ! family=7: probit  
               goto 1                                         ! family=8: cloglog  
          end if
          aux = min(vmax, muhat)
          aux = max(vmin, aux)
     
          temp1 = aux*(1.0 - aux)
          aux = temp1*temp**2
          weight = w/aux
1    continue
     elseif (family.eq.4) then !Gamma
          weight = w/(temp*temp*muhat*muhat)
     elseif (family.eq.5) then  !Poisson
          temp1 = muhat
          if(temp1.lt.0.0001) then
               temp1 = 0.0001
          end if
          weight = w/(temp*temp*temp1)
     else
          weight = w/(temp*temp)
     end if
     return
     end
!    ************************************************************
!    ************************************************************
     double precision function DIRIV(muhat,linc)
     implicit none
     double precision muhat,linc
     double precision,external::dirvlo,dirvlt,dirvin,dirvid,dirvpr,dirvcll
     if(linc.eq.2) then          !Identity
          DIRIV = 1
     elseif(linc.eq.1) then      !Logit
          DIRIV=dirvlt(muhat)
     elseif(linc.eq.5) then      !Logarithm
          DIRIV = dirvlo(muhat)
     elseif(linc.eq.4) then      !Inverse
          DIRIV = dirvin(muhat)
     elseif(linc.eq.7) then      !Probit
          DIRIV = dirvpr(muhat)
     elseif(linc.eq.8) then      !cloglog
          DIRIV = dirvcll(muhat)
     else
          DIRIV=1  
     end if
     return
     end
!    ************************************************************
!    ************************************************************
!                              Derivatives
!    ************************************************************
!    ************************************************************
     double precision function DIRVLT(muhat) !Logit
     implicit none
     double precision muhat, pr, vmax, vmin
          vmax = 0.9999
          vmin = 0.0001       
          pr = min(vmax, muhat)
          pr = max(vmin, muhat)
          pr=pr*(1.0-pr)
          dirvlt=1.0/pr
     end
!    ************************************************************
!    ************************************************************
     double precision function DIRVPR(muhat) !Probit
     implicit none
     double precision muhat, pr, vmax, vmin
     double precision, external :: normdev
          vmax = 0.9999
          vmin = 0.0001
          pr = min(vmax, muhat)
          pr = max(vmin, pr)
          pr = normdev(pr)
          pr = exp(-0.5*(pr**2))/sqrt(2*3.141593)
          dirvpr = 1.0/pr
     end
!    ************************************************************
!    ************************************************************
     double precision function DIRVCLL(muhat) !cloglog
     implicit none
     double precision muhat,etahat,aux, vmax, vmin
     double precision, external :: LINCCLL
          vmax = 700.0
          vmin = 0.0001
          etahat = LINCCLL(muhat)
          etahat = min(etahat, vmax)
          aux = max(exp(etahat)*exp(-exp(etahat)), vmin)
          dirvcll = 1.0/aux
     end
!    ************************************************************
!    ************************************************************       
     double precision function DIRVLO(muhat) !Logarithm
     implicit none
     double precision muhat
          if(muhat.le.0.0001) then
               DIRVLO=1.0/0.0001
          else
               DIRVLO=1.0/muhat
          end if
     return
     end
!    ************************************************************
!    ************************************************************       
     double precision function DIRVIN(muhat) !Inverse
     implicit none
     double precision muhat
          DIRVIN=-1.0/(muhat*muhat)
          return
     end

!    ******************************************************************
!    ******************************************************************
     double precision function dnorm2(n,Y1,Y2,W)
     implicit none
     integer n,i
     double precision Y1(n),Y2(n),w(n),wtot,wsum
     wsum=0e0
     wtot=0e0
     do i=1,n
          wsum=wsum+w(i)*(Y2(i)-Y1(i))**2
          wtot=wtot+w(i)
     end do
     dnorm2=0
     if(wtot.gt.0) dnorm2=wsum/wtot
     end    
!    ******************************************************************
!    ******************************************************************
     subroutine sampleBinning(X, n, W, Xb, nb, wb)
     implicit none
     integer n, nb, i, j
     double precision X(n), W(n) , Xb(nb), Wb(nb), Area(2), dis1, dis2

     !Binning sample
     Wb=0     
     do i=1,n     
          if (W(i).gt.0) then
               if (X(i).lt.Xb(1)) then
                    Wb(1)=wb(1)+W(i)
               elseif (X(i).gt.Xb(nb)) then
                    Wb(nb)=wb(nb)+W(i)
               else
                    do j=1,nb-1
                         if (Xb(j).le.X(i).and.X(i).le.Xb(j+1)) then
                              dis1=X(i)-Xb(j)
                              dis2=Xb(j+1)-X(i)
                              Area(1)=dis2/(dis1+dis2)
                              Area(2)=dis1/(dis1+dis2)
                              Wb(j)=Wb(j)+W(i)*Area(1)
                              Wb(j+1)=Wb(j+1)+W(i)*Area(2)
                         end if
                    end do
               end if
          end if
     end do
     end
!    ******************************************************************
!    ******************************************************************
     subroutine LocScaleGam(X,Y,W,n,nvar, mode0, &
           nparm,IIm, nparv, IIv, hm, hv, pm, pv, kbin, &
          M,V,Xp,Mp,Vp,np)
     implicit none
     integer n, nvar, mode0(nvar), &
          nparm, nparv, IIm(2,nparm), IIv (2,nparv), &
          pm(nparm), pv(nparv), np, naux, &
          kbin,i,j
     double precision X(n,nvar), Y(n), W(n), &
          hm(nparm), hv(nparv), &
          M(n), V(n), Xp(np,nvar), Mp(np), Vp(np), &
          Err(n), LnErr(n), sumv, theta, vmin
     double precision,external::Reg_0
     double precision,allocatable::LnV(:), LnVp(:), hmcall(:,:), &
          hvcall(:,:), Fm(:,:), Fmp(:, :), Fv(:,:), Fvp(:, :), eLnV(:), &
          eLnVp(:), coeff(:)
     allocate (LnV(n), LnVp(np), hmcall(n,nparm), hvcall(n,nparv), &
          Fm(n,nparm), Fmp(np, nparm), Fv(n,nparv), Fvp(np, nparv), &
          eLnV(n), eLnVp(np), coeff(20))

     naux = n
     vmin = 0.000001
     !******************************
     !     Mean
     !******************************
     sumv=0
     do i=1,n
          do j=1,nparm
               hmcall(i,j)=hm(j)
          end do
          do j=1, nparv
               hvcall(i,j)=hv(j)
               sumv= sumv + hv(j)
          end do     
     end do
     call GAM(naux, nvar+1, nparm, Mode0, IIm, X, Y, W, &
          hmcall, kbin, pm, 2, Fm, coeff, M, Xp, Fmp, Mp, np)
     !******************************
     !     Variance
     !******************************
     do i =1,naux
          Err(i)=(Y(i)-M(i))**2
          LnErr(i)=log(max(Err(i),vmin))
     end do
     if(nparv.eq.1.or.sumv.eq.0) then
          pv = 0
          call GAM(naux, nvar+1, nparv, Mode0, IIv, X, Err, W, &
            hvcall, kbin, pv, 2, Fv, coeff, V, Xp, Fvp, Vp, np)
     else
          call GAM(naux, nvar+1, nparv, Mode0, IIv, X, LnErr, W, &
               hvcall, kbin, pv, 2, Fv, coeff, LnV, Xp, Fvp, LnVp, np)
           do i=1,naux
               eLnV(i) = exp(LnV(i))
          end do
          do i=1,np
               eLnVp(i) = exp(LnVp(i))
          end do
          theta = Reg_0(eLnv,Err,naux)
          do i=1,naux
               V(i) =  theta*eLnV(i)
          end do
          do i=1,np
               Vp(i) = theta*eLnVp(i)
          end do
     end if
     do i=1,naux
          V(i) = max(V(i), vmin)
     end do
     do i=1,np
          Vp(i) = max(Vp(i), vmin)
     end do
     deallocate (LnV, LnVp, hmcall, hvcall, Fm, Fmp, Fv, &
          Fvp, eLnV, eLnVp, coeff)
     end subroutine
!    ******************************************************************
!    ******************************************************************
     subroutine Param(n,p,Y,W,M,F,Mp,Fp,np,vpar,nvpar)
     Use Data
     implicit none
     integer i, j, l, n, np, p, k, nvpar, vpar(nvpar)
     double precision Y(n), W(n), M(n), &
          F(n,p),Mp(np),Fp(np,p),Med, w2(n)
     do i=1,n
          w2(i)=sqrt(w(i))
     end do
     ! Only the intercept
     if(nvpar.eq.0) then
          call Mean(Y,W2,n,Med)
          M=Med
          Mp=Med
          goto 1
     end if
     !Weighted linear model fitting
     call Regl(Xpar,Y,W2,n,q,B,M)
     !Prediction data
     do i=1,np
          Mp(i)=B(1)
          do j=1,q
               Mp(i) = Mp(i)+ Xppar(i,j)*B(j+1)
          end do
     end do
     !Partial functions
     k=0
     do j=1,nvpar
          do i=1,n
               F(i,vpar(j)) = 0
               do l=1,nf(vpar(j))
                    F(i,vpar(j)) = F(i,vpar(j)) + B(k+l+1)*Xpar(i,k+l)
               end do
          end do
          do i=1,np
               Fp(i,vpar(j)) = 0
               do l=1,nf(vpar(j))
                    Fp(i,vpar(j)) = Fp(i,vpar(j)) + B(k+l+1)*Xppar(i,k+l)
               end do
          end do
          k=k+nf(vpar(j))
     end do
1     continue
     end subroutine
!    ******************************************************************
!    ******************************************************************
     double precision function generateRV(t,ROC,nt)
     implicit none
     integer nt, i
     double precision t(nt), ROC(nt), ru, min
     call random_number(ru)
     min = 1.0
     do i=1,nt
          if(ROC(i).ge.ru) then
               min=t(i)
               goto 2
          end if
     end do
2    continue
     generateRV = min
     end function
!    ******************************************************************
!    ******************************************************************
!                                   Direct ROC regression
!    ******************************************************************
!    ******************************************************************
     subroutine ROCDirectab(Z0,X0,W0,n0,Z1,X1,W1,n1, &     ! Input
          t, nt, &                                         ! FPF
          family, &                                        ! logit or probit
          nvarz, mode0z, &                                 ! Number of vbles
          nparm, IIm, nparv, IIv, nparr, IIr, &            ! Partials in Mean, Variance and ROC
          npart, IIt, &                                    ! Fit args
          kbin,pm,pv,hm,hv,hroc,iopt, &                    ! Fit args
          parz, &                                          ! Test
          tb,ntb, &                                        ! Prediction
          pvalue)                                          ! Output
     implicit none
     integer, parameter :: nboot=400, ntaux=1
     integer n0, n1, kbin, nt, family, ntb, i, iopt, &
          iboot,parz, & 
          nvarz, mode0z(nvarz), &
          nparm, IIm(2,nparm), &
          nparv, IIv(2,nparv), &
          nparr, IIr(2,nparr), &
          npart, IIt(npart), &
          pm(nparm), pv(nparv), &
          ir0(n0)

     double precision Z0(n0,nvarz), X0(n0), W0(n0), Z1(n1,nvarz), &
          X1(n1), W1(n1), &  
          hm(nparm), hv(nparv), t(nt), tb(ntb), &
          hroc(nparr+npart+1), &
          pvalue(2), &
          taux(ntaux)
     
     double precision,external:: cInvROC,QQ,generateRV

     double precision,external::SD, cAUC,Med, Var
     
     double precision, allocatable:: M0(:), V0(:), ROC(:,:), AUC(:), &
          Fp(:,:),&
          M00(:), V00(:), ROC0(:,:), AUC0(:), &
          Fp0(:,:), &
          hroc0(:), t0(:), &
          ROCs(:,:), Fps(:,:), &
          hrocs(:),Test(:), Testb(:), &
          hrocs2(:), M10(:), V10(:), &
          Err0(:), X0b(:), X1b(:), SH(:), AUCs(:), &
          coeffb(:)

     allocate (M0(n0),V0(n0), &
          M00(n0), V00(n0), & 
          ROC0(ntb,n1), AUC0(n1), &
          ROC(ntaux,n1), AUC(n1), &
          ROCs(ntaux,n1), AUCs(n1), &
          Fp(ntaux*n1,nparr+npart+1), &
          Fps(ntaux*n1,nparr+npart+1), &
          Fp0(n1*ntb,nparr+npart+1), &
          hroc0(nparr+npart+1), &
          hrocs(nparr+npart+1), &
          hrocs2(nparr+npart+1), &
          Test(20),Testb(20), &
          t0(n1), &
          M10(n1), V10(n1), &
          Err0(n0), X0b(n0), X1b(n1), SH(n1), &
          coeffb(20))

     hroc0=hroc 
     hrocs2=hroc

     Test=0
     taux(1)=0.0
!    ****************************************************************
!     Original model
!    ****************************************************************
          call ROCDirecta(Z0,X0,W0,n0,Z1,X1,W1,n1,t,nt,family, &
                    nvarz,mode0z, &
                    nparm,IIm,nparv,IIv,nparr,IIr, &
                    npart,IIt, &
                    kbin,pm,pv,hm,hv,hroc,iopt, &
                    Z1,n1,taux,ntaux, &
                    M0,V0,ROC,AUC,Fp,coeffb)
          ! Tests
          do i=1,ntaux*n1
               Test(1) = Test(1)+ (Fp(i,parz))**2
               Test(2) = Test(2)+ abs(Fp(i,parz))
          end do
     !****************************************************************
     !      Under H0
     !****************************************************************          
          hroc0(parz)=-3
          call ROCDirecta(Z0,X0,W0,n0,Z1,X1,W1,n1, t,nt, family, & 
                    nvarz,mode0z, &
                    nparm,IIm,nparv,IIv,nparr,IIr, &
                    npart,IIt, &
                    kbin,pm,pv,hm,hv,hroc0,iopt, &
                    Z1,n1,tb,ntb, &
                    M00,V00,ROC0,AUC0,Fp0,coeffb)
          call LocScaleGam(Z0,X0,W0,n0, &
               nvarz,mode0z, nparm,IIm, nparv, IIv, hm, hv, pm, pv, &
               kbin, M0, V0, Z1, M10, V10, n1)
     !****************************************************************
     !     Resampling under H0
     !****************************************************************
     do i=1, n0
          Err0(i)=(X0(i)-M0(i))/sqrt(V0(i))
     end do     
     do iboot=1,nboot
          ! Healthy individuals
          call sample_int(n0,n0,ir0)
          do i=1,n0
               X0b(i)=M0(i) + sqrt(V0(i))*Err0(int(ir0(i)))
          end do
          ! Diseased individuals
          do i=1,n1
               t0(i) = generateRV(tb, ROC0(1,i),ntb)
          end do
          call SH_(t0, n1, Err0, W0, n0, SH)
          do i=1,n1
               X1b(i) = M10(i) + sqrt(V10(i))*SH(i)
          end do
          hrocs = hrocs2
          call ROCDirecta(Z0,X0b,W0,n0,Z1,X1b,W1,n1, t,nt, family, & 
                    nvarz,mode0z, &
                    nparm,IIm,nparv,IIv,nparr,IIr,  &
                    npart,IIt, & 
                    kbin,pm,pv,hm,hv,hrocs,iopt, &
                    Z1,n1,taux,ntaux, &
                    M00,V00,ROCs,AUCs,Fps,coeffb)

          Testb = 0
          do i=1,ntaux*n1
               Testb(1) = Testb(1) + (Fps(i,parz))**2
               Testb(2) = Testb(2) + abs(Fps(i,parz))
          end do
          !p-valor
          if (Testb(1).ge.Test(1)) then
               pvalue(1)=pvalue(1)+1.0/nboot
          end if
          if (Testb(2).ge.Test(2)) then
               pvalue(2)=pvalue(2)+1.0/nboot
          end if
     end do
     deallocate (M0, V0, ROC, AUC, Fp,  & 
          M00, V00, ROC0, AUC0, Fp0,  &
          hroc0, t0, ROCs, Fps,  & 
          hrocs, Test, Testb, hrocs2, M10, V10,  &
          Err0, X0b, X1b, SH, AUCs,  &
          coeffb)
     end subroutine
!    ******************************************************************
!    ******************************************************************
     subroutine ROCDirecta(Z0,X0,W0,n0,Z1,X1,W1,n1, &
          t, nt, family, &
          nvarz, mode0z, &
          nparm, IIm, nparv, IIv, nparr, IIr, &
          npart, IIt, &
          kbin,pm,pv, hm, hv, hroc,iopt, &
          Zb,nb,tb,ntb, &
          M0, V0, ROC, AUC, Fp, coeff)
     implicit none
     integer n0, n1, kbin, nb, nt, ntb, family, &
          i, j, it, iopt, &
          nvarz, mode0z(nvarz), mode0(nvarz+1), &
          nparm, IIm(2,nparm), &
          nparv, IIv(2,nparv), &
          nparr, IIr(2,nparr), &
          npart, IIt(npart), &
          pm(nparm), pv(nparv),&
          II(2,nparr+npart+1)

     double precision Z0(n0,nvarz), X0(n0), W0(n0), Z1(n1,nvarz), &
          X1(n1), W1(n1), hm(nparm), hv(nparv), &
          M0(n0), V0(n0), Zb(nb,nvarz), tb(ntb), t(nt), &
          ROC(ntb,nb), AUC(nb), M10(n1), V10(n1), &
          Err0(n0), &     
          pPV(n1), &
          Fp(nb*ntb, nparr+npart+1), &
          hroc(nparr+npart+1), aux, &
          coeff(20), vmin, vmax

     double precision,external:: SD, cAUC,Med, Var, normdev

     double precision, allocatable:: ZROC(:,:), XROC(:), WROC(:)     
     allocate (ZROC(n1*nt,nvarz+1), XROC(n1*nt), WROC(n1*nt))
     
     vmin = 0.0001
     vmax = 0.9999

     !*******************************
     !     Healthy individuals
     !*******************************
     call LocScaleGam(Z0,X0,W0,n0, &
               nvarz,mode0z, nparm,IIm, nparv, IIv, hm, hv, pm, pv, &
               kbin, M0, V0, Z1, M10, V10, n1)
     !*******************************
     !      Placement values
     !*******************************
     do i=1,n0          
          Err0(i)=(X0(i)-M0(i))/sqrt(V0(i))
     end do
     do i=1,n1
          pPV(i) = SD((X1(i)-M10(i))/sqrt(V10(i)),Err0, W0,n0)
          do it=1,nt
               do j=1,nvarz
                    ZROC((i-1)*nt+it,j)=Z1(i,j)
               end do
               if(iopt.eq.0) then
                    ZROC((i-1)*nt+it,nvarz+1)=t(it)
               else
                    aux=t(it)
                    aux = max(vmin, aux)
                    aux = min(vmax, aux)
                    ZROC((i-1)*nt+it,nvarz+1)=normdev(aux)
               end if
               if (pPV(i).gt.t(it)) then
                    XROC((i-1)*nt+it)= 0     
               else
                    XROC((i-1)*nt+it)= 1
               end if
               WROC((i-1)*nt+it)= W1(i)
          end do
     end do
     !*******************************
     !     ROC-GAM fit
     !*******************************
     do i =1,nvarz
          mode0(i)=mode0z(i)
     end do
     mode0(nvarz+1) = 5     ! FPF
     
     ! Covariates
     do i=1,nparr
          do j=1,2
               II(j,i)=IIr(j,i)
          end do
     end do

     ! Main effect (FPF)
     II(1,nparr+1)= -1
     II(2,nparr+1)= nvarz+1

     ! Possible interactions with FPF (not implemented)
     if(npart.gt.0) then
          do i=1,npart
               II(1,nparr+i+1)=nvarz+1
               II(2,nparr+i+1)=IIt(i)
          end do          
     end if
     call GAMROC(ZROC, XROC, WROC, nt*n1, family, &
          hroc, nvarz+1, nparr+npart+1,II, mode0, iopt,kbin, &
          Zb,nb,tb,ntb, &
          ROC, Fp, coeff)
     do i=1,nb
          AUC(i) = cAUC(ROC(1,i),tb,ntb)
     end do
     deallocate (ZROC, XROC, WROC)
     end
!    ******************************************************************
!    ******************************************************************
     subroutine GAMROC(ZROC, XROC, WROC, nroc, &
          family, hroc, nvar, npar, II, mode0, &
          iopt,kbin, Zb, nb, t, nt, &     
          ROC, Fp, coeff)
     implicit none
     integer i,j,it, nvar, npar, nroc, kbin, nb, nt,II(2,npar), &
          mode0(nvar), iopt, proc(npar), family 
     double precision ZROC(nroc,nvar),XROC(nroc), WROC(nroc), &
          Zb(nb,nvar-1), t(nt), ROC(nt,nb), hroc(npar), &
          Fp(nb*nt,npar), aux, coeff(20), vmax, vmin
     double precision, external :: normdev
     double precision,allocatable::F(:,:), muhat(:), muhatp(:), &
          Zp(:,:), hrocn(:,:)
     allocate (F(nroc,npar), muhat(nroc), muhatp(nb*nt), &
          Zp(nb*nt, nvar), hrocn(nroc,npar))
     ! Initialize
     proc=1
     vmax = 0.9999
     vmin = 0.0001
     do i=1,nroc
          do j=1,npar
               hrocn(i,j) = hroc(j)
          end do
     end do
     ! Prediction matrix
     do i=1,nb
          do it=1,nt
               do j=1,nvar-1 ! "Proper" covariates
                    Zp((i-1)*nt+it,j)=Zb(i,j)
               end do
               if(iopt.eq.0) then ! FPF (smooth)
                    Zp((i-1)*nt+it,nvar)=t(it)
               else
                    aux = t(it)
                    aux = max(vmin, aux)
                    aux = min(vmax, aux)
                    Zp((i-1)*nt+it,nvar)=normdev(aux)
               end if
          enddo
     end do
     call GAM(nroc, nvar+1, npar, mode0, II, ZROC, XROC, WROC, &
           hrocn, kbin, proc, family, F, coeff, muhat, Zp, Fp, &
          muhatp, nb*nt)
     ! ROC curve
     do i = 1,nb
          do it = 1,nt
               ROC(it,i) = muhatp((i-1)*nt+it)
          end do
     end do
     deallocate(F,muhat,Zp, muhatp)
     end
!    ******************************************************************
!    ******************************************************************
!                                   ROC INDUCED
!    ******************************************************************
!    ******************************************************************
     subroutine ROCInducedb(Z0,X0,W0,n0,Z1,X1,W1,n1, &
          kbin,p,h,t,nt,pvalue)
     implicit none
     integer, parameter::nboot=400, nb=100
     integer nt,n0,n1,i,j,kbin,p(2,2),iboot,ir0(n0), &
          ir1(n1)
     logical optaccROC, optaccYI, optacc(2),opttest
     double precision X0(n0),X1(n1),Z0(n0),Z1(n1),W0(n0),W1(n1), &
          zmin,zmax, h(2,2),t(nt), &
          Test(20),Tb(20), &
          pvalue(20), &
          Zp(n0+n1), Wp(n0+n1), Wb(nb)

     double precision, allocatable:: ROC(:,:), AUC(:), &
          M0(:),M1(:),V0(:),V1(:), &
          M0p(:),M1p(:),V0p(:),V1p(:), &
          Err0(:),Err1(:), &
          SH(:),M1c(:),V1c(:),M0c(:),V0c(:), &
          M0g(:),V0g(:),M1g(:),V1g(:), Err0g(:), Err1g(:), &
          X0b(:), X1b(:),AROC(:), &
          YI(:),TH(:), Zb(:)

     double precision,external:: SD, cAUC, minimum, maximum

     allocate (ROC(nt,nb),AUC(nb), &
          M0(n0),M1(n1),V0(n0),V1(n1), &
          M0p(nb),M1p(nb),V0p(nb),V1p(nb), &
          Err0(n0),Err1(n1), &
          SH(nt),M1c(n1),V1c(n1),M0c(n1),V0c(n1), &
          M0g(n0),V0g(n0),M1g(n1),V1g(n1), Err0g(n0), Err1g(n1), &
          X0b(n0), X1b(n1),AROC(nt), &
          YI(nb),TH(nb),Zb(nb))

     ! Initialize options
     optacc(1) = .FALSE.
     optacc(2) = .FALSE.
     optaccROC = .FALSE. 
     optaccYI = .FALSE. 
     opttest = .TRUE.

     zmin = min(minimum(Z0,n0),minimum(Z1,n1))
     zmax = max(maximum(Z0,n0),maximum(Z1,n1))
     do i=1,nb          
          Zb(i)=zmin+(i-1)*(zmax-zmin)/(nb-1)
     end do
     call ROCInduced(Z0,X0,W0,n0,Z1,X1,W1,n1,kbin,p,h, &
          Zb,nb,t,nt, &
          M0,V0,M1,V1, &
          M0p,V0p,M1p,V1p, &
          ROC,AUC,AROC, &
          optaccROC, optaccYI, optacc, YI, TH, &
          opttest,M1c,V1c)
     !*****************************
     ! Weights
     !*****************************
     do i=1,n0
          Zp(i)=Z0(i)
          Wp(i)=W0(i)
     end do
     do i=1,n1
          Zp(i+n0)=Z1(i)
          Wp(i+n0)=W1(i)
     end do
     call sampleBinning(Zp,n0+n1,Wp,Zb,nb,Wb)
     !*****************************
     ! Test
     !*****************************
     Test=0
     do i=1,nb
          do j=1,nt
               Test(1)=Test(1)+Wb(i)*abs(ROC(j,i)-AROC(j))/nt     
          end do
     end do
     !*****************************
     ! Pilot estimates
     !*****************************
     M0g=M0
     M1g=M1c
     V0g=V0
     V1g=V1c
     !*****************************
     ! Residuals under H0     
     !*****************************
     Err0g=(X0-M0g)/sqrt(V0g)
     Err1g=(X1-M1g)/sqrt(V1g)

     pvalue=0
     do iboot=1,nboot
          opttest=.FALSE.
          call sample_int(n0,n0,ir0)
          do i=1,n0
               X0b(i)=M0g(i)+sqrt(V0g(i))*Err0g(ir0(i))
          end do
          call sample_int(n1,n1,ir1)
          do i=1,n1
               X1b(i)=M1g(i)+sqrt(V1g(i))*Err1g(ir1(i))
          end do
          call ROCInduced(Z0,X0b,W0,n0,Z1,X1b,W1,n1, &
          kbin,p,h, &
          Zb,nb,t,nt, &
          M0,V0,M1,V1, &
          M0p,V0p,M1p,V1p, &
          ROC,AUC,AROC, &
          optaccROC,optaccYI,optacc,YI,TH, &
          opttest,M1c,V1c)
          Tb=0
          do i=1,nb
               do j=1,nt
                    Tb(1)=Tb(1)+Wb(i)*abs(ROC(j,i)-AROC(j))/nt     
               end do
          end do
          ! p-value
          if (Tb(1).ge.Test(1)) then
               pvalue(1)=pvalue(1)+1.0/nboot
          end if
     end do

     deallocate (ROC, AUC, &
          M0,M1,V0,V1, &
          M0p, M1p, V0p, V1p, &
          Err0, Err1, &
          SH,M1c,V1c,M0c,V0c, &
          M0g,V0g,M1g,V1g, Err0g, Err1g, &
          X0b, X1b,AROC, &
          YI,TH,Zb)
     end
!    ******************************************************************
!    ******************************************************************
     subroutine ROCInduced( &
          Z0,X0,W0,n0,Z1,X1,W1,n1,kbin,p,h, &
          Zp, np,t,nt, &
          M0,V0,M1,V1, &
          M0p,V0p,M1p,V1p, &
          ROC,AUC, AROC, &
          optaccROC, &
          optaccYI, &
          optacc,YI,TH, &
          opttest,M1c,V1c)
     implicit none     
     integer, parameter::ntyi = 1000, nvar = 1, npar = 1
     integer n0,n1,i,j,kbin, p(2,2), np, nt, maxt
     logical optacc(2), opttest, optaccROC, &
          optaccYI
     double precision X0(n0),X1(n1),Z0(n0),Z1(n1),W0(n0),W1(n1), &
          Zp(np),t(nt),ROC(nt,np),AUC(np),h(2,2), &
          M1c(n1),V1c(n1),AROC(nt), &
          YI(np),TH(np), A, B, &
          temp,temp1, &
          M0(n0),M1(n1),V0(n0),V1(n1), &
          M0p(np),M1p(np),V0p(np),V1p(np)
     integer, allocatable :: II0(:,:), mode0(:), II1(:,:)
     double precision, allocatable:: Err0(:),Err1(:), &
          SH(:),M0c(:),V0c(:), &
          M01(:),V01(:), &
          ROC_YI(:,:),SH_YI(:),AROCTH(:), &
          tyi(:)
     double precision,external::csval, SD, cAUC, minimum, maximum
     allocate (Err0(n0),Err1(n1), &
          SH(nt),M01(n1),V01(n1), &
          II0(2,npar),mode0(npar),II1(2,npar))

     mode0 = 5
     !******************************************************************
     !     Healthy individuals
     !******************************************************************     
     II0(1,1)= -1
     II0(2,1)= 1
     call LocScaleGam(Z0, X0, W0, n0, nvar, mode0, &
          npar, II0, npar, II0, &
          h(1,1), h(2,1), p(1,1), p(2,1), &
          kbin, M0, V0, Zp, M0p, V0p, np)     
     !******************************************************************
     !     Diseased individuals
     !******************************************************************
     II1(1,1)= -1
     II1(2,1)= 1
     call LocScaleGam(Z1, X1, W1, n1, nvar, mode0, &
          npar, II0, npar, II0, &
          h(1,2), h(2,2), p(1,2), p(2,2), &
          kbin, M1, V1, Zp, M1p, V1p, np) 
     !******************************************************************
     !     Standardized residuals
     !******************************************************************
      do i=1,n0
          Err0(i)=(X0(i)-M0(i))/sqrt(V0(i))
     end do
      do i=1,n1
          Err1(i)=(X1(i)-M1(i))/sqrt(V1(i))
     end do     
     !******************************************************************
     !     ROC curve and AUC
     !******************************************************************
     call SH_(t,nt,Err0,W0,n0,SH)
     do i=1,np
          call cROC(M0p(i),M1p(i),V0p(i), V1p(i), Err1, W1, n1, &
                nt, SH, ROC(1,i))
          AUC(i) = cAUC(ROC(1,i),t,nt)
     end do
     !******************************************************************
     !     Adjusted ROC curve
     !******************************************************************          
     AROC = 0
     call LocScaleGam(Z0, X0, W0, n0, nvar, mode0, &
          npar, II0, npar, II0, &
          h(1,1), h(2,1), p(1,1), p(2,1), &
          kbin, M0, V0, Z1, M01, V01, n1) 
     do i=1,nt
          AROC(i)=0.0
          do j=1,n1
               if((X1(j)-M01(j))/sqrt(V01(j)).gt.SH(i)) then
                    AROC(i)=AROC(i)+1.0
               end if
          end do
          AROC(i) = AROC(i)/(1.0*n1)
     end do          
     !******************************************************************
     !     YI, EQ and associated cut-points
     !******************************************************************
     YI = 0
     TH = 0
     if(optacc(1).or.optacc(2)) then
          allocate(SH_YI(ntyi),tyi(ntyi))
          do i=1,ntyi
               tyi(i) = (i-1)*1.0/(ntyi-1)
          end do
          call SH_(tyi,ntyi,Err0,W0,n0,SH_YI)
          
          if (optaccROC) then
               allocate(ROC_YI(ntyi,np))
               do i=1,np
                    call cROC(M0p(i),M1p(i),V0p(i), V1p(i), Err1, W1, n1, &
                         ntyi, SH_YI, ROC_YI(1,i))          
                    if(optaccYI) then
                         maxt = 1
                         YI(i) = 0
                         do j = 1, ntyi
                              temp = abs(ROC_YI(j,i)-tyi(j))
                              if(temp .gt. YI(i)) then 
                                   YI(i) = temp
                                   maxt = j
                              end if
                         end do
                    else
                         maxt=1
                         temp1=1
                         do j = 1, ntyi
                              temp=abs(ROC_YI(j,i)-1+tyi(j))
                              if(temp.lt.temp1) then
                                   temp1=temp
                                   maxt=j
                              end if
                         end do
                         YI(i)=1.0-tyi(maxt)
                    end if
                    if(optacc(2)) then
                         TH(i) = M0p(i) + sqrt(V0p(i))*SH_YI(maxt)
                    end if
               end do
               deallocate(ROC_YI)
          else
               allocate(AROCTH(ntyi))
               AROCTH = 0
               do i=1,ntyi
                    do j=1,n1
                         if((X1(j)-M01(j))/sqrt(V01(j)).gt.SH_YI(i)) then
                              AROCTH(i) = AROCTH(i) + 1.0
                         end if
                    end do
                    AROCTH(i) = AROCTH(i)/(1.0*n1)
               end do
               if(optaccYI) then
                    maxt=1
                    temp1=0
                    do j = 1, ntyi
                         temp = abs(AROCTH(j)-tyi(j))
                         if(temp.gt.temp1) then 
                              temp1 = temp
                              maxt = j
                         end if
                    end do
                    YI = temp1
               else
                    maxt=1
                    temp1=1
                    do j = 1, ntyi
                         temp=abs(AROCTH(j)-1+tyi(j))
                         if(temp.lt.temp1) then
                              temp1=temp
                              maxt=j
                         end if
                    end do
                    YI = 1.0 - tyi(maxt)
               end if
               if(optacc(2)) then
                    do  i=1,np
                         TH(i) = M0p(i)+sqrt(V0p(i))*SH_YI(maxt)
                    end do
               end if
               deallocate(AROCTH)     
          end if
          deallocate(SH_YI,tyi)
     end if 
     !******************************************************************
     ! Test
     !******************************************************************
     M1c = 0
     V1c = 0
     B = 1.0
     A = 1.0          
     if(opttest) then
          allocate(M0c(n1),V0c(n1))
          do i=1,n1
               M0c(i) = M01(i)
               V0c(i) = V01(i)
          end do
          do i=1,n1
               V1c(i) = (B**2)*V0c(i)
          end do
          do i=1,n1
               M1c(i) = M0c(i)-A*sqrt(V1c(i))
          end do
          deallocate(M0c, V0c)
     end if
     deallocate (Err0, Err1, SH,& 
          M01, V01, &
          II0, II1, mode0)
     end
!    ******************************************************************
!    ******************************************************************
     MODULE lsq
     !  Module for unconstrained linear least-squares calculations.
     !  The algorithm is suitable for updating LS calculations as more
     !  data are added.   This is sometimes called recursive estimation.
     !  Only one dependent variable is allowed.
     !  Based upon Applied Statistics algorithm AS 274.
     !  Translation from Fortran 77 to Fortran 90 by Alan Miller.
     !  A function, VARPRD, has been added for calculating the variances
     !  of predicted values, and this uses a subroutine BKSUB2.

     !  Version 1.14, 19 August 2002 - ELF90 compatible version
     !  Author: Alan Miller
     !  e-mail : amiller @ bigpond.net.au
     !  WWW-pages: http://www.ozemail.com.au/~milleraj
     !             http://users.bigpond.net.au/amiller/

     !  Bug fixes:
     !  1. In REGCF a call to TOLSET has been added in case the user had
     !     not set tolerances.
     !  2. In SING, each time a singularity is detected, unless it is in the
     !     variables in the last position, INCLUD is called.   INCLUD assumes
     !     that a new observation is being added and increments the number of
     !     cases, NOBS.   The line:  nobs = nobs - 1 has been added.
     !  3. row_ptr was left out of the DEALLOCATE statement in routine startup
     !     in version 1.07.
     !  4. In COV, now calls SS if rss_set = .FALSE.  29 August 1997
     !  5. In TOLSET, correction to accomodate negative values of D.  19 August 2002

     !  Other changes:
     !  1. Array row_ptr added 18 July 1997.   This points to the first element
     !     stored in each row thus saving a small amount of time needed to
     !     calculate its position.
     !  2. Optional parameter, EPS, added to routine TOLSET, so that the user
     !     can specify the accuracy of the input data.
     !  3. Cosmetic change of lsq_kind to dp (`Double precision')
     !  4. Change to routine SING to use row_ptr rather than calculate the position
     !     of first elements in each row.

     !  The PUBLIC variables are:
     !  dp       = a KIND parameter for the floating-point quantities calculated
     !             in this module.   See the more detailed explanation below.
     !             This KIND parameter should be used for all floating-point
     !             arguments passed to routines in this module.

     !  nobs    = the number of observations processed to date.
     !  ncol    = the total number of variables, including one for the constant,
     !            if a constant is being fitted.
     !  r_dim   = the dimension of array r = ncol*(ncol-1)/2
     !  vorder  = an integer vector storing the current order of the variables
     !            in the QR-factorization.   The initial order is 0, 1, 2, ...
     !            if a constant is being fitted, or 1, 2, ... otherwise.
     !  initialized = a logical variable which indicates whether space has
     !                been allocated for various arrays.
     !  tol_set = a logical variable which is set when subroutine TOLSET has
     !            been called to calculate tolerances for use in testing for
     !            singularities.
     !  rss_set = a logical variable indicating whether residual sums of squares
     !            are available and usable.
     !  d()     = array of row multipliers for the Cholesky factorization.
     !            The factorization is X = Q.sqrt(D).R where Q is an ortho-
     !            normal matrix which is NOT stored, D is a diagonal matrix
     !            whose diagonal elements are stored in array d, and R is an
     !            upper-triangular matrix with 1's as its diagonal elements.
     !  rhs()   = vector of RHS projections (after scaling by sqrt(D)).
     !            Thus Q'y = sqrt(D).rhs
     !  r()     = the upper-triangular matrix R.   The upper triangle only,
     !            excluding the implicit 1's on the diagonal, are stored by
     !            rows.
     !  tol()   = array of tolerances used in testing for singularities.
     !  rss()   = array of residual sums of squares.   rss(i) is the residual
     !            sum of squares with the first i variables in the model.
     !            By changing the order of variables, the residual sums of
     !            squares can be found for all possible subsets of the variables.
     !            The residual sum of squares with NO variables in the model,
     !            that is the total sum of squares of the y-values, can be
     !            calculated as rss(1) + d(1)*rhs(1)^2.   If the first variable
     !            is a constant, then rss(1) is the sum of squares of
     !            (y - ybar) where ybar is the average value of y.
     !  sserr   = residual sum of squares with all of the variables included.
     !  row_ptr() = array of indices of first elements in each row of R.
     !
     !--------------------------------------------------------------------------
     !     General declarations
     !--------------------------------------------------------------------------
     IMPLICIT NONE
     INTEGER, SAVE                :: nobs, ncol, r_dim
     INTEGER, ALLOCATABLE, SAVE   :: vorder(:), row_ptr(:)
     LOGICAL, SAVE                :: initialized = .false., &
                                               tol_set = .false., rss_set = .false.
     ! Note. dp is being set to give at least 12 decimal digit
     !       representation of floating point numbers.   This should be adequate
     !       for most problems except the fitting of polynomials.   dp is
     !       being set so that the same code can be run on PCs and Unix systems,
     !       which will usually represent floating-point numbers in `double
     !       precision', and other systems with larger word lengths which will
     !       give similar accuracy in `single precision'.

     INTEGER, PARAMETER  :: dp = SELECTED_REAL_KIND(12,60)
     double precision, ALLOCATABLE :: d(:), rhs(:), r(:), tol(:), & 
                                                  rss(:)
     double precision              :: zero = 0.0_dp, one = 1.0_dp, &
                                                  vsmall
     double precision              :: sserr, toly

     PUBLIC                 :: dp, nobs, ncol, r_dim, vorder, row_ptr, &
                                 initialized, tol_set, rss_set, &
                                 d, rhs, r, tol, rss, sserr
     PRIVATE                :: zero, one, vsmall
     CONTAINS
!    ***********************************************************************************
!    ***********************************************************************************
     SUBROUTINE startup(nvar, fit_const)
     !     Allocates dimensions for arrays and initializes to zero
     !     The calling program must set nvar = the number of variables, and
     !     fit_const = .true. if a constant is to be included in the model,
     !     otherwise fit_const = .false.
     !
     !--------------------------------------------------------------------------
     IMPLICIT NONE
     INTEGER, INTENT(IN)  :: nvar
     LOGICAL, INTENT(IN)  :: fit_const

     !     Local variable
     INTEGER   :: i

     vsmall = 10. * TINY(zero)

     nobs = 0
     IF (fit_const) THEN
       ncol = nvar + 1
     ELSE
       ncol = nvar
     END IF

     IF (initialized) DEALLOCATE(d, rhs, r, tol, rss, vorder, row_ptr)
     r_dim = ncol * (ncol - 1)/2
     ALLOCATE( d(ncol), rhs(ncol), r(r_dim), tol(ncol), rss(ncol), &
               vorder(ncol), row_ptr(ncol))

     d = zero
     rhs = zero
     r = zero
     sserr = zero

     IF (fit_const) THEN
       DO i = 1, ncol
         vorder(i) = i-1
       END DO
     ELSE
       DO i = 1, ncol
         vorder(i) = i
       END DO
     END IF ! (fit_const)

     ! row_ptr(i) is the position of element R(i,i+1) in array r().

     row_ptr(1) = 1
     DO i = 2, ncol-1
       row_ptr(i) = row_ptr(i-1) + ncol - i + 1
     END DO
     row_ptr(ncol) = 0

     initialized = .true.
     tol_set = .false.
     rss_set = .false.

     RETURN
     END SUBROUTINE startup
!    ***********************************************************************************
!    ***********************************************************************************
     SUBROUTINE endup()
          IF (initialized) DEALLOCATE(d, rhs, r, tol, rss, &
                                        vorder, row_ptr)
          initialized = .FALSE.
     RETURN
     END SUBROUTINE endup
!    ***********************************************************************************
!    ***********************************************************************************
     SUBROUTINE includ(weight, xrow, yelem)

     !     ALGORITHM AS75.1  APPL. STATIST. (1974) VOL.23, NO. 3

     !     Calling this routine updates D, R, RHS and SSERR by the
     !     inclusion of xrow, yelem with the specified weight.

     !     *** WARNING  Array XROW is overwritten.

     !     N.B. As this routine will be called many times in most applications,
     !          checks have been eliminated.
     !
     !--------------------------------------------------------------------------
     IMPLICIT NONE
     double precision,INTENT(IN)                    :: weight, yelem
     double precision, DIMENSION(:), INTENT(IN OUT) :: xrow

     !     Local variables

     INTEGER     :: i, k, nextr
     double precision   :: w, y, xi, di, wxi, dpi, cbar, sbar, xk

     nobs = nobs + 1
     w = weight
     y = yelem
     rss_set = .false.
     nextr = 1
     DO i = 1, ncol

     !     Skip unnecessary transformations.   Test on exact zeroes must be
     !     used or stability can be destroyed.

       IF (ABS(w) < vsmall) RETURN
       xi = xrow(i)
       IF (ABS(xi) < vsmall) THEN
         nextr = nextr + ncol - i
       ELSE
         di = d(i)
         wxi = w * xi
         dpi = di + wxi*xi
         cbar = di / dpi
         sbar = wxi / dpi
         w = cbar * w
         d(i) = dpi
         DO k = i+1, ncol
               xk = xrow(k)
               xrow(k) = xk - xi * r(nextr)
               r(nextr) = cbar * r(nextr) + sbar * xk
               nextr = nextr + 1
         END DO
         xk = y
         y = xk - xi * rhs(i)
         rhs(i) = cbar * rhs(i) + sbar * xk
       END IF
     END DO ! i = 1, ncol

     !     Y * SQRT(W) is now equal to the Brown, Durbin & Evans recursive
     !     residual.

     sserr = sserr + w * y * y

     RETURN
     END SUBROUTINE includ
!    ***********************************************************************************
!    ***********************************************************************************
     SUBROUTINE regcf(beta, nreq, ifault)

     !     ALGORITHM AS274  APPL. STATIST. (1992) VOL.41, NO. 2

     !     Modified version of AS75.4 to calculate regression coefficients
     !     for the first NREQ variables, given an orthogonal reduction from
     !     AS75.1.
     !
     !--------------------------------------------------------------------------
     IMPLICIT NONE
     INTEGER, INTENT(IN)                  :: nreq
     INTEGER, INTENT(OUT)                 :: ifault
     double precision, DIMENSION(:), INTENT(OUT) :: beta

     !     Local variables

     INTEGER   :: i, j, nextr

     !     Some checks.

     ifault = 0
     IF (nreq < 1 .OR. nreq > ncol) ifault = ifault + 4
     IF (ifault /= 0) RETURN

     IF (.NOT. tol_set) CALL tolset()

     DO i = nreq, 1, -1
       IF (SQRT(d(i)) < tol(i)) THEN
         beta(i) = zero
         d(i) = zero
         ifault = -i
       ELSE
         beta(i) = rhs(i)
         nextr = row_ptr(i)
         DO j = i+1, nreq
               beta(i) = beta(i) - r(nextr) * beta(j)
               nextr = nextr + 1
         END DO ! j = i+1, nreq
       END IF
     END DO ! i = nreq, 1, -1

     RETURN
     END SUBROUTINE regcf
!    ***********************************************************************************
!    ***********************************************************************************
     SUBROUTINE tolset(eps)

     !     ALGORITHM AS274  APPL. STATIST. (1992) VOL.41, NO. 2

     !     Sets up array TOL for testing for zeroes in an orthogonal
     !     reduction formed using AS75.1.
     double precision, INTENT(IN), OPTIONAL :: eps

     !     Unless the argument eps is set, it is assumed that the input data are
     !     recorded to full machine accuracy.   This is often not the case.
     !     If, for instance, the data are recorded to `single precision' of about
     !     6-7 significant decimal digits, then singularities will not be detected.
     !     It is suggested that in this case eps should be set equal to
     !     10.0 * EPSILON(1.0)
     !     If the data are recorded to say 4 significant decimals, then eps should
     !     be set to 1.0E-03
     !     The above comments apply to the predictor variables, not to the
     !     dependent variable.

     !     Correction - 19 August 2002
     !     When negative weights are used, it is possible for an alement of D
     !     to be negative.
     !--------------------------------------------------------------------------
     !     Local variables
     !--------------------------------------------------------------------------
     INTEGER    :: col, row, pos
     double precision  :: eps1, ten = 10.0, total, work(ncol)

     !     EPS is a machine-dependent constant.

     IF (PRESENT(eps)) THEN
       eps1 = MAX(ABS(eps), ten * EPSILON(ten))
     ELSE
       eps1 = ten * EPSILON(ten)
     END IF
     work = SQRT(ABS(d))
     DO col = 1, ncol
       pos = col - 1
       total = work(col)
       DO row = 1, col-1
         total = total + ABS(r(pos)) * work(row)
         pos = pos + ncol - row - 1
       END DO
       tol(col) = eps1 * total
     END DO

     tol_set = .TRUE.
     RETURN
     END SUBROUTINE tolset
!    ***********************************************************************************
!    ***********************************************************************************
     SUBROUTINE sing(lindep, ifault)

     !     ALGORITHM AS274  APPL. STATIST. (1992) VOL.41, NO. 2

     !     Checks for singularities, reports, and adjusts orthogonal
     !     reductions produced by AS75.1.

     !     Correction - 19 August 2002
     !     When negative weights are used, it is possible for an alement of D
     !     to be negative.

     !     Auxiliary routines called: INCLUD, TOLSET
     !
     !--------------------------------------------------------------------------
     INTEGER, INTENT(OUT)                :: ifault
     LOGICAL, DIMENSION(:), INTENT(OUT)  :: lindep

     !     Local variables

     double precision  :: temp, x(ncol), work(ncol), y, weight
     INTEGER    :: pos, row, pos2

     ifault = 0

     work = SQRT(ABS(d))
     IF (.NOT. tol_set) CALL tolset()

     DO row = 1, ncol
       temp = tol(row)
       pos = row_ptr(row)         ! pos = location of first element in row

     !     If diagonal element is near zero, set it to zero, set appropriate
     !     element of LINDEP, and use INCLUD to augment the projections in
     !     the lower rows of the orthogonalization.

       lindep(row) = .FALSE.
       IF (work(row) <= temp) THEN
         lindep(row) = .TRUE.
         ifault = ifault - 1
         IF (row < ncol) THEN
          pos2 = pos + ncol - row - 1
          x = zero
          x(row+1:ncol) = r(pos:pos2)
          y = rhs(row)
          weight = d(row)
          r(pos:pos2) = zero
          d(row) = zero
          rhs(row) = zero
          CALL includ(weight, x, y)
          ! INCLUD automatically increases the number
          ! of cases each time it is called.
          nobs = nobs - 1
         ELSE
          sserr = sserr + d(row) * rhs(row)**2
         END IF ! (row < ncol)
       END IF ! (work(row) <= temp)
     END DO ! row = 1, ncol

     RETURN
     END SUBROUTINE sing
!    ***********************************************************************************
!    ***********************************************************************************
     SUBROUTINE ss()

     !     ALGORITHM AS274  APPL. STATIST. (1992) VOL.41, NO. 2

     !     Calculates partial residual sums of squares from an orthogonal
     !     reduction from AS75.1.
     !
     !--------------------------------------------------------------------------

     !     Local variables
     INTEGER    :: i
     double precision  :: total

     total = sserr
     rss(ncol) = sserr
     DO i = ncol, 2, -1
       total = total + d(i) * rhs(i)**2
       rss(i-1) = total
     END DO

     rss_set = .TRUE.
     RETURN
     END SUBROUTINE ss
!    ***********************************************************************************
!    ***********************************************************************************
     SUBROUTINE cov(nreq, var, covmat, dimcov, sterr, ifault)

     !     ALGORITHM AS274  APPL. STATIST. (1992) VOL.41, NO. 2

     !     Calculate covariance matrix for regression coefficients for the
     !     first nreq variables, from an orthogonal reduction produced from
     !     AS75.1.

     !     Auxiliary routine called: INV
     !
     !--------------------------------------------------------------------------
     INTEGER, INTENT(IN)                   :: nreq, dimcov
     INTEGER, INTENT(OUT)                  :: ifault
     double precision, INTENT(OUT)                :: var
     double precision, DIMENSION(:), INTENT(OUT)  :: covmat, sterr

     !     Local variables.

     INTEGER::dim_rinv, pos, row, start, pos2, col, pos1, k
     double precision              :: total
     double precision, ALLOCATABLE :: rinv(:)

     !     Check that dimension of array covmat is adequate.

     IF (dimcov < nreq*(nreq+1)/2) THEN
       ifault = 1
       RETURN
     END IF

     !     Check for small or zero multipliers on the diagonal.

     ifault = 0
     DO row = 1, nreq
       IF (ABS(d(row)) < vsmall) ifault = -row
     END DO
     IF (ifault /= 0) RETURN

     !     Calculate estimate of the residual variance.

     IF (nobs > nreq) THEN
       IF (.NOT. rss_set) CALL ss()
       var = rss(nreq) / (nobs - nreq)
     ELSE
       ifault = 2
       RETURN
     END IF

     dim_rinv = nreq*(nreq-1)/2
     ALLOCATE ( rinv(dim_rinv) )

     CALL INV(nreq, rinv)
     pos = 1
     start = 1
     DO row = 1, nreq
       pos2 = start
       DO col = row, nreq
         pos1 = start + col - row
         IF (row == col) THEN
          total = one / d(col)
         ELSE
          total = rinv(pos1-1) / d(col)
         END IF
         DO K = col+1, nreq
               total = total + rinv(pos1) * rinv(pos2) / d(k)
               pos1 = pos1 + 1
               pos2 = pos2 + 1
         END DO ! K = col+1, nreq
         covmat(pos) = total * var
         IF (row == col) sterr(row) = SQRT(covmat(pos))
         pos = pos + 1
       END DO ! col = row, nreq
       start = start + nreq - row
     END DO ! row = 1, nreq

     DEALLOCATE(rinv)
     RETURN
     END SUBROUTINE cov
!    ***********************************************************************************
!    ***********************************************************************************
     SUBROUTINE inv(nreq, rinv)

     !     ALGORITHM AS274  APPL. STATIST. (1992) VOL.41, NO. 2

     !     Invert first nreq rows and columns of Cholesky factorization
     !     produced by AS 75.1.
     !
     !--------------------------------------------------------------------------
     INTEGER, INTENT(IN)                  :: nreq
     double precision, DIMENSION(:), INTENT(OUT) :: rinv

     !     Local variables.

     INTEGER    :: pos, row, col, start, k, pos1, pos2
     double precision  :: total

     !     Invert R ignoring row multipliers, from the bottom up.

     pos = nreq * (nreq-1)/2
     DO row = nreq-1, 1, -1
       start = row_ptr(row)
       DO col = nreq, row+1, -1
         pos1 = start
         pos2 = pos
         total = zero
         DO k = row+1, col-1
               pos2 = pos2 + nreq - k
               total = total - r(pos1) * rinv(pos2)
               pos1 = pos1 + 1
         END DO ! k = row+1, col-1
         rinv(pos) = total - r(pos1)
         pos = pos - 1
       END DO ! col = nreq, row+1, -1
     END DO ! row = nreq-1, 1, -1

     RETURN
     END SUBROUTINE inv
!    ***********************************************************************************
!    ***********************************************************************************
     SUBROUTINE partial_corr(in, cormat, dimc, ycorr, ifault)

     !     Replaces subroutines PCORR and COR of:
     !     ALGORITHM AS274  APPL. STATIST. (1992) VOL.41, NO. 2

     !     Calculate partial correlations after the variables in rows
     !     1, 2, ..., IN have been forced into the regression.
     !     If IN = 1, and the first row of R represents a constant in the
     !     model, then the usual simple correlations are returned.

     !     If IN = 0, the value returned in array CORMAT for the correlation
     !     of variables Xi & Xj is:
     !       sum ( Xi.Xj ) / Sqrt ( sum (Xi^2) . sum (Xj^2) )

     !     On return, array CORMAT contains the upper triangle of the matrix of
     !     partial correlations stored by rows, excluding the 1's on the diagonal.
     !     e.g. if IN = 2, the consecutive elements returned are:
     !     (3,4) (3,5) ... (3,ncol), (4,5) (4,6) ... (4,ncol), etc.
     !     Array YCORR stores the partial correlations with the Y-variable
     !     starting with YCORR(IN+1) = partial correlation with the variable in
     !     position (IN+1).
     !
     !--------------------------------------------------------------------------
     INTEGER, INTENT(IN)                  :: in, dimc
     INTEGER, INTENT(OUT)                 :: ifault
     double precision, DIMENSION(:), INTENT(OUT) :: cormat, ycorr

     !     Local variables.

     INTEGER    :: base_pos, pos, row, col, col1, col2, pos1, pos2
     double precision  :: rms(in+1:ncol), sumxx, sumxy, sumyy, &
                         work(in+1:ncol)
     !     Some checks.

     ifault = 0
     IF (in < 0 .OR. in > ncol-1) ifault = ifault + 4
     IF (dimc < (ncol-in)*(ncol-in-1)/2) ifault = ifault + 8
     IF (ifault /= 0) RETURN

     !     Base position for calculating positions of elements in row (IN+1) of R.

     base_pos = in*ncol - (in+1)*(in+2)/2

     !     Calculate 1/RMS of elements in columns from IN to (ncol-1).

     IF (d(in+1) > zero) rms(in+1) = one / SQRT(d(in+1))
     DO col = in+2, ncol
       pos = base_pos + col
       sumxx = d(col)
       DO row = in+1, col-1
         sumxx = sumxx + d(row) * r(pos)**2
         pos = pos + ncol - row - 1
       END DO ! row = in+1, col-1
       IF (sumxx > zero) THEN
         rms(col) = one / SQRT(sumxx)
       ELSE
         rms(col) = zero
         ifault = -col
       END IF ! (sumxx > zero)
     END DO ! col = in+1, ncol-1

     !     Calculate 1/RMS for the Y-variable

     sumyy = sserr
     DO row = in+1, ncol
       sumyy = sumyy + d(row) * rhs(row)**2
     END DO ! row = in+1, ncol
     IF (sumyy > zero) sumyy = one / SQRT(sumyy)

     !     Calculate sums of cross-products.
     !     These are obtained by taking dot products of pairs of columns of R,
     !     but with the product for each row multiplied by the row multiplier
     !     in array D.

     pos = 1
     DO col1 = in+1, ncol
       sumxy = zero
       work(col1+1:ncol) = zero
       pos1 = base_pos + col1
       DO row = in+1, col1-1
         pos2 = pos1 + 1
         DO col2 = col1+1, ncol
               work(col2) = work(col2) + d(row) * r(pos1) * r(pos2)
               pos2 = pos2 + 1
         END DO ! col2 = col1+1, ncol
         sumxy = sumxy + d(row) * r(pos1) * rhs(row)
         pos1 = pos1 + ncol - row - 1
       END DO ! row = in+1, col1-1

     !     Row COL1 has an implicit 1 as its first element (in column COL1)

       pos2 = pos1 + 1
       DO col2 = col1+1, ncol
         work(col2) = work(col2) + d(col1) * r(pos2)
         pos2 = pos2 + 1
         cormat(pos) = work(col2) * rms(col1) * rms(col2)
         pos = pos + 1
       END DO ! col2 = col1+1, ncol
       sumxy = sumxy + d(col1) * rhs(col1)
       ycorr(col1) = sumxy * rms(col1) * sumyy
     END DO ! col1 = in+1, ncol-1

     ycorr(1:in) = zero

     RETURN
     END SUBROUTINE partial_corr
!    ***********************************************************************************
!    ***********************************************************************************
     SUBROUTINE vmove(from, to, ifault)

     !     ALGORITHM AS274 APPL. STATIST. (1992) VOL.41, NO. 2

     !     Move variable from position FROM to position TO in an
     !     orthogonal reduction produced by AS75.1.
     !
     !--------------------------------------------------------------------------
     INTEGER, INTENT(IN)    :: from, to
     INTEGER, INTENT(OUT)   :: ifault

     !     Local variables

     double precision  :: d1, d2, x, d1new, d2new, cbar, sbar, y
     INTEGER    :: m, first, last, inc, m1, m2, mp1, col, pos, row

     !     Check input parameters

     ifault = 0
     IF (from < 1 .OR. from > ncol) ifault = ifault + 4
     IF (to < 1 .OR. to > ncol) ifault = ifault + 8
     IF (ifault /= 0) RETURN

     IF (from == to) RETURN

     IF (.NOT. rss_set) CALL ss()

     IF (from < to) THEN
       first = from
       last = to - 1
       inc = 1
     ELSE
       first = from - 1
       last = to
       inc = -1
     END IF

     DO m = first, last, inc

     !     Find addresses of first elements of R in rows M and (M+1).

       m1 = row_ptr(m)
       m2 = row_ptr(m+1)
       mp1 = m + 1
       d1 = d(m)
       d2 = d(mp1)

     !     Special cases.

       IF (d1 < vsmall .AND. d2 < vsmall) GO TO 40
       x = r(m1)
       IF (ABS(x) * SQRT(d1) < tol(mp1)) THEN
         x = zero
       END IF
       IF (d1 < vsmall .OR. ABS(x) < vsmall) THEN
         d(m) = d2
         d(mp1) = d1
         r(m1) = zero
         DO col = m+2, ncol
          m1 = m1 + 1
          x = r(m1)
          r(m1) = r(m2)
          r(m2) = x
          m2 = m2 + 1
         END DO ! col = m+2, ncol
         x = rhs(m)
         rhs(m) = rhs(mp1)
         rhs(mp1) = x
         GO TO 40
       ELSE IF (d2 < vsmall) THEN
         d(m) = d1 * x**2
         r(m1) = one / x
         r(m1+1:m1+ncol-m-1) = r(m1+1:m1+ncol-m-1) / x
         rhs(m) = rhs(m) / x
         GO TO 40
       END IF

     !     Planar rotation in regular case.

       d1new = d2 + d1*x**2
       cbar = d2 / d1new
       sbar = x * d1 / d1new
       d2new = d1 * cbar
       d(m) = d1new
       d(mp1) = d2new
       r(m1) = sbar
       DO col = m+2, ncol
         m1 = m1 + 1
         y = r(m1)
         r(m1) = cbar*r(m2) + sbar*y
         r(m2) = y - x*r(m2)
         m2 = m2 + 1
       END DO ! col = m+2, ncol
       y = rhs(m)
       rhs(m) = cbar*rhs(mp1) + sbar*y
       rhs(mp1) = y - x*rhs(mp1)

     !     Swap columns M and (M+1) down to row (M-1).
40       pos = m
       DO row = 1, m-1
         x = r(pos)
         r(pos) = r(pos-1)
         r(pos-1) = x
         pos = pos + ncol - row - 1
       END DO ! row = 1, m-1

     !     Adjust variable order (VORDER), the tolerances (TOL) and
     !     the vector of residual sums of squares (RSS).

       m1 = vorder(m)
       vorder(m) = vorder(mp1)
       vorder(mp1) = m1
       x = tol(m)
       tol(m) = tol(mp1)
       tol(mp1) = x
       rss(m) = rss(mp1) + d(mp1) * rhs(mp1)**2
     END DO

     RETURN
     END SUBROUTINE vmove
!    ***********************************************************************************
!    ***********************************************************************************
     SUBROUTINE reordr(list, n, pos1, ifault)

     !     ALGORITHM AS274  APPL. STATIST. (1992) VOL.41, NO. 2

     !     Re-order the variables in an orthogonal reduction produced by
     !     AS75.1 so that the N variables in LIST start at position POS1,
     !     though will not necessarily be in the same order as in LIST.
     !     Any variables in VORDER before position POS1 are not moved.

     !     Auxiliary routine called: VMOVE
     !
     !--------------------------------------------------------------------------
     INTEGER, INTENT(IN)               :: n, pos1
     INTEGER, DIMENSION(:), INTENT(IN) :: list
     INTEGER, INTENT(OUT)              :: ifault

     !     Local variables.

     INTEGER    :: next, i, l, j

     !     Check N.

     ifault = 0
     IF (n < 1 .OR. n > ncol+1-pos1) ifault = ifault + 4
     IF (ifault /= 0) RETURN

     !     Work through VORDER finding variables which are in LIST.

     next = pos1
     i = pos1
10     l = vorder(i)
     DO j = 1, n
       IF (l == list(j)) GO TO 40
     END DO
30     i = i + 1
     IF (i <= ncol) GO TO 10

     !     If this point is reached, one or more variables in LIST has not
     !     been found.

     ifault = 8
     RETURN

     !     Variable L is in LIST; move it up to position NEXT if it is not
     !     already there.

40     IF (i > next) CALL vmove(i, next, ifault)
     next = next + 1
     IF (next < n+pos1) GO TO 30

     RETURN
     END SUBROUTINE reordr
!    ***********************************************************************************
!    ***********************************************************************************
     SUBROUTINE hdiag(xrow, nreq, hii, ifault)

     !     ALGORITHM AS274  APPL. STATIST. (1992) VOL.41, NO. 2
     !
     !                         *1           *1
     ! The hat matrix H = x(X'X) x' = x(R'DR) x' = z'Dz

     !              *1
     ! where z = x'R

     ! Here we only calculate the diagonal element hii corresponding to one
     ! row (xrow).   The variance of the i-th least-squares residual is (1 - hii).
     !--------------------------------------------------------------------------
     INTEGER, INTENT(IN)                  :: nreq
     INTEGER, INTENT(OUT)                 :: ifault
     double precision, DIMENSION(:), INTENT(IN)  :: xrow
     double precision, INTENT(OUT)               :: hii

     !     Local variables

     INTEGER    :: col, row, pos
     double precision  :: total, wk(ncol)

     !     Some checks

     ifault = 0
     IF (nreq > ncol) ifault = ifault + 4
     IF (ifault /= 0) RETURN

     !     The elements of xrow.inv(R).sqrt(D) are calculated and stored in WK.

     hii = zero
     DO col = 1, nreq
       IF (SQRT(d(col)) <= tol(col)) THEN
         wk(col) = zero
       ELSE
         pos = col - 1
         total = xrow(col)
         DO row = 1, col-1
          total = total - wk(row)*r(pos)
          pos = pos + ncol - row - 1
         END DO ! row = 1, col-1
         wk(col) = total
         hii = hii + total**2 / d(col)
       END IF
     END DO ! col = 1, nreq

     RETURN
     END SUBROUTINE hdiag
!    ***********************************************************************************
!    ***********************************************************************************
     FUNCTION varprd(x, nreq) RESULT(fn_val)

     !     Calculate the variance of x'b where b consists of the first nreq
     !     least-squares regression coefficients.
     !
     !--------------------------------------------------------------------------
     INTEGER, INTENT(IN)                  :: nreq
     double precision, DIMENSION(:), INTENT(IN)  :: x
     double precision                            :: fn_val

     !     Local variables

     INTEGER    :: ifault, row
     double precision  :: var, wk(nreq)

     !     Check input parameter values

     fn_val = zero
     ifault = 0
     IF (nreq < 1 .OR. nreq > ncol) ifault = ifault + 4
     IF (nobs <= nreq) ifault = ifault + 8
     IF (ifault /= 0) THEN
!       WRITE(*, '(1x, a, i4)') 'Error in function VARPRD: 
!     *       ifault =', ifault
       RETURN
     END IF

     !     Calculate the residual variance estimate.

     var = sserr / (nobs - nreq)

     !     Variance of x'b = var.x'(inv R)(inv D)(inv R')x
     !     First call BKSUB2 to calculate (inv R')x by back-substitution.

     CALL BKSUB2(x, wk, nreq)
     DO row = 1, nreq
       IF(d(row) > tol(row)) fn_val = fn_val + wk(row)**2 / d(row)
     END DO

     fn_val = fn_val * var

     RETURN
     END FUNCTION varprd
!    ***********************************************************************************
!    ***********************************************************************************
     SUBROUTINE bksub2(x, b, nreq)

     !     Solve x = R'b for b given x, using only the first nreq rows and
     !     columns of R, and only the first nreq elements of R.
     !
     !--------------------------------------------------------------------------
     INTEGER, INTENT(IN)                  :: nreq
     double precision, DIMENSION(:), INTENT(IN)  :: x
     double precision, DIMENSION(:), INTENT(OUT) :: b

     !     Local variables

     INTEGER    :: pos, row, col
     double precision  :: temp

     !     Solve by back-substitution, starting from the top.

     DO row = 1, nreq
       pos = row - 1
       temp = x(row)
       DO col = 1, row-1
         temp = temp - r(pos)*b(col)
         pos = pos + ncol - col - 1
       END DO
       b(row) = temp
     END DO

     RETURN
     END SUBROUTINE bksub2
     END MODULE lsq
!    ***********************************************************************************
!    ***********************************************************************************
     subroutine WRegresion(X,Y,W,n,nvar,beta,sterr,se,r2,iopt)
     USE lsq
     IMPLICIT NONE
     INTEGER             :: i, ier, j, m, n,nvar,iopt
     double precision    :: x(n,nvar), y(n),W(n), xrow(0:nvar+1), &
          beta(0:nvar+1),var, covmat(231), sterr(0:nvar+1), &
          totalSS,se,r2
     LOGICAL             :: fit_const = .TRUE., lindep(0:20)

     ! Least-squares calculations
     m=nvar
     CALL startup(m, fit_const)
     DO i = 1, n
       xrow(0) = 1.0_dp
       DO j = 1, m
         xrow(j) = x(i,j)
       END DO
       CALL includ(W(i), xrow, y(i))
     END DO

     if (iopt.gt.0) then
          CALL sing(lindep, ier)
     end if

     ! Calculate progressive residual sums of squares
     CALL ss()
     var = rss(m+1) / (n - m - 1)

     ! Calculate least-squares regn. coeffs.
     CALL regcf(beta, m+1, ier)

     if (iopt.gt.0) then
          ! Calculate covariance matrix, and hence std. errors of coeffs.
          CALL cov(m+1, var, covmat, 231, sterr, ier)
          se=SQRT(var)
          totalSS = rss(1)
          r2=(totalSS - rss(m+1))/totalSS
     end if
     call endup()
     END
!    ***********************************************************************************
!    ***********************************************************************************
     subroutine qsortd(x,ind,n)
     ! Code converted using TO_F90 by Alan Miller
     ! Date: 2002-12-18  Time: 11:55:47
     implicit none
     integer, parameter  :: dp = SELECTED_REAL_KIND(12, 60)
     integer n,ind(n)
     double precision x(n)
     !***************************************************************************

     !                                                         ROBERT RENKA
     !                                                 OAK RIDGE NATL. LAB.

     ! THIS SUBROUTINE USES AN ORDER N*LOG(N) QUICK SORT TO SORT A double precision
     ! ARRAY X INTO INCREASING ORDER.  THE ALGORITHM IS AS FOLLOWS.  IND IS
     ! INITIALIZED TO THE ORDERED SEQUENCE OF INDICES 1,...,N, AND ALL INTERCHANGES
     ! ARE APPLIED TO IND.  X IS DIVIDED INTO TWO PORTIONS BY PICKING A CENTRAL
     ! ELEMENT T.  THE FIRST AND LAST ELEMENTS ARE COMPARED WITH T, AND
     ! INTERCHANGES ARE APPLIED AS NECESSARY SO THAT THE THREE VALUES ARE IN
     ! ASCENDING ORDER.  INTERCHANGES ARE THEN APPLIED SO THAT ALL ELEMENTS
     ! GREATER THAN T ARE IN THE UPPER PORTION OF THE ARRAY AND ALL ELEMENTS
     ! LESS THAN T ARE IN THE LOWER PORTION.  THE UPPER AND LOWER INDICES OF ONE
     ! OF THE PORTIONS ARE SAVED IN LOCAL ARRAYS, AND THE PROCESS IS REPEATED
     ! ITERATIVELY ON THE OTHER PORTION.  WHEN A PORTION IS COMPLETELY SORTED,
     ! THE PROCESS BEGINS AGAIN BY RETRIEVING THE INDICES BOUNDING ANOTHER
     ! UNSORTED PORTION.

     ! INPUT PARAMETERS -   N - LENGTH OF THE ARRAY X.

     !                      X - VECTOR OF LENGTH N TO BE SORTED.

     !                    IND - VECTOR OF LENGTH >= N.

     ! N AND X ARE NOT ALTERED BY THIS ROUTINE.

     ! OUTPUT PARAMETER - IND - SEQUENCE OF INDICES 1,...,N PERMUTED IN THE SAME
     !                          FASHION AS X WOULD BE.  THUS, THE ORDERING ON
     !                          X IS DEFINED BY Y(I) = X(IND(I)).

     !*********************************************************************

     ! NOTE -- IU AND IL MUST BE DIMENSIONED >= LOG(N) WHERE LOG HAS BASE 2.

     !*********************************************************************

     INTEGER   :: iu(21), il(21)
     INTEGER   :: m, i, j, k, l, ij, it, itt, indx
     double precision     :: r
     double precision :: t

     ! LOCAL PARAMETERS -

     ! IU,IL =  TEMPORARY STORAGE FOR THE UPPER AND LOWER
     !            INDICES OF PORTIONS OF THE ARRAY X
     ! M =      INDEX FOR IU AND IL
     ! I,J =    LOWER AND UPPER INDICES OF A PORTION OF X
     ! K,L =    INDICES IN THE RANGE I,...,J
     ! IJ =     RANDOMLY CHOSEN INDEX BETWEEN I AND J
     ! IT,ITT = TEMPORARY STORAGE FOR INTERCHANGES IN IND
     ! INDX =   TEMPORARY INDEX FOR X
     ! R =      PSEUDO RANDOM NUMBER FOR GENERATING IJ
     ! T =      CENTRAL ELEMENT OF X

     IF (n <= 0) RETURN

     ! INITIALIZE IND, M, I, J, AND R

     DO  i = 1, n
       ind(i) = i
     END DO
     m = 1
     i = 1
     j = n
     r = .375

     ! TOP OF LOOP

20   IF (i >= j) GO TO 70
     IF (r <= .5898437) THEN
       r = r + .0390625
     ELSE
       r = r - .21875
     END IF

     ! INITIALIZE K

30   k = i

     ! SELECT A CENTRAL ELEMENT OF X AND SAVE IT IN T

     ij = floor(i + r*(j-i))
     it = ind(ij)
     t = x(it)

     ! IF THE FIRST ELEMENT OF THE ARRAY IS GREATER THAN T,
     !   INTERCHANGE IT WITH T

     indx = ind(i)
     IF (x(indx) > t) THEN
       ind(ij) = indx
       ind(i) = it
       it = indx
       t = x(it)
     END IF

     ! INITIALIZE L

     l = j

     ! IF THE LAST ELEMENT OF THE ARRAY IS LESS THAN T,
     !   INTERCHANGE IT WITH T

     indx = ind(j)
     IF (x(indx) >= t) GO TO 50
     ind(ij) = indx
     ind(j) = it
     it = indx
     t = x(it)

     ! IF THE FIRST ELEMENT OF THE ARRAY IS GREATER THAN T,
     !   INTERCHANGE IT WITH T

     indx = ind(i)
     IF (x(indx) <= t) GO TO 50
     ind(ij) = indx
     ind(i) = it
     it = indx
     t = x(it)
     GO TO 50

     ! INTERCHANGE ELEMENTS K AND L

40   itt = ind(l)
     ind(l) = ind(k)
     ind(k) = itt

     ! FIND AN ELEMENT IN THE UPPER PART OF THE ARRAY WHICH IS
     !   NOT LARGER THAN T

50   l = l - 1
     indx = ind(l)
     IF (x(indx) > t) GO TO 50

     ! FIND AN ELEMENT IN THE LOWER PART OF THE ARRAY WHCIH IS NOT SMALLER THAN T

60   k = k + 1
     indx = ind(k)
     IF (x(indx) < t) GO TO 60

     ! IF K <= L, INTERCHANGE ELEMENTS K AND L

     IF (k <= l) GO TO 40

     ! SAVE THE UPPER AND LOWER SUBSCRIPTS OF THE PORTION OF THE
     !   ARRAY YET TO BE SORTED

     IF (l-i > j-k) THEN
       il(m) = i
       iu(m) = l
       i = k
       m = m + 1
       GO TO 80
     END IF

     il(m) = k
     iu(m) = j
     j = l
     m = m + 1
     GO TO 80

     ! BEGIN AGAIN ON ANOTHER UNSORTED PORTION OF THE ARRAY

70   m = m - 1
     IF (m == 0) RETURN
     i = il(m)
     j = iu(m)

80   IF (j-i >= 11) GO TO 30
     IF (i == 1) GO TO 20
     i = i - 1

     ! SORT ELEMENTS I+1,...,J.  NOTE THAT 1 <= I < J AND J-I < 11.

90   i = i + 1
     IF (i == j) GO TO 70
     indx = ind(i+1)
     t = x(indx)
     it = indx
     indx = ind(i)
     IF (x(indx) <= t) GO TO 90
     k = i

100  ind(k+1) = ind(k)
     k = k - 1
     indx = ind(k)
     IF (t < x(indx)) GO TO 100

     ind(k+1) = it
     GO TO 90
     end subroutine qsortd
!    ******************************************************************
!    ******************************************************************
     FUNCTION normal(z) RESULT(prob)
     IMPLICIT NONE

     INTEGER, PARAMETER    :: dp = SELECTED_REAL_KIND(12, 60)
     REAL (dp), INTENT(IN) :: z
     REAL (dp)             :: prob

     CALL nprob(z, p=prob)

     RETURN

     CONTAINS
!    ******************************************************************
!    ******************************************************************
     SUBROUTINE nprob(z, p, q, pdf)
     ! For a given number (z) of standard deviations from the mean, the
     ! probabilities to the left (p) and right (q) of z are calculated,
     ! and the probability density (pdf).

     ! Programmer - Alan J. Miller
     ! Latest revision - 7 August 1997
     ! This version is compatible with the ELF90 subset of Fortran 90.

     ! REFERENCE: ADAMS,A.G. AREAS UNDER THE NORMAL CURVE,
     ! ALGORITHM 39, COMPUTER J., VOL. 12, 197-8, 1969.

     IMPLICIT NONE

     INTEGER, PARAMETER :: doubleprec = SELECTED_REAL_KIND(15, 60)
     REAL(doubleprec), INTENT(IN)            :: z
     REAL(doubleprec), INTENT(OUT), OPTIONAL :: p, q, pdf

     ! Local variables
     REAL(doubleprec) :: a0 = 0.5D0, &
                              a1 = 0.398942280444D0, &
                              a2 = 0.399903438504D0, &
                              a3 = 5.75885480458D0, & 
                              a4 = 29.8213557808D0, &
                              a5 = 2.62433121679D0, &
                              a6 = 48.6959930692D0, &
                              a7 = 5.92885724438D0, &
                              b0 = 0.398942280385D0,&
                              b1 = 3.8052D-8, &
                              b2 = 1.00000615302D0, &
                              b3 = 3.98064794D-4,&
                              b4 = 1.98615381364D0, &
                              b5 = 0.151679116635D0,&
                              b6 = 5.29330324926D0, &
                              b7 = 4.8385912808D0, &
                              b8 = 15.1508972451D0,& 
                              b9 = 0.742380924027D0,&
                              b10 = 30.789933034D0, &
                              b11 = 3.99019417011D0, &
                              zero = 0.D0, one = 1.D0,& 
                              zabs, y, pp, qq, ppdf


     zabs = ABS(z)
     IF (zabs < 12.7D0) THEN
       y = a0*z*z
       ppdf = EXP(-y)*b0
       IF (PRESENT(pdf)) pdf = ppdf

     !     Z BETWEEN -1.28 AND +1.28

       IF (zabs < 1.28) THEN
         qq = a0 - zabs*(a1-a2*y/(y+a3-a4/(y+a5+a6/(y+a7))))
         IF(z < zero) THEN
          pp = qq
          qq = one - pp
         ELSE
          pp = one - qq
         END IF
         IF (PRESENT(p)) p = pp
         IF (PRESENT(q)) q = qq
         RETURN
       END IF

     !     ZABS BETWEEN 1.28 AND 12.7

       qq = ppdf/(zabs-b1+b2/(zabs+b3+b4/(zabs-b5+b6/(zabs+b7-b8/ &
                      (zabs+b9+b10/(zabs+b11))))))
       IF(z < zero) THEN
         pp = qq
         qq = one - pp
       ELSE
         pp = one - qq
       END IF
       IF (PRESENT(p)) p = pp
       IF (PRESENT(q)) q = qq
       RETURN

     !     Z FAR OUT IN TAIL

     ELSE
       ppdf = zero
       IF(z < zero) THEN
         pp = zero
         qq = one
       ELSE
         pp = one
         qq = zero
       END IF
       IF (PRESENT(p)) p = pp
       IF (PRESENT(q)) q = qq
       IF (PRESENT(pdf)) pdf = ppdf
       RETURN
     END IF

     RETURN
     END SUBROUTINE nprob

     END FUNCTION normal
!    ******************************************************************
!    ******************************************************************
     FUNCTION normdev(p) RESULT(fn_val)
     IMPLICIT NONE
     INTEGER, PARAMETER    :: dp = SELECTED_REAL_KIND(12, 60)
     REAL (dp), INTENT(IN) :: p
     REAL (dp)             :: fn_val

     ! Local variable
     INTEGER :: ifault

     CALL ppnd16(p, fn_val, ifault)
     IF (ifault /= 0) THEN
          !WRITE(*, *) 'Error in ppnd16: ifault =', ifault
     ENDIF

     RETURN

     CONTAINS
!    ******************************************************************
!    ******************************************************************
     SUBROUTINE ppnd16 (p, normal_dev, ifault)

     ! ALGORITHM AS241  APPL. STATIST. (1988) VOL. 37, NO. 3

     ! Produces the normal deviate Z corresponding to a given lower
     ! tail area of P; Z is accurate to about 1 part in 10**16.

     ! The hash sums below are the sums of the mantissas of the
     ! coefficients.   They are included for use in checking
     ! transcription.

     ! This ELF90-compatible version by Alan Miller - 20 August 1996
     ! N.B. The original algorithm is as a function; this is a subroutine

     IMPLICIT NONE

     INTEGER, PARAMETER      :: doub_prec = SELECTED_REAL_KIND(15, 60)
     REAL (doub_prec), INTENT(IN)  :: p
     INTEGER, INTENT(OUT)          :: ifault
     REAL (doub_prec), INTENT(OUT) :: normal_dev

     ! Local variables

     REAL (doub_prec) :: zero = 0.d0, one = 1.d0, half = 0.5d0, &
          split1 = 0.425d0, split2 = 5.d0, const1 = 0.180625d0, &
          const2 = 1.6d0, q, r

     ! Coefficients for P close to 0.5

     REAL (doub_prec) :: a0 = 3.3871328727963666080D0, &
                         a1 = 1.3314166789178437745D+2, &
                         a2 = 1.9715909503065514427D+3, &
                         a3 = 1.3731693765509461125D+4, &
                         a4 = 4.5921953931549871457D+4,&
                         a5 = 6.7265770927008700853D+4,&
                         a6 = 3.3430575583588128105D+4,&
                         a7 = 2.5090809287301226727D+3,&
                         b1 = 4.2313330701600911252D+1,&
                         b2 = 6.8718700749205790830D+2,&
                         b3 = 5.3941960214247511077D+3,&
                         b4 = 2.1213794301586595867D+4,&
                         b5 = 3.9307895800092710610D+4,&
                         b6 = 2.8729085735721942674D+4,&
                         b7 = 5.2264952788528545610D+3 
     ! HASH SUM AB           55.8831928806149014439

     ! Coefficients for P not close to 0, 0.5 or 1.

     REAL (doub_prec) :: c0 = 1.42343711074968357734D0, &
                         c1 = 4.63033784615654529590D0, &
                         c2 = 5.76949722146069140550D0, &
                         c3 = 3.64784832476320460504D0, &
                         c4 = 1.27045825245236838258D0, &
                         c5 = 2.41780725177450611770D-1, &
                         c6 = 2.27238449892691845833D-2, &
                         c7 = 7.74545014278341407640D-4, &
                         d1 = 2.05319162663775882187D0, &
                         d2 = 1.67638483018380384940D0, &
                         d3 = 6.89767334985100004550D-1, &
                         d4 = 1.48103976427480074590D-1, &
                         d5 = 1.51986665636164571966D-2, &
                         d6 = 5.47593808499534494600D-4, &
                         d7 = 1.05075007164441684324D-9
     ! HASH SUM CD           49.33206503301610289036

     ! Coefficients for P near 0 or 1.

     REAL (doub_prec) :: e0 = 6.65790464350110377720D0, &
                         e1 = 5.46378491116411436990D0, &
                         e2 = 1.78482653991729133580D0, &
                         e3 = 2.96560571828504891230D-1, &
                         e4 = 2.65321895265761230930D-2, &
                         e5 = 1.24266094738807843860D-3, &
                         e6 = 2.71155556874348757815D-5, & 
                         e7 = 2.01033439929228813265D-7, &
                         f1 = 5.99832206555887937690D-1, &
                         f2 = 1.36929880922735805310D-1, &
                         f3 = 1.48753612908506148525D-2, &
                         f4 = 7.86869131145613259100D-4, &
                         f5 = 1.84631831751005468180D-5, &
                         f6 = 1.42151175831644588870D-7, &
                         f7 = 2.04426310338993978564D-15 
     ! HASH SUM EF           47.52583317549289671629

     ifault = 0
     q = p - half
     IF (ABS(q) <= split1) THEN
       r = const1 - q * q
       normal_dev = q * (((((((a7*r + a6)*r + a5)*r + a4)*r + a3)*r &
          + a2)*r + a1)*r + a0) / &
           (((((((b7*r + b6)*r + b5)*r + b4)*r + b3)*r + b2)*r &
          + b1)*r + one)
       RETURN
     ELSE
       IF (q < zero) THEN
         r = p
       ELSE
         r = one - p
       END IF
       IF (r <= zero) THEN
         ifault = 1
         normal_dev = zero
         RETURN
       END IF
       r = SQRT(-LOG(r))
       IF (r <= split2) THEN
         r = r - const2
         normal_dev = (((((((c7*r + c6)*r + c5)*r + c4)*r + c3)*r &
               + c2)*r + c1)*r + c0) / & 
               (((((((d7*r + d6)*r + d5)*r + d4)*r + d3)*r + d2)*r &
               + d1)*r + one)
       ELSE
         r = r - split2
         normal_dev = (((((((e7*r + e6)*r + e5)*r + e4)*r + e3)*r &
               + e2)*r + e1)*r + e0) / &
               (((((((f7*r + f6)*r + f5)*r + f4)*r + f3)*r + f2)*r  &
               + f1)*r + one)
       END IF
       IF (q < zero) normal_dev = - normal_dev
       RETURN
     END IF

     RETURN
     END SUBROUTINE ppnd16

     END FUNCTION normdev
!    *******************************************************************************
!    *******************************************************************************
     subroutine DLLROCDirect ( &
          Z, X, W, Status, n,tagh, &      ! Input
          nt, &                    
          nvar, mode, &   
          nparm, IIm, hm, &               ! Input: Mean
          nparv, IIv, hv, &               ! Input: Variance
          nparr, IIr, hr, &               ! Input: ROC
          family, &                       ! Probit or logit
          Zb,nb, ntb, &                   ! Prediction
          p,kbin,coutcome,nboot, &        ! Fitting arguments
          seed, &                         ! Seed for the bootstrap
          cifit, level, &                 ! CI
          tpartial,npartial, &            ! Test
          Fp_, &                          ! Output
          coeff_, &                       ! Output
          ROC_, AUC_, &                   ! Output
          pvalue_)                        ! Output
     implicit none
     integer, parameter :: npart = 0, ntbc = 500
     integer n(3),nb,nt,ntb,p,kbin, nboot,tagh, Status(n(1)), &
          i,j,k,iboot,auxh, auxd, &
          nvar, mode(nvar), &
          nparv, nparm, nparr, &
          IIm(2,nparm), IIv(2,nparv), IIr(2,nparr), &
          IIt(npart), family, &
          n0b, n1b, &
          npartial, tpartial(npartial), seed
     logical coutcome, cifit
     double precision, external:: QQ
     ! Argumentos de entrada
     double precision Z(n(1),nvar), X(n(1)), W(n(1)), &
          Zb(nb,nvar), ROC_(ntb,nb), AUC_(nb,3), &
          hm(nparm), hv(nparv), hr(nparr), &
          Fp_(nb*ntb, nparr+npart+1, 3), &
          level, pvalue_(2, npartial), &
          coeff_(20)

     integer, allocatable::pm(:), pv(:)
     double precision, allocatable:: &
          AUC(:), Fp(:,:), &
          ROC(:,:),  &
          tb(:), t(:), &
          hcallr(:), &
          hcallrb(:), &
          hmb(:), hvb(:), &
          Z0(:,:),X0(:),W0(:), &
          Z1(:,:),X1(:),W1(:), &
          M0(:), V0(:), &
          X0b(:),X1b(:),Z0b(:,:),Z1b(:,:), &
          W0b(:), W1b(:), &
          M0b(:), V0b(:), &
          hopt(:,:), &
          AUCs(:,:), &
          Z0b2(:,:), &
          Z1b2(:,:), &
          Fps(:,:,:), &
          tbc(:),pvalue(:), &
          coeffb(:)
     allocate (AUC(nb), &
          ROC(ntb,nb), &
          Fp(nb*ntb, nparr+npart+1), &
          t(nt),tb(ntb), &
          pm(nparm),pv(nparv), &
          hcallr(nparr+npart+1), &
          hcallrb(nparr+npart+1), &
          hmb(nparm), hvb(nparv), &
          Z0(n(2),nvar),X0(n(2)),W0(n(2)), &
          Z1(n(3),nvar),X1(n(3)),W1(n(3)), &
          M0(n(2)), V0(n(2)), &
          M0b(n(1)), V0b(n(1)), &
          hopt(2,2), &
          AUCs(nboot,nb), & 
          Fps(nboot, nb*ntb, nparr+npart+1), & 
          tbc(ntbc), pvalue(2), coeffb(20))

     !FPF grid
     ! Estimation
     do i=1,nt
          t(i)=(i-1)*1.0/(nt-1)
     end do
     ! Prediction
     do i=1,ntb
          tb(i)=(i-1)*1.0/(ntb-1)
     end do
     ! Test
     do i=1,ntbc
          tbc(i)=(i-1)*1.0/(ntbc-1)
     end do

     !Load healthy and diseased individuals
     auxh = 0
     auxd = 0
     do i = 1, n(1)
          if (Status(i).eq.tagh) then
               auxh = auxh+1
               W0(auxh) = W(i)
               X0(auxh) = X(i)
               do j=1,nvar
                    Z0(auxh,j)=Z(i,j)
               end do
          else
               auxd = auxd+1
               W1(auxd) = W(i)
               X1(auxd) = X(i)
               do j = 1,nvar
                    Z1(auxd,j)=Z(i,j)
               end do
          end if
     end do

     ! Polynomial order: mean and variace
     pm = p
     pv = p
     ! Bandwithds
     do i=1,nparr
          hcallr(i)= hr(i)
     end do
     ! FPR
     hcallr(nparr+1)=-1
     ! Original bandwithds
     hmb = hm
     hvb = hv
     hcallrb = hcallr
     ! Estimation
     call ROCDirecta(Z0,X0,W0,n(2),Z1,X1,W1,n(3), &
          t, nt, &
          family, &
          nvar, mode, & 
          nparm, IIm, & 
          nparv, IIv, &
          nparr, IIr, & 
          npart, IIt, &
          kbin, pm, pv, &
          hm, hv, hcallr, 0, &
          Zb,nb,tb,ntb, &
          M0, V0, ROC_, AUC, Fp, coeff_)
     do i=1,nb
          AUC_(i,1)= AUC(i)
          do j=1,ntb
               do k=1,nparr+npart+1
                    Fp_((i-1)*ntb+j,k,1)=Fp((i-1)*ntb+j,k)
               end do     
          end do     
     end do
     !Set seed
     call init_random_seed(seed)
     !Confidence intervals
     if(cifit) then
          allocate(X0b(n(1)),X1b(n(1)), &
          Z0b(n(1),nvar), &  
          Z1b(n(1),nvar), &
          W0b(n(1)),W1b(n(1)))
          
          do iboot = 1,nboot
               call sampleROC(Z0, X0, W0, n(2), Z1, X1, W1, n(3), nvar, &
                    coutcome, Z0b, X0b, W0b, n0b, Z1b, X1b, W1b, n1b)
               hm = hmb
               hv = hvb
               hcallr = hcallrb
               allocate (Z0b2(n0b,nvar), Z1b2(n1b,nvar))
               do i=1,n0b
                    do j=1,nvar
                         Z0b2(i,j) = Z0b(i,j)
                    end do
               end do
               do i=1,n1b
                    do j=1,nvar
                         Z1b2(i,j) = Z1b(i,j)
                    end do
               end do
               call ROCDirecta(Z0b2,X0b,W0b,n0b, &
                    Z1b2,X1b,W1b,n1b, &
                    t, nt, &
                    family, &
                    nvar, mode, & 
                    nparm, IIm, & 
                    nparv, IIv, &
                    nparr, IIr, & 
                    npart, IIt, &
                    kbin, pm, pv, &
                    hm, hv, hcallr,0, &
                    Zb,nb,tb,ntb, &
                    M0b, V0b, ROC, AUC, & 
                    Fp, coeffb)
               deallocate (Z0b2, Z1b2)
               do i=1,nb
                    AUCs(iboot,i) = AUC(i)
                    do j=1,ntb
                         do k=1,nparr+npart+1
                              Fps(iboot,(i-1)*ntb+j,k)=Fp((i-1)*ntb+j,k)
                         end do
                    end do
               end do
          end do
          do i=1,nb
               ! AUC
               AUC_(i,2) = QQ(AUCs(1,i), nboot, (1-level)/2)
               AUC_(i,3) = QQ(AUCs(1,i), nboot, 1-((1-level)/2))
               do j=1,ntb
                    do k=1,nparr+npart+1
                              Fp_((i-1)*ntb+j,k,2)=QQ(Fps(1,(i-1)*ntb+j,k), &
                              nboot, (1-level)/2)
                         Fp_((i-1)*ntb+j,k,3)=QQ(Fps(1,(i-1)*ntb+j,k), &
                              nboot, 1-((1-level)/2))
                    end do     
               end do
          end do
          deallocate (X0b, X1b, Z0b, Z1b, W0b, W1b)
     endif
     ! Covariate effect
     if(npartial.gt.0) then
          do i=1,npartial
               hm = hmb
               hv = hvb
               hcallr = hcallrb
               pvalue = 0     
               call ROCDirectab(Z0,X0,W0,n(2), &
                    Z1,X1,W1,n(3), &
                    t,nt, &
                    family, &
                    nvar, mode, & 
                    nparm, IIm, & 
                    nparv, IIv, &
                    nparr, IIr, & 
                    npart, IIt, &
                    kbin, pm, pv, &
                    hm, hv, hcallr,0, &
                    tpartial(i), &
                    tbc,ntbc, &
                    pvalue)
               pvalue_(1,i) = pvalue(1) 
               pvalue_(2,i) = pvalue(2)
          end do
     end if
     deallocate(AUC, ROC, Fp, &
          t,tb, &
          pm,pv, &
          hcallr, &
          hcallrb, &   
          hmb, hvb, &
          Z0,X0,W0, &
          Z1,X1,W1, &
          M0, V0, &
          M0b, V0b, &
          hopt, &
          AUCs, &
          Fps,tbc,pvalue,coeffb)   
     end
!    ******************************************************************************
!    ******************************************************************************
     subroutine DLLROCInduced( &
          Z,X,W, Status, n,tagh, &                    ! Input
          Zb, nb, nt, &                               ! Prediction
          p, h, kbin, coutcome, nboot, seed, &        ! Fit args
          cifit, level, test, accuracy, &             ! Fit args
          M0_, V0_, M1_, V1_,     &                   ! Output
          ROC_, AUC_, AROC_, YI_, TH_, pvalue_)       ! Output
     
     implicit none
     integer n(3),nb,nt,p,kbin, nboot,tagh, Status(n(1)), &
          i,j,boot,ir0(n(2)), ir1(n(3)),auxh, auxd, &
          ptemp(2,2), seed, n0b, n1b
     logical coutcome, cifit, test, accuracy(4)
     double precision, external:: QQ
     double precision Z(n(1)), X(n(1)), W(n(1)), &
          Zb(nb), h(2,2), hb(2,2), &
          ROC_(nt,nb), AUC_(nb,3), AROC_(nt, 3), YI_(nb,3), TH_(nb,3), &
          M0_(nb,3),V0_(nb,3),M1_(nb,3),V1_(nb,3), &  
          pvalue_, level
     double precision, allocatable:: &
          t(:), & 
          Z0(:),X0(:),W0(:), &
          Z1(:),X1(:),W1(:), &
          M0(:), V0(:), &
          M1(:), V1(:), &
          M0b(:), V0b(:), &
          M1b(:), V1b(:), &
          M0p(:),V0p(:),M1p(:),V1p(:), &
          ROC(:,:),AUC(:),AROC(:), &
          YI(:),TH(:), &
          M1c(:), V1c(:), &     
          pvalue(:), &          
          X0b(:),X1b(:),Z0b(:,:),Z1b(:,:), &
          W0b(:), W1b(:), &
          Err0(:),Err1(:), &
          Err0b(:),Err1b(:), &
          Z0b2(:), Z1b2(:), &
          M1cb(:), V1cb(:), &
          hopt(:,:), &
          M0s(:,:), V0s(:,:), &
          M1s(:,:), V1s(:,:), &
          AUCs(:,:), AROCs(:,:), &
          YIs(:,:), THs(:,:)
     allocate (t(nt),  &
          Z0(n(2)),X0(n(2)),W0(n(2)), &
          Z1(n(3)),X1(n(3)),W1(n(3)), &
          M0(n(2)), V0(n(2)), &
          M1(n(3)), V1(n(3)), &
          M0p(nb),V0p(nb),M1p(nb),V1p(nb), &
          ROC(nt,nb),AUC(nb),AROC(nt), &
          YI(nb),TH(nb), &
          M1c(n(3)), V1c(n(3)),  &
          pvalue(20), &
          hopt(2,2))

     !FPF
     do i=1,nt
          t(i)=(i-1)*1.0/(nt-1)
     end do
     !Load healthy and diseased observations
     auxh=0
     auxd=0
     do i=1,n(1)
          if (Status(i).eq.tagh) then
               auxh=auxh+1
               W0(auxh)=W(i)
               Z0(auxh)=Z(i)
               X0(auxh)=X(i)
          else
               auxd=auxd+1
               W1(auxd)=W(i)
               Z1(auxd)=Z(i)
               X1(auxd)=X(i)
          end if          
     end do
     !Original bandwitdhs
     hb = h          
     !Estimation
     ptemp(1,1)=p
     ptemp(2,1)=0
     ptemp(1,2)=p
     ptemp(2,2)=0
     call  ROCInduced(Z0,X0,W0,n(2),Z1,X1,W1,n(3), &
          kbin,ptemp,h, &
          Zb,nb,t,nt, &
          M0,V0,M1,V1, &
          M0p,V0p,M1p,V1p, &
          ROC,AUC,AROC, &
          accuracy(1),accuracy(2),accuracy(3:4),YI,TH, &
          .FALSE.,M1c,V1c)
     do i=1,nb
          AUC_(i,1)= AUC(i)
          YI_(i,1) = YI(i)
          TH_(i,1) = TH(i)
          M0_(i,1) = M0p(i)
          V0_(i,1) = V0p(i)
          M1_(i,1) = M1p(i)
          V1_(i,1) = V1p(i)
          do j=1,nt
               AROC_(j,1)=AROC(j)
               ROC_(j,i) = ROC(j,i)
          end do
     end do
     !Set seed
     call init_random_seed(seed)
     !Inference about the effect
     if (test) then
          h = hb
          call ROCInducedb(Z0,X0,W0,n(2), &
               Z1,X1,W1,n(3), &
               kbin, ptemp, h,t,nt,pvalue)
          pvalue_=pvalue(1)
     end if
     !Confidence intervals
     if (cifit) then
          allocate(X0b(n(1)),X1b(n(1)), &
          Z0b(n(1),1), &  
          Z1b(n(1),1), &
          W0b(n(1)),W1b(n(1)), &
          Err0(n(2)),Err1(n(3)),& 
          Err0b(n(2)),Err1b(n(3)), &
          M0s(nboot,nb), V0s(nboot,nb), &
          M1s(nboot,nb), V1s(nboot,nb), &
          AUCs(nboot,nb), AROCs(nboot,nt), &
          YIs(nboot, nb), THs(nboot,nb))
          
          ! If conditional on the true disease status
          Err0=(X0-M0)/sqrt(V0)
          Err1=(X1-M1)/sqrt(V1)
          do boot = 1, nboot
               if(coutcome) then
                   allocate (Z0b2(n(2)), Z1b2(n(3)), &
                             M1cb(n(3)), V1cb(n(3)), &  
                             M0b(n(2)), V0b(n(2)), &
                             M1b(n(3)), V1b(n(3)))
                   call sample_int(n(2),n(2),ir0)     
                   do i=1,n(2)
                        X0b(i) = M0(i) + sqrt(V0(i))*Err0(ir0(i))
                        Z0b2(i) = Z0(i)
                        W0b(i) = W0(i)                    
                   end do
                   n0b = n(2)
                   call sample_int(n(3),n(3),ir1)
                   do i=1,n(3)
                        X1b(i) = M1(i) + sqrt(V1(i))*Err1(ir1(i))
                        Z1b2(i) = Z1(i)
                        W1b(i) = W1(i)  
                   end do
                   n1b = n(3)
               else
                   call sampleROC(Z0, X0, W0, n(2), &
                                  Z1, X1, W1, n(3),&
                                   1, .FALSE., &
                                  Z0b, X0b, W0b, n0b, &
                                  Z1b, X1b, W1b, n1b)
                   allocate (Z0b2(n0b), Z1b2(n1b), &
                             M1cb(n1b), V1cb(n1b), &
                             M0b(n0b), V0b(n0b), &
                             M1b(n1b), V1b(n1b))                             
                   do i=1,n0b
                        Z0b2(i) = Z0b(i,1)
                   end do
                   do i=1,n1b
                        Z1b2(i) = Z1b(i,1)
                   end do
               end if
               h = hb
               call  ROCInduced(Z0b2,X0b,W0b,n0b, &
                    Z1b2,X1b,W1b,n1b, &
                    kbin,ptemp,h, &
                    Zb,nb,t,nt, &
                    M0b,V0b,M1b,V1b, &
                    M0p,V0p,M1p,V1p, &
                    ROC,AUC,AROC, &
                    accuracy(1),accuracy(2),accuracy(3:4),YI,TH, &
                    .FALSE.,M1cb,V1cb)
               deallocate (Z0b2, Z1b2, M1cb, V1cb, M0b, V0b, M1b, V1b)
               ! Store the results
                    do i=1,nb
                    AUCs(boot,i)= AUC(i)
                    YIs(boot,i) = YI(i)
                    THs(boot,i) = TH(i)
                    M0s(boot,i) = M0p(i)
                    V0s(boot,i) = V0p(i)
                    M1s(boot,i) = M1p(i)
                    V1s(boot,i) = V1p(i)
               end do
               do i=1,nt
                    AROCs(boot,i)=AROC(i)
               end do
          end do
          ! Obtain the quantiles
          do i=1,nb
               ! AUC
               AUC_(i,2) = QQ(AUCs(1,i), nboot, (1-level)/2)
               AUC_(i,3) = QQ(AUCs(1,i), nboot, 1-((1-level)/2))
               !TH
               TH_(i,2) = QQ(THs(1,i), nboot, (1-level)/2)
               TH_(i,3) = QQ(THs(1,i), nboot, 1-((1-level)/2))
               !YI
               YI_(i,2) = QQ(YIs(1,i), nboot, (1-level)/2)
               YI_(i,3) = QQ(YIs(1,i), nboot, 1-((1-level)/2))
               ! Healthy
               M0_(i,2) = QQ(M0s(1,i), nboot, (1-level)/2)
               M0_(i,3) = QQ(M0s(1,i), nboot, 1-((1-level)/2))
               V0_(i,2) = QQ(V0s(1,i), nboot, (1-level)/2)
               V0_(i,3) = QQ(V0s(1,i), nboot, 1-((1-level)/2))
               ! Diseased                                        
               M1_(i,2) = QQ(M1s(1,i), nboot, (1-level)/2)
               M1_(i,3) = QQ(M1s(1,i), nboot, 1-((1-level)/2))
               V1_(i,2) = QQ(V1s(1,i), nboot, (1-level)/2)
               V1_(i,3) = QQ(V1s(1,i), nboot, 1-((1-level)/2))
          end do
          do i=1,nt
               AROC_(i,2) = QQ(AROCs(1,i), nboot, (1-level)/2)
               AROC_(i,3) = QQ(AROCs(1,i), nboot, 1-((1-level)/2))
          end do
          deallocate (X0b, X1b, &
                      Z0b, Z1b, &
                      W0b, W1b, &
                      Err0, Err1, &
                      Err0b, Err1b,&          
                      M0s,V0s, &
                      M1s,V1s, &
                      AUCs,AROCs, &
                      YIs,THs)    
     end if
     deallocate (t,  &
          Z0,X0,W0, &
          Z1,X1,W1,  &
          M0, V0, &
          M1, V1, &
          M0p,V0p,M1p,V1p, &
          ROC,AUC,AROC, &
          YI,TH, &
          M1c, V1c, &
          pvalue, &
          hopt)
     end
!    ******************************************************************
!    ******************************************************************     
     SUBROUTINE spline(x,y,n,y2)
     INTEGER n,NMAX
     double precision x(n),y(n),y2(n)
     PARAMETER (NMAX=500)
     INTEGER i,k
     double precision p,qn,sig,un,u(NMAX)
     
     y2(1)=0.
     u(1)=0.
     do i=2,n-1
          sig=(x(i)-x(i-1))/(x(i+1)-x(i-1))
          p=sig*y2(i-1)+2.
          y2(i)=(sig-1.)/p
          u(i)=(6.*((y(i+1)-y(i))/(x(i+1)-x(i))-(y(i)-y(i-1)) &
               /(x(i)-x(i-1)))/(x(i+1)-x(i-1))-sig*u(i-1))/p
     enddo

     qn=0.
     un=0.
     y2(n)=(un-qn*u(n-1))/(qn*y2(n-1)+1.)
     do k= n-1,1,-1
          y2(k)=y2(k)*y2(k+1)+u(k)
     enddo
     return
     END
!    ******************************************************************************
!    ******************************************************************************     
     subroutine splint(xa,ya,y2a,n,x,y)
     integer n,k,khi,klo
     double precision x,y,xa(n),y2a(n),ya(n)
     double precision a,b,h
     klo=1
     khi=n
1     if (khi-klo.gt.1) then
          k=(khi+klo)/2
          if(xa(k).gt.x)then
               khi=k
          else
               klo=k
          endif
          goto 1
     endif
     h = xa(khi) - xa(klo)
     a=(xa(khi)-x)/h
     b=(x-xa(klo))/h
     y=a*ya(klo)+b*ya(khi)+ &
          ((a**3-a)*y2a(klo)+(b**3-b)*y2a(khi))*(h**2)/6.
     return
     end