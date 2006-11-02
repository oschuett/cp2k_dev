SUBROUTINE hcth(iparset,rho,drho,exc,vxc,vxcg)

! Purpose: Calculate the gradient-corrected exchange energy and potential
!          of Hamprecht, Cohen, Tozer, and Handy (HCTH) for a closed shell
!          density.

! Literature: - F. A. Hamprecht, A. J. Cohen, D. J. Tozer, and N. C. Handy,
!               J. Chem. Phys. 109, 6264 (1998) -> HCTH/93
!             - A. D. Boese, N. L. Doltsinis, N. C. Handy, and M. Sprik,
!               J. Chem. Phys. 112, 1670 (2000) -> HCTH/120 and HCTH/147
!             - A. D. Boese and N. C. Handy,
!               J. Chem. Phys. 114, 5497 (2001) -> HCTH/407
!             - J. P. Perdew and Y. Wang,
!               Phys. Rev. B 45, 13244 (1992) -> PW92

! History: - Creation (27.07.2001, Matthias Krack)

! ***************************************************************************

  USE basic_data_types, ONLY: dp

  IMPLICIT NONE

  INTEGER, INTENT(IN)   :: iparset
  REAL(dp), INTENT(IN)  :: drho,rho
  REAL(dp), INTENT(OUT) :: exc,vxc,vxcg

! *** Local parameters ***

  REAL(dp), PARAMETER :: eps_rho = 1.0E-12_dp,&
                         f13 = 1.0_dp/3.0_dp,&
                         f83 = 8.0_dp*f13,&
                         cx_vwn_e = -0.738558766382022405884230032681_dp,&
                         cx_vwn_v = 4.0_dp*f13*cx_vwn_e,&
                         pi = 3.14159265358979323846264338328_dp,&
                         rsfac = 0.620350490899400016668006812048_dp,&
                         two13 = 1.25992104989487316476721060728_dp

! *** LSDA correlation parametrisation (PW92) ***

  REAL(dp), PARAMETER :: a0 = 0.031091_dp,&
                         a1 = 0.015545_dp,&
                         alpha0 = 0.21370_dp,&
                         alpha1 = 0.20548_dp

  REAL(dp), DIMENSION(4), PARAMETER :: beta0 = (/ 7.59570_dp,&
                                                  3.58760_dp,&
                                                  1.63820_dp,&
                                                  0.49294_dp/),&
                                       beta1 = (/14.11890_dp,&
                                                  6.19770_dp,&
                                                  3.36620_dp,&
                                                  0.62517_dp/)

! *** GGA parametrisation (HCTH/iparset) ***

  REAL(dp), PARAMETER :: gamma_cab = 0.006_dp,&
                         gamma_css = 0.200_dp,&
                         gamma_xss = 0.004_dp

! *** Local variables ***

  REAL(dp) :: dgcabddrho,dgcabdrho,dgcabds,dgcssddrho,dgcssdrho,dgcssds,&
              dgdrs,dgxssddrho,dgxssdrho,dgxssds,drsdrho,drhos,ecab,ecss,&
              exss,g,gcab,gcss,gs2,gxss,p,q,rho13,rho43,rhos13,rhos43,rhos,&
              rs,rs12,s,s2,u,vcab,vcss,vxss,x,y

  REAL(dp), DIMENSION(0:4) :: ccab,ccss,cxss

! ---------------------------------------------------------------------------

  IF (rho < eps_rho) THEN

    exc = 0.0_dp
    vxc = 0.0_dp
    vxcg = 0.0_dp

  ELSE

!   *** Load the HCTH parameter set HCTH/iparset ***

    SELECT CASE (iparset)
    CASE (93)
      cxss(0) =  0.109320E+01_dp
      ccss(0) =  0.222601E+00_dp
      ccab(0) =  0.729974E+00_dp
      cxss(1) = -0.744056E+00_dp
      ccss(1) = -0.338622E-01_dp
      ccab(1) =  0.335287E+01_dp
      cxss(2) =  0.559920E+01_dp
      ccss(2) = -0.125170E-01_dp
      ccab(2) = -0.115430E+02_dp
      cxss(3) = -0.678549E+01_dp
      ccss(3) = -0.802496E+00_dp
      ccab(3) =  0.808564E+01_dp
      cxss(4) =  0.449357E+01_dp
      ccss(4) =  0.155396E+01_dp
      ccab(4) = -0.447857E+01_dp
    CASE (120)
      cxss(0) =  0.109163E+01_dp
      ccss(0) =  0.489508E+00_dp
      ccab(0) =  0.514730E+00_dp
      cxss(1) = -0.747215E+00_dp
      ccss(1) = -0.260699E+00_dp
      ccab(1) =  0.692982E+01_dp
      cxss(2) =  0.507833E+01_dp
      ccss(2) =  0.432917E+00_dp
      ccab(2) = -0.247073E+02_dp
      cxss(3) = -0.410746E+01_dp
      ccss(3) = -0.199247E+01_dp
      ccab(3) =  0.231098E+02_dp
      cxss(4) =  0.117173E+01_dp
      ccss(4) =  0.248531E+01_dp
      ccab(4) = -0.113234E+02_dp
    CASE (147)
      cxss(0) =  0.109025E+01_dp
      ccss(0) =  0.562576E+00_dp
      ccab(0) =  0.542352E+00_dp
      cxss(1) = -0.799194E+00_dp
      ccss(1) =  0.171436E-01_dp
      ccab(1) =  0.701464E+01_dp
      cxss(2) =  0.557212E+01_dp
      ccss(2) = -0.130636E+01_dp
      ccab(2) = -0.283822E+02_dp
      cxss(3) = -0.586760E+01_dp
      ccss(3) =  0.105747E+01_dp
      ccab(3) =  0.350329E+02_dp
      cxss(4) =  0.304544E+01_dp
      ccss(4) =  0.885429E+00_dp
      ccab(4) = -0.204284E+02_dp
    CASE (407)
      cxss(0) =  0.108184E+01_dp
      ccss(0) =  0.118777E+01_dp
      ccab(0) =  0.589076E+00_dp
      cxss(1) = -0.518339E+00_dp
      ccss(1) = -0.240292E+01_dp
      ccab(1) =  0.442374E+01_dp
      cxss(2) =  0.342562E+01_dp
      ccss(2) =  0.561741E+01_dp
      ccab(2) = -0.192218E+02_dp
      cxss(3) = -0.262901E+01_dp
      ccss(3) = -0.917923E+01_dp
      ccab(3) =  0.425721E+02_dp
      cxss(4) =  0.228855E+01_dp
      ccss(4) =  0.624798E+01_dp
      ccab(4) = -0.420052E+02_dp
    CASE DEFAULT
      STOP "in SUBROUTINE hcth: Invalid HCTH parameter set requested"
    END SELECT

!   *** rho_sigma = rho/2 = rho_alpha = rho_beta (same for |nabla rho|) ***

    rhos = 0.5_dp*rho
    drhos = 0.5_dp*drho

    rhos13 = rhos**f13
    rhos43 = rhos13*rhos

    rho13 = two13*rhos13
    rho43 = rho13*rho

!   *** LSDA exchange part (VWN) ***

    exss = cx_vwn_e*rho43
    vxss = cx_vwn_v*rho13

!   *** LSDA correlation part (PW92) ***

!   *** G(rho_sigma,0) => spin polarisation zeta = 1 ***

    rs = rsfac/rhos13
    rs12 = SQRT(rs)
    q = 2.0_dp*a1*(beta1(1) + (beta1(2) + (beta1(3) +&
                   beta1(4)*rs12)*rs12)*rs12)*rs12
    p = 1.0_dp + 1.0_dp/q
    x = -2.0_dp*a1*(1.0_dp + alpha1*rs)
    y = LOG(p)
    g = x*y
    dgdrs = -2.0_dp*a1*alpha1*y -&
             x*a1*(beta1(1)/rs12 + 2.0_dp*beta1(2) +&
                   3.0_dp*beta1(3)*rs12 + 4.0_dp*beta1(4)*rs)/(p*q*q)
    drsdrho = -f13*rs/rho
    ecss = rho*g
    vcss = g + rho*dgdrs*drsdrho

!   *** G(rho_alpha,rho_beta) => spin polarisation zeta = 0 ***

    rs = rsfac/rho13
    rs12 = SQRT(rs)
    q = 2.0_dp*a0*(beta0(1) + (beta0(2) + (beta0(3) +&
                   beta0(4)*rs12)*rs12)*rs12)*rs12
    p = 1.0_dp + 1.0_dp/q
    x = -2.0_dp*a0*(1.0_dp + alpha0*rs)
    y = LOG(p)
    g = x*y
    dgdrs = -2.0_dp*a0*alpha0*y -&
             x*a0*(beta0(1)/rs12 + 2.0_dp*beta0(2) +&
                   3.0_dp*beta0(3)*rs12 + 4.0_dp*beta0(4)*rs)/(p*q*q)
    drsdrho = -f13*rs/rho
    ecab = rho*g - ecss
    vcab = g + rho*dgdrs*drsdrho - vcss

!   *** GGA part (HCTH) ***

    s = drhos/rhos43
    s2 = s*s
    x = -f83/rho
    y = 2.0_dp/(drho*drho)

!   *** g_x(rho_sigma,rho_sigma) ***

    gs2 = gamma_xss*s2
    q = 1.0_dp/(1.0_dp + gs2)
    u = gs2*q
    gxss = cxss(0) + (cxss(1) + (cxss(2) + (cxss(3) + cxss(4)*u)*u)*u)*u
    dgxssds = q*(cxss(1) + (2.0_dp*cxss(2) + (3.0_dp*cxss(3) +&
                 4.0_dp*cxss(4)*u)*u)*u)*u
    dgxssdrho = x*dgxssds
    dgxssddrho = y*dgxssds

!   *** g_c(rho_sigma,rho_sigma) ***

    gs2 = gamma_css*s2
    q = 1.0_dp/(1.0_dp + gs2)
    u = gs2*q
    gcss = ccss(0) + (ccss(1) + (ccss(2) + (ccss(3) + ccss(4)*u)*u)*u)*u
    dgcssds = q*(ccss(1) + (2.0_dp*ccss(2) + (3.0_dp*ccss(3) +&
                 4.0_dp*ccss(4)*u)*u)*u)*u
    dgcssdrho = x*dgcssds
    dgcssddrho = y*dgcssds

!   *** g_c(rho_alpha,rho_beta) ***

    gs2 = gamma_cab*s2
    q = 1.0_dp/(1.0_dp + gs2)
    u = gs2*q
    gcab = ccab(0) + (ccab(1) + (ccab(2) + (ccab(3) + ccab(4)*u)*u)*u)*u
    dgcabds = q*(ccab(1) + (2.0_dp*ccab(2) + (3.0_dp*ccab(3) +&
                 4.0_dp*ccab(4)*u)*u)*u)*u
    dgcabdrho = x*dgcabds
    dgcabddrho = y*dgcabds

!   *** Finally collect all contributions ***

    exc = exss*gxss + ecss*gcss + ecab*gcab
    vxc = vxss*gxss + exss*dgxssdrho +&
          vcss*gcss + ecss*dgcssdrho +&
          vcab*gcab + ecab*dgcabdrho
    vxcg = exss*dgxssddrho + ecss*dgcssddrho + ecab*dgcabddrho

  END IF

END SUBROUTINE hcth
