
!****************************************************************
!*                                                              *
!*   module MO__SALSA_DYNAMICS                               *
!*                                                              *
!*   Contains subroutines and functions that are used           *
!*   to calculate aerosol dynamics                              *
!*                                                              *
!****************************************************************
!----------------------------------------------------------------------------
!!$   Copyright 2014 Atmospheric Research Centre of Eastern Finland,
!!$         Finnish Meteorological Institute, Kuopio, Finland
!!$
!!$   Licensed under the Apache License, Version 2.0 (the "License");
!!$   you may not use this file except in compliance with the License.
!!$   You may obtain a copy of the License at
!!$
!!$       http://www.apache.org/licenses/LICENSE-2.0
!!$
!!$   Unless required by applicable law or agreed to in writing, software
!!$   distributed under the License is distributed on an "AS IS" BASIS,
!!$   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
!!$   See the License for the specific language governing permissions and
!!$   limitations under the License.
!----------------------------------------------------------------------------
MODULE mo_salsa_dynamics


CONTAINS

  ! this calculated for empty bins too!!!
  ! fxm: test well, esp. self-coagulation (but other bits too!)
  ! AL_note: Diagnostic variables of cond and nucl mass
  !********************************************************************
  !
  ! subroutine COAGULATION(kproma,kbdim,klev, &
  !       pnaero,pvols,pdwet, &
  !       pcore, ptstep)
  !
  !********************************************************************
  !
  ! Purpose:
  ! --------
  ! Calculates particle loss and change in size distribution
  !  due to (Brownian) coagulation
  !
  !
  ! Method:
  ! -------  
  ! Semi-implicit, non-iterative method:
  !  Volume concentrations of the smaller colliding particles
  !  added to the bin of the larger colliding particles.
  !  Start from first bin and use the updated number and volume
  !  for calculation of following bins. NB! Our bin numbering
  !  does not follow particle size in regime 2.
  !
  !Schematic for bin numbers in different regimes:
  !        	 1             			2             
  !    +-------------------------------------------+
  !  a | 1 | 2 | 3 || 4 | 5 | 6 | 7 |  8 |  9 | 10||
  !  b |           ||11 |12 |13 |14 | 15 | 16 | 17||
  !    +-------------------------------------------+
  !
  ! Exact coagulation coefficients for each pressure level
  !  are calculated in subroutine SET_COAGC (in mo_salsa_init) 
  !  which is called once at the beginning of the simulation 
  !  from model driver. In subroutine COAGULATION, these exact 
  !  coefficients are scaled according to current particle wet size
  !  (linear scaling).
  !  
  ! Juha: Now also considers coagulation between hydrometeors,
  !       and hydrometeors and aerosols.
  !
  !       Since the bins are organized in terms of the dry size of
  !       of the condensation nucleus, while coagulation kernell is
  !       calculated with the actual hydrometeor size, some assumptions
  !       are laid out:
  !                 1. Cloud droplets from each size bin are lost by 
  !                    coagulation with other cloud droplets that have 
  !                    larger condensation nucleus.
  !
  !                 2. Cloud droplets from each size bin are lost by 
  !                    coagulation with all drizzle bins, regardless of
  !                    the nucleus size in the latter (collection of cloud
  !                    droplets by rain).
  !
  !                 3. Coagulation between drizzle bins acts like 1.
  !
  !       ISSUES:
  !           Process selection should be made smarter - now just lots of IFs 
  !           inside loops. Bad.
  !
  !
  ! Interface:
  ! ----------
  ! Called from main aerosol model
  ! 
  !
  ! Coded by:
  ! ---------
  ! Hannele Korhonen (FMI) 2005 
  ! Harri Kokkola (FMI) 2006
  ! Tommi Bergman (FMI) 2012
  ! Matti Niskanen(FMI) 2012
  ! Anton Laakso  (FMI) 2013
  ! Juha Tonttila (FMI) 2014
  !
  !---------------------------------------------------------------------


  SUBROUTINE coagulation(kproma, kbdim,  klev,    &
                         paero,  pcloud, pprecp,  &
                         ptstep, ptemp,  ppres    )

    USE mo_submctl, ONLY:        &
         t_parallelbin, t_section,   & ! Datatypes for the cloud bin representation
         in1a, fn1a,                 & ! size bin indices
         in2a, fn2a,                 &
         in2b, fn2b,                 & 
         ica,fca,icb,fcb,            &
         ncld, nprc, ira,fra,        &
         debug,                      &
         ncldbin,                    & ! Sizes of cloud and drizzle droplet regimes
         pi6,                        &
         rhosu,rhooc,rhono,rhonh,    &
         rhobc,rhodu,rhoss,rhowa,    &
         nlim,prlim,                       &
         lscgaa, lscgcc, lscgca, lscgpp, lscgpa, lscgpc
    !USE mo_aero_mem_salsa, ONLY :    &
    !     d_cond_so4,                 & ! diagnostic variable for condensated mass of so4
    !     d_nuc_so4                     !diagnostic variable for nucleated mass of so4
    !USE mo_salsa_init, only: coagc
    !USE mo_salsa_init, ONLY : coagc
    USE mo_kind, ONLY : dp

    IMPLICIT NONE


    !-- Input and output variables -------------
    INTEGER, INTENT(IN) ::          &
         kproma,                    & ! number of horiz. grid kproma 
         kbdim,                     & ! dimension for arrays 
         klev                         ! number of vertical klev 

    TYPE(t_section), INTENT(inout) :: &
         pcloud(kbdim,klev,ncld),     &  ! Hydrometeor properties
         paero(kbdim,klev,fn2b),      &  ! Aerosol properties
         pprecp(kbdim,klev,nprc)

    REAL(dp), INTENT(IN) ::         &
         ptstep,                    & ! time step [s]
         ptemp(kbdim,klev),         &
         ppres(kbdim,klev)
    !-- Local variables ------------------------
    INTEGER ::                      &
         ii,jj,kk,ll,mm,nn,cc,      & ! loop indices 
         index_2a, index_2b,        & ! corresponding bin in regime 2a/2b
         index_cd                     ! corresponding bin in cloud droplet regime

    REAL(dp) ::                     &
         zntemp(kbdim,klev),        & ! variable for saving pnaero(fn2b) temporarily
         zcc(fn2b,fn2b),            & ! updated coagulation coefficients [m3/s]
         zcccc(ncld,ncld),          & ! - '' - for collision-coalescence [m3/s]
         zccca(fn2b,ncld),          & ! - '' - for cloud collection of aerosols [m3/s]
         zccpc(ncld,nprc),          & ! - '' - for collection of cloud droplets by precip [m3/s]
         zccpa(fn2b,nprc),          & ! - '' - for collection of aerosols by precip
         zccpp(nprc,nprc),          & ! - '' - for collitions between precip particles (neglected?)
         zminusterm,                & ! coagulation loss in a bin [1/s] 
         zplusterm(8)                 ! coagulation gain in a bin [fxm/s]
                                      ! (for each chemical compound)
         !zcc(kbdim,klev,fn2b,fn2b), & ! updated coagulation coefficients [m3/s]

    !REAL(dp) ::                     &
    !     zrho(8)                     ! Table for densities of different species in hydrometeor distribution
    REAL(dp) :: &
         zmpart(fn2b),  & ! approximate mass of particles [kg]
         zmcloud(ncld), &    ! approximate mass of cloud droplets [kg]
         zmprecp(nprc)   ! Approximate mass for rain drops [kg]

    REAL(dp) :: &
         temppi,pressi,pdmm,pdnn
    REAL(dp) :: zsec(nprc) ! Security coefficient

    !REAL(dp) :: &
    !     zhvol(kbdim,klev,ncld),  &       ! Total volume concentration of hydrometeors
    !     zavol(kbdim,klev,fn2b)           ! Total volume concentration of aerosols


    REAL(dp) :: t1,t2

    IF (debug) WRITE(*,*) 'coagulation init'


    !-----------------------------------------------------------------------------
    !-- 1) Coagulation to coarse mode calculated in a simplified way: ------------
    !      CoagSink ~ Dp in continuum regime, thus we calculate
    !      'effective' number concentration of coarse particles

    zntemp = paero(:,:,fn2b)%numc

    !-- 2) Updating coagulation coefficients -------------------------------------

    IF (debug) WRITE(*,*) 'start kernels'

     DO jj = 1,klev      ! vertical grid
        DO ii = 1,kproma ! horizontal kproma in the slab

           
           ! NÄIHIN KANNATTAISI EHKÄ LASKEA SUORAAN TILAVUUSKONSENTRAATIOISTA MASSA; DWET EI VÄLTTÄMÄTTÄ AJANTASAINEN
           !-- particle mass; density of 1500 kg/m3 assumed [kg] 
           zmpart = pi6*(MIN(paero(ii,jj,1:fn2b)%dwet, 30.e-6_dp)**3)*1500._dp

           !-- Cloud mass; Assume water density
           zmcloud(1:ncld) = pi6*(pcloud(ii,jj,1:ncld)%dwet**3)*rhowa

           !-- Precipitation mass
           zmprecp(1:nprc) = pi6*(MIN(pprecp(ii,jj,1:nprc)%dwet, 2.e-3_dp)**3)*rhowa

           temppi=ptemp(ii,jj)
           pressi=ppres(ii,jj)
           zcc = 0._dp
           zcccc = 0._dp
           zccca = 0._dp
           zccpp = 0._dp
           zccpc = 0._dp
           zccpa = 0._dp

           zsec(:) = MERGE(1._dp,0._dp,pprecp(ii,jj,:)%numc > 1.e-10_dp) ! Typical "physical" minimum number concentration for LES
           ! Aero-aero coagulation
           IF (lscgaa) THEN
              pdmm = 0._dp
              pdnn = 0._dp
              DO mm = 1,fn2b         ! smaller colliding particle
                 DO nn = mm,fn2b            ! larger colliding particle 
                    pdmm=MIN(paero(ii,jj,mm)%dwet, 30.e-6_dp)
                    pdnn=MIN(paero(ii,jj,nn)%dwet, 30.e-6_dp)
                    zcc(mm,nn) = coagc(pdmm,pdnn,zmpart(mm),zmpart(nn),temppi,pressi,1)
                    zcc(nn,mm) = zcc(mm,nn)
                 END DO
              END DO
           END IF

           ! Collision-coalescence between cloud droplets
           IF (lscgcc .AND. ANY(pcloud(ii,jj,:)%numc > nlim)) THEN
              pdmm = 0._dp
              pdnn = 0._dp
              DO mm = 1,ncld
                 DO nn = mm,ncld 
                    pdmm = pcloud(ii,jj,mm)%dwet
                    pdnn = pcloud(ii,jj,nn)%dwet
                    zcccc(mm,nn) = coagc(pdmm,pdnn,zmcloud(mm),zmcloud(nn),temppi,pressi,2)
                    zcccc(nn,mm) = zcccc(mm,nn)
                 END DO
              END DO
           END IF

           ! Self-collection of rain drops
           IF (lscgpp .AND. ANY(pprecp(ii,jj,:)%numc > prlim)) THEN
              pdmm = 0._dp
              pdnn = 0._dp
              DO mm = 1,nprc
                 DO nn = mm,nprc
                    pdmm = MIN(pprecp(ii,jj,mm)%dwet,2.e-3_dp)
                    pdnn = MIN(pprecp(ii,jj,nn)%dwet,2.e-3_dp)
                    zccpp(mm,nn) =  coagc(pdmm,pdnn,zmprecp(mm),zmprecp(nn),temppi,pressi,2)*zsec(mm)*zsec(nn)
                    zccpp(nn,mm) = zccpp(mm,nn)
                 END DO
              END DO
           END IF

           ! Cloud collection of aerosols
           IF (lscgca .AND. ANY(pcloud(ii,jj,:)%numc > nlim)) THEN
              DO mm = 1,fn2b
                 DO nn = 1,ncld
                    pdmm = MIN(paero(ii,jj,mm)%dwet, 30.e-6_dp)
                    pdnn = pcloud(ii,jj,nn)%dwet
                    zccca(mm,nn) = coagc(pdmm,pdnn,zmpart(mm),zmcloud(nn),temppi,pressi,2)
                 END DO
              END DO
           END IF 
           
           ! Collection of aerosols by rain
           IF (lscgpa .AND. ANY(pprecp(ii,jj,:)%numc > prlim)) THEN
              DO mm = 1,fn2b
                 DO nn = 1,nprc
                    pdmm = MIN(paero(ii,jj,mm)%dwet,30.e-6_dp)
                    pdnn = MIN(pprecp(ii,jj,nn)%dwet, 2.e-3_dp)
                    zccpa(mm,nn) = coagc(pdmm,pdnn,zmpart(mm),zmprecp(nn),temppi,pressi,2)*zsec(nn)
                 END DO
              END DO
           END IF

           ! Collection of cloud droplets by rain
           IF (lscgpc .AND. (ANY(pcloud(ii,jj,:)%numc > nlim) .AND. ANY(pprecp(ii,jj,:)%numc > prlim)) ) THEN 
              DO mm = 1,ncld
                 DO nn = 1,nprc
                    pdmm = pcloud(ii,jj,mm)%dwet
                    pdnn = MIN(pprecp(ii,jj,nn)%dwet,2.e-3_dp)
                    zccpc(mm,nn) = coagc(pdmm,pdnn,zmcloud(mm),zmprecp(nn),temppi,pressi,2)!*zsec(nn)
                  END DO
              END DO
           END IF

           !-- 3) New particle and volume concentrations after coagulation -------------
           
           ! Aerosols in regime 1a
           ! --------------------------------
           DO kk = in1a,fn1a

              zminusterm = 0._dp
              zplusterm(:) = 0._dp
              ! Particles lost by coagulation with larger aerosols
              DO ll = kk+1,fn2b
                 zminusterm = zminusterm + zcc(kk,ll)*paero(ii,jj,ll)%numc
              END DO

              ! Particles lost by cloud collection
              DO ll = 1,ncld
                 zminusterm = zminusterm + zccca(kk,ll)*pcloud(ii,jj,ll)%numc
              END DO

              ! particles lost by rain collection
              DO ll = 1,nprc
                 zminusterm = zminusterm + zccpa(kk,ll)*pprecp(ii,jj,ll)%numc
              END DO
              
              DO ll = in1a,kk-1
                 zplusterm(1:2) = zplusterm(1:2) + zcc(ll,kk)*paero(ii,jj,ll)%volc(1:2)
                 zplusterm(6:7) = zplusterm(6:7) + zcc(ll,kk)*paero(ii,jj,ll)%volc(6:7)
                 zplusterm(8) = zplusterm(8) + zcc(ll,kk)*paero(ii,jj,ll)%volc(8)
              END DO

              !-- Volume and number concentrations after coagulation update [fxm]
              paero(ii,jj,kk)%volc(1:2) = ( paero(ii,jj,kk)%volc(1:2)+ptstep*zplusterm(1:2) * &
                   paero(ii,jj,kk)%numc ) / (1._dp + ptstep*zminusterm)
              
              paero(ii,jj,kk)%volc(6:7) = ( paero(ii,jj,kk)%volc(6:7)+ptstep*zplusterm(6:7) * &
                   paero(ii,jj,kk)%numc ) / (1._dp + ptstep*zminusterm)

              paero(ii,jj,kk)%volc(8) = ( paero(ii,jj,kk)%volc(8)+ptstep*zplusterm(8) * &
                   paero(ii,jj,kk)%numc ) / ( 1._dp + ptstep*zminusterm)

              paero(ii,jj,kk)%numc = paero(ii,jj,kk)%numc/(1. + ptstep*zminusterm  + &
                   0.5_dp*ptstep*zcc(kk,kk)*paero(ii,jj,kk)%numc)
 
           END DO

           ! Aerosols in regime 2a
           ! ---------------------------------
           DO kk = in2a,fn2a
              
              zminusterm = 0._dp
              zplusterm(:) = 0._dp
              
              ! Find corresponding size bin in subregime 2b
              index_2b = kk - in2a + in2b
              
              ! Particles lost by larger particles in 2a
              DO ll = kk+1, fn2a
                 zminusterm = zminusterm + zcc(kk,ll)*paero(ii,jj,ll)%numc ! 2a
              END DO

              ! Particles lost by larger particles in 2b
              DO ll = index_2b+1, fn2b
                 zminusterm = zminusterm + zcc(kk,ll)*paero(ii,jj,ll)%numc ! 2b 
              END DO

              ! Particles lost by cloud collection
              DO ll = 1,ncld
                 zminusterm = zminusterm + zccca(kk,ll)*pcloud(ii,jj,ll)%numc
              END DO

              ! Particles lost by collection by rain
              DO ll = 1,nprc
                 zminusterm = zminusterm + zccpa(kk,ll)*pprecp(ii,jj,ll)%numc
              END DO

              ! Particle volume gained from smaller particles in regimes 1, 2a and 2b
              DO ll = in1a, kk-1
                 zplusterm(1:2) = zplusterm(1:2) + zcc(ll,kk)*paero(ii,jj,ll)%volc(1:2) 
                 zplusterm(6:7) = zplusterm(6:7) + zcc(ll,kk)*paero(ii,jj,ll)%volc(6:7) ! NO + NH
                 zplusterm(8) = zplusterm(8) + zcc(ll,kk)*paero(ii,jj,ll)%volc(8) ! for h2o
              END DO
                   
              ! Particle volume gained from smaller particles in 2a
              ! (Note, for components not included in the previous loop!)
              DO ll = in2a, kk-1
                 ! Don't do water twice!
                 zplusterm(3:5) = zplusterm(3:5) + zcc(ll,kk)*paero(ii,jj,ll)%volc(3:5)             
              END DO
              
              ! Particle volume gained from smaller (and equal) particles in 2b
              DO ll = in2b, index_2b
                 zplusterm(1:8) = zplusterm(1:8) + zcc(ll,kk)*paero(ii,jj,ll)%volc(1:8) ! 2b
              END DO

              !-- Volume and number concentrations after coagulation update [fxm]
              paero(ii,jj,kk)%volc(1:8) = ( paero(ii,jj,kk)%volc(1:8)+ptstep*zplusterm(1:8) *  &
                   paero(ii,jj,kk)%numc ) / (1._dp + ptstep*zminusterm)
              
              paero(ii,jj,kk)%numc = paero(ii,jj,kk)%numc/(1. + ptstep*zminusterm  + &
                   0.5_dp*ptstep*zcc(kk,kk)*paero(ii,jj,kk)%numc)

           END DO

           ! Aerosols in regime 2b
           ! ---------------------------------
           DO kk = in2b,fn2b

              zminusterm = 0._dp
              zplusterm(:) = 0._dp

              !-- Find corresponding size bin in subregime 2a
              index_2a = kk - in2b + in2a

              ! Particles lost to larger particles in regimes 2b
              DO ll = kk+1, fn2b
                 zminusterm = zminusterm + zcc(kk,ll)*paero(ii,jj,ll)%numc ! 2b                     
              END DO

              ! Particles lost to larger and equal particles in 2a
              DO ll = index_2a, fn2a                       
                 zminusterm = zminusterm + zcc(kk,ll)*paero(ii,jj,ll)%numc                 
              END DO
                    
              ! Particles lost by cloud collection
              DO ll = 1,ncld
                 zminusterm = zminusterm + zccca(kk,ll)*pcloud(ii,jj,ll)%numc
              END DO
              
              ! Particles lost by collection by rain
              DO ll = 1,nprc
                 zminusterm = zminusterm + zccpa(kk,ll)*pprecp(ii,jj,ll)%numc
              END DO

              ! Particle volume gained from smaller particles in 1/2a
              DO ll = in1a, index_2a-1
                 zplusterm(1:2) = zplusterm(1:2) + zcc(ll,kk)*paero(ii,jj,ll)%volc(1:2)
                 zplusterm(6:7) = zplusterm(6:7) + zcc(ll,kk)*paero(ii,jj,ll)%volc(6:7)
                 zplusterm(8) = zplusterm(8) + zcc(ll,kk)*paero(ii,jj,ll)%volc(8)
              END DO
              DO ll = in2a, index_2a-1
                 ! Don't do water twice!
                 zplusterm(3:5) = zplusterm(3:5) + zcc(ll,kk)*paero(ii,jj,ll)%volc(3:5)
              END DO

              ! Particle volume gained from smaller particles in 2b 
              DO ll = in2b, kk-1
                 zplusterm(1:8) = zplusterm(1:8) + zcc(ll,kk)*paero(ii,jj,ll)%volc(1:8)
              END DO

              !-- Volume and number concentrations after coagulation update [fxm]
              paero(ii,jj,kk)%volc(1:8) = ( paero(ii,jj,kk)%volc(1:8)+ptstep*zplusterm(1:8) *  &
                   paero(ii,jj,kk)%numc ) / (1._dp + ptstep*zminusterm)
              
              paero(ii,jj,kk)%numc = paero(ii,jj,kk)%numc/(1. + ptstep*zminusterm  + &
                   0.5_dp*ptstep*zcc(kk,kk)*paero(ii,jj,kk)%numc)
              
           END DO

           ! Cloud droplets, regime a
           ! ------------------------------------------------
           IF ( ANY(pcloud(ii,jj,:)%numc > nlim) ) THEN
              DO cc = ica%cur,fca%cur
              
                 zminusterm = 0._dp
                 zplusterm(:) = 0._dp

                 ! corresponding index for regime b cloud droplets
                 kk = MAX(cc-fca%cur+ncld,icb%cur) ! Regime a has more bins than b: 
                                                   ! Set this at minimum to beginnign of b.

                 ! Droplets lost by those with larger nucleus in regime a
                 DO ll = cc+1,fca%cur
                    zminusterm = zminusterm + zcccc(cc,ll)*pcloud(ii,jj,ll)%numc
                 END DO

                 ! Droplets lost by those with larger nucleus in regime b
                 DO ll = kk+1,fcb%cur
                    zminusterm = zminusterm + zcccc(cc,ll)*pcloud(ii,jj,ll)%numc
                 END DO

                 ! Droplets lost by collection by rain drops
                 DO ll = 1,nprc
                    zminusterm = zminusterm + zccpc(cc,ll)*pprecp(ii,jj,ll)%numc
                 END DO
              
                 ! Volume gained from cloud collection of aerosols
                 DO ll = in1a,fn2b
                    zplusterm(1:8) = zplusterm(1:8) + zccca(ll,cc)*paero(ii,jj,ll)%volc(1:8)
                 END DO

                 ! Volume gained from smaller droplets in a
                 DO ll = ica%cur,cc-1
                    zplusterm(1:8) = zplusterm(1:8) + zcccc(ll,cc)*pcloud(ii,jj,ll)%volc(1:8)
                 END DO

                 ! Volume gained from smaller or equal droplets in b
                 DO ll = icb%cur,kk
                    zplusterm(1:8) = zplusterm(1:8) + zcccc(ll,cc)*pcloud(ii,jj,ll)%volc(1:8)
                 END DO

                 ! Update the hydrometeor volume concentrations
                 pcloud(ii,jj,cc)%volc(1:8) = ( pcloud(ii,jj,cc)%volc(1:8) +  &
                      ptstep*zplusterm(1:8)*pcloud(ii,jj,cc)%numc ) /         &
                      (1._dp + ptstep*zminusterm)
                    
                 ! Update the hydrometeor number concentration (Removal by coagulation with lrger bins and self)      
                 pcloud(ii,jj,cc)%numc = pcloud(ii,jj,cc)%numc/( 1._dp + ptstep*zminusterm +  &
                      0.5_dp*ptstep*zcccc(cc,cc)*pcloud(ii,jj,cc)%numc )

              END DO

              ! Cloud droplets, regime b
              ! -----------------------------------------
              DO cc = icb%cur,fcb%cur
              
                 zminusterm = 0._dp
                 zplusterm(:) = 0._dp
                 
                 ! corresponding index for regime a cloud droplets
                 kk = cc - ncld + fca%cur  

                 ! Droplets lost by those with larger nucleus in regime b
                 DO ll = cc+1,fcb%cur
                    zminusterm = zminusterm + zcccc(cc,ll)*pcloud(ii,jj,ll)%numc
                 END DO

                 ! Droplets lost by those with larger nucleus in regime a
                 DO ll = kk+1,fca%cur
                    zminusterm = zminusterm + zcccc(cc,ll)*pcloud(ii,jj,ll)%numc
                 END DO

                 ! Droplets lost by collection by rain drops
                 DO ll = 1,nprc
                    zminusterm = zminusterm + zccpc(cc,ll)*pprecp(ii,jj,ll)%numc
                 END DO
              
                 ! Volume gained from cloud collection of aerosols
                 DO ll = in1a,fn2b
                    zplusterm(1:8) = zplusterm(1:8) + zccca(ll,cc)*paero(ii,jj,ll)%volc(1:8)
                 END DO

                 ! Volume gained from smaller droplets in b
                 DO ll = icb%cur,cc-1
                    zplusterm(1:8) = zplusterm(1:8) + zcccc(ll,cc)*pcloud(ii,jj,ll)%volc(1:8)
                 END DO

                 ! Volume gained from smaller or equal droplets in a
                 DO ll = ica%cur,kk
                    zplusterm(1:8) = zplusterm(1:8) + zcccc(ll,cc)*pcloud(ii,jj,ll)%volc(1:8)
                 END DO

                 ! Update the hydrometeor volume concentrations
                 pcloud(ii,jj,cc)%volc(1:8) = ( pcloud(ii,jj,cc)%volc(1:8) +  &
                      ptstep*zplusterm(1:8)*pcloud(ii,jj,cc)%numc ) /         &
                      (1._dp + ptstep*zminusterm)
                    
                 ! Update the hydrometeor number concentration (Removal by coagulation with lrger bins and self)      
                 pcloud(ii,jj,cc)%numc = pcloud(ii,jj,cc)%numc/( 1._dp + ptstep*zminusterm +  &
                      0.5_dp*ptstep*zcccc(cc,cc)*pcloud(ii,jj,cc)%numc )

              END DO
           END IF ! nlim

           ! Rain drops
           ! -----------------------------------
           IF ( ANY(pprecp(ii,jj,:)%numc > prlim) ) THEN
              DO cc = 1,nprc
              
                 zminusterm = 0._dp
                 zplusterm(:) = 0._dp

                 ! Drops lost by coagulation with larger drops
                 DO ll = cc+1,nprc
                    zminusterm = zminusterm + zccpp(cc,ll)*pprecp(ii,jj,ll)%numc
                 END DO

                 ! Volume gained by collection of aerosols
                 DO ll = in1a,fn2b
                    zplusterm(1:8) = zplusterm(1:8) + zccpa(ll,cc)*paero(ii,jj,ll)%volc(1:8) 
                 END DO
              
                 ! Volume gained by collection of cloud droplets
                 DO ll = 1,ncld
                    zplusterm(1:8) = zplusterm(1:8) + zccpc(ll,cc)*pcloud(ii,jj,ll)%volc(1:8)
                 END DO
               
                 ! Volume gained from smaller drops
                 IF (cc > 1) THEN
                    DO ll = 1,cc-1
                       zplusterm(1:8) = zplusterm(1:8) + zccpp(ll,cc)*pprecp(ii,jj,ll)%volc(1:8)
                    END DO
                 END IF

                 ! Update the hydrometeor volume concentrations
                 pprecp(ii,jj,cc)%volc(1:8) = ( pprecp(ii,jj,cc)%volc(1:8) +  &
                      ptstep*zplusterm(1:8)*pprecp(ii,jj,cc)%numc ) /         &
                      (1._dp + ptstep*zminusterm)
                 
                 ! Update the hydrometeor number concentration (Removal by coagulation with lrger bins and self)      
                 pprecp(ii,jj,cc)%numc = pprecp(ii,jj,cc)%numc/( 1._dp + ptstep*zminusterm +  &
                      0.5_dp*ptstep*zccpp(cc,cc)*pprecp(ii,jj,cc)%numc )
                 
              END DO
           END IF  !prlim

        END DO ! kbdim
     END DO ! klev
     ! fxm: here we approximate that the sea salt regime 2b particles have
     ! gained by coagulation can be treated as sulphate
     paero(:,:,in2b:fn2b)%volc(1) = paero(:,:,in2b:fn2b)%volc(1) + paero(:,:,in2b:fn2b)%volc(5)
     paero(:,:,in2b:fn2b)%volc(5) = 0._dp
     
     DO cc = 1,8
        paero(:,:,in2b:fn2b)%volc(cc) = MAX(paero(:,:,in2b:fn2b)%volc(cc), 0._dp)
     END DO

  END SUBROUTINE coagulation


  ! fxm: calculated for empty bins too
  ! fxm: same diffusion coefficients and mean free paths used for sulphuric acid
  !      and organic vapours (average values? 'real' values for each?)
  !********************************************************************
  !
  ! subroutine CONDENSATION(kproma, kbdim,  klev,        &
  !                         pnaero, pvols,  pdwet, plwc, &
  !                         pcsa,   pcocnv, pcocsv,      &
  !                         ptemp,  ppres,  ptstep)
  !
  !********************************************************************
  !
  ! Purpose:
  ! --------
  ! Calculates the increase in particle volume and 
  !  decrease in gas phase concentrations due to condensation 
  !  of sulphuric acid and two organic compounds (non-volatile
  !  and semivolatile)
  !
  !
  ! Method:
  ! -------
  ! Regime 3 particles only act as a sink for condensing vapours
  !  while their size and composition does not change.
  ! Exception: Soluble fraction of regime 3c particles can change
  !  and thus they can be moved to regime 3b 
  !
  ! New gas and aerosol phase concentrations calculated according
  !  to Jacobson (1997): Numerical techniques to solve 
  !  condensational and dissolutional growth equations 
  !  when growth is coupled to reversible reactions, 
  !  Aerosol Sci. Tech., 27, pp 491-498.
  !
  ! fxm: one should really couple with vapour production and loss terms as well
  !      should nucleation be coupled here as well????
  !
  ! Juha: Now does the condensation of water vapour on hydrometeors as well,
  !       + the condensation of semivolatile aerosol species on hydromets.
  !       Modified for the new aerosol datatype. LWC is obtained from %volc(8)
  !
  !
  ! Interface:
  ! ----------
  ! Called from main aerosol model
  ! 
  !
  ! Coded by:
  ! ---------
  ! Hannele Korhonen (FMI) 2005 
  ! Harri Kokkola (FMI) 2006
  ! Juha Tonttila (FMI) 2014
  !
  !---------------------------------------------------------------
  !
  ! Following parameterization has been used:
  ! ------------------------------------------
  !
  ! Molecular diffusion coefficient of condensing vapour [m2/s]
  !  (Reid et al. (1987): Properties of gases and liquids,
  !   McGraw-Hill, New York.)
  !
  ! D = {1.d-7*sqrt(1/M_air + 1/M_gas)*T^1.75} / &
  !  {p_atm/p_stand * (d_air^(1/3) + d_gas^(1/3))^2 }
  !   
  ! M_air = 28.965 : molar mass of air [g/mol]
  ! d_air = 19.70  : diffusion volume of air
  ! M_h2so4 = 98.08  : molar mass of h2so4 [g/mol]
  ! d_h2so4 = 51.96  : diffusion volume of h2so4
  !
  !---------------------------------------------------------------

  SUBROUTINE condensation(kproma,  kbdim,  klev,   krow,      &
                          paero,   pcloud, pprecp, pcsa,      &
                          pcocnv,  pcocsv, pchno3, pcnh3,     &
                          prv,prs, ptemp,  ppres,  ptstep,    &
                          ppbl,    prtcl                      )

    USE mo_salsa_nucleation

    USE mo_submctl,    ONLY :   &
         pi,                        & 
         pi6,                       & ! pi/6 
         in1a, in2a,                & ! size bin indices
         fn1a,                       &
         fn2a, fn2b,                & 
         nbin,                      & ! number of size bins in each regime
         nbins,                     & ! Total number of aerosol bins

         t_section,                 & ! Data type for the cloud bin representation
         ica,fca,icb,fcb,           & ! bin indices for cloud/rain bins
         ira,fra,                   &
         ncld,                      &
         nprc,                      &
         debug,                     &
         lscndgas,                  &
         !lscndh2o,                  &
         avog,                      &
         nlim,                      &
         prlim,                     &
         rhowa,                     & ! density of water (kg/m3)
         rhosu,                     & ! density of sulphate (kg/m3)
         rhooc,                     & ! density of organic carbon (kg/m3)
         rhoss,                     & ! density of sea salt (kg/m3)
         rhono,                     & ! density of nitric acid (kg/m3)
         rhonh,                     & ! density of ammonia (kg/m3)
         rhobc,                     &
         rhodu,                      &
         boltz,                     & ! Boltzmann constant [J/K]
         rg,                        & ! molar gas constant [J/(mol K)]
         pstand,                    & ! standard pressure [Pa]
         msu,                       & ! molar mass of sulphate [kg/mol]
         moc,                       & !       "       organic carbon
         mss,                       & !       "       sea salt
         mno,                       & !       "       nitrate
         mnh,                       & !       "       ammonium
         mbc,                       & ! 
         mdu,                       &
         mwa,                       & !               water
         mair,                      & !       "       air
         mvsu, mvoc,                & ! molecular volumes of sulphate and OC [m3]
         mvnh, mvno, mvwa,          & ! molecular volumes of HNO3 and NH3,H20 [m3]
         d_sa,                      & ! diameter of H2SO4 molecule [m]
         
         epsoc,                     & ! soluble fraction of organics (scaled to sulphate)
         massacc,                   & ! mass accomodation coefficients in each bin
         n3,                        & ! number of molecules in one 3 nm particle [1]
         nsnucl,                    & ! nucleation
         surfw0                       ! surface tension of water

    USE mo_kind, ONLY : dp

    USE class_componentIndex, ONLY : ComponentIndex,IsUsed

    !USE mo_time_control,   ONLY: delta_time,time_step_len

    USE mo_constants,      ONLY: g, avo, alv, rv

   IMPLICIT NONE

    !-- Input and output variables ----------
    INTEGER, INTENT(IN) ::          &
         kproma,                    & ! number of horiz. grid kproma 
         kbdim,                     & ! dimension for arrays 
         klev,                      & ! number of vertical klev 
         krow

    REAL(dp), INTENT(IN) ::         &  
         ptemp(kbdim,klev),         & ! ambient temperature [K]
         ppres(kbdim,klev),         & ! ambient pressure [Pa]
         ptstep,                    & ! timestep [s]
         prs(kbdim,klev)              ! Water vapor saturation mixing ratio

    TYPE(ComponentIndex), INTENT(in) :: prtcl  ! Keeps track which substances are used

    !LOGICAL, INTENT(in) :: pactmask(kbdim,klev)

    INTEGER :: ppbl(kbdim)           ! Planetary boundary layer top level

    REAL(dp), INTENT(INOUT) ::     &
         !prh(kbdim,klev),          & ! Juha: Moved from above
         prv(kbdim,klev),          & ! Water vapor mixing ratio
         pcsa(kbdim,klev),         & ! sulphuric acid concentration [#/m3]
         pcocnv(kbdim,klev),       & ! non-volatile organic concentration [#/m3]
         pcocsv(kbdim,klev),       & ! semivolatile organic concentration [#/m3]
         pchno3(kbdim,klev),       & ! nitric acid concentration [#/m3]
         pcnh3(kbdim,klev)           ! ammonia concentration [#/m3]

    TYPE(t_section), INTENT(inout) :: &
         pcloud(kbdim,klev,ncld),     & ! Hydrometeor properties
         paero(kbdim,klev,fn2b),      & ! Aerosol properties
         pprecp(kbdim,klev,nprc)



    !-- Local variables ----------------------
    INTEGER :: ii, jj, kk, ll, mm, cc    ! loop indices

    REAL(dp) ::                      &
         zvisc,                      & ! viscosity of air [kg/(m s)]
         zdfvap,                     & ! air diffusion coefficient [m2/s]
         zmfp,                       & ! mean free path of condensing vapour [m]
         zcs_tot,                    & ! total condensation sink [1/s] (gases)
         zcs_ocsv,                   & ! condensation sink for semivolatile organics [1/s]
         zcs_su,                     & ! condensation sink for sulfate [1/s]
         zcs_ocnv,                   & ! condensation sink for nonvolatile organics [1/s]
         zcs_no,                     & ! condensation sink for NO3 [1/s]
         zcs_nh,                     & ! condensation sink for NH3 [1/s]
         zcs_h2o,                    & ! Condensation sink for H2O [1/s]
                                       ! vapour concentration after time step [#/m3]
         zcvap_new1,                 & ! sulphuric acid
         zcvap_new2,                 & ! nonvolatile organics
         zcvap_new3,                 & ! semivolatile organics
         zcvap_newno,                & ! HNO3
         zcvap_newnh,                & ! NH3
                                       ! change in vapour concentration [#/m3]
         zdvap1,                     & ! sulphuric acid
         zdvap2,                     & ! nonvolatile organics
         zdvap3,                     & ! semivolatile organics
         zdvapno,                    & ! HNO3
         zdvapnh,                    & ! NH3
         
         zdfpart(in1a+1),            & ! particle diffusion coefficient
         zdfh2o,                     & ! diffusion coefficient for h2o vapour in air

         zknud(fn2b),                & ! particle Knudsen number
         zknaw(fn2b),                & ! Kundsen number for aerosols and water vapour
         zkncw(ncld),                & ! Knudsen number for cloud droplets and water vapour
         zknca(ncld),                & ! Knudsen number for cloud droplets and aerosol vapours
         zknpw(nprc),                & ! Knudsen number for rain drops and water vapour
         zknpa(nprc),                & ! Knudsen number for rain drops and aerosol vapours

         zbeta(fn2b),                & ! transitional correction factor for aerosols
         zbetaaw(fn2b),              & ! - '' - for water condensing on aerosols 
         zbetacw(ncld),              & ! transitional correction for condensing water vapour on clouds (is this even needed?)
         zbetaca(ncld),              & ! - '' - for condensing aerosol vapours on clouds (is this needed?)
         zbetapw(nprc),              & ! - '' - for condensing water on rain drops
         zbetapa(nprc),              & ! - '' - for condensing aerosol vapours on rain drops

         zcolrate(fn2b),             & ! collision rate of molecules to particles [1/s]
         zcolrate_ocnv(fn2b),        & ! collision rate of organic molecules to particles [1/s]
         zcolrateaw(fn2b),           & ! Collision rate of water vapour on aerosols
         zcolratecw(ncld),           & ! Collision rate of water vapour to cloud drops
         zcolrateca(ncld),           & ! Collision rate of aerosol vapour molecules to cloud drops
         zcolratepw(nprc),           & ! Collision rate of water vapour to rain drops
         zcolratepa(nprc),           & ! Collision rate of gases to rain drops


         zmtae(fn2b),                & ! Mass transfer coefficients for aerosols
         zmtcd(ncld),                & ! Mass transfer coefficients for cloud droplets
         zmtpd(nprc),                & ! Mass transfer coefficients for rain drops
         zmttot,                     & ! Total mass transfer coefficient
         
         zdvolsa(fn2b),              & ! change of sulphate volume in each bin [fxm]
         zdvoloc(fn2b),              & !    - " - organics 

         zj3n3(kbdim,klev,2),        & ! Formation massrate of molecules in nucleation, [molec/m3s].  (kbdim,klev,1) for H2SO4 and (kbdim,klev,2) for Organic vapor
         zn_vs_c,                    & ! ratio of nucleation of all mass transfer in the smallest bin
         zxsa(kbdim,klev),           & ! ratio of sulphuric acid and organic vapor in 3nm particles 
         zxocnv(kbdim,klev),         & !
         zmaxdvol(fn2b),             & ! 
         zcno_tot,                   & ! total available NO3
         zcnh3_tot,                  & ! total available NH3         

         zthcond               ! thermal conductivity of air
         !zmfph2o                   ! Mean free path for h2o vapour 

    REAL(dp) :: zns, znw    ! Number of moles of solute and water in particles

    REAL(dp) :: zhlp1,zhlp2,zhlp3  ! Helper variables
    REAL(dp) :: zhlp1a(ncld), zhlp2a(ncld), zhlp3a(ncld)
    REAL(dp) :: zhlp1b(nbins), zhlp2b(nbins), zhlp3b(nbins)

    real(dp)::ztmst,zqtmst,        &
         ztstep2, ztstep_pre,zthno3,ztnh3,zph2o,zphno3,zpnh3,zaw

    REAL(dp) :: zrh(kbdim,klev)
    
   ! REAL(dp) :: zbetann(kbdim,klev,nbins) ! Correction factor for nitrate calculations

                      ! FROM/TO ISOROPIA
         REAL(dp) :: zwi_S(5)          ! Total species concentrations in moles/m**3 air
         REAL(dp) :: zcntrl_s(2)       ! nug for different types of problems solved
                                       ! different state of aerosols (deliquescent or
                                       ! metastable)
         REAL(dp) :: zwt_s(1)          ! ?
         REAL(dp) :: zaerliq_S(12)     ! Aqueous-phase concentration array in moles/m**3air
         REAL(dp) :: zaersld_S(9)      ! Solid-phase concentration array in moles/m**3 air  
         REAL(dp) :: zother_s(6)       ! Solution information array
         REAL(dp) :: zgas_s(3)         ! Gas-phase concentration array in moles/m**3 air
         CHARACTER*15 scasi_s          ! Returns the subcase which the input corresponds to
    INTEGER, PARAMETER :: NCmax=3, NAmax=5, NNmax=1
    REAL(dp) :: MOLAL(-NAmax:NCmax)    !-- Initializations
    INTEGER :: tt




    zj3n3 = 0._dp
    zrh(1:kproma,:) = prv(1:kproma,:)/prs(1:kproma,:)
    
    !------------------------------------------------------------------------------

    IF (nsnucl > 0) CALL nucleation(kproma, kbdim,  klev,   krow,   &
                                    paero,  ptemp,  zrh,    ppres,  &
                                    pcsa,   pcocnv, ptstep, zj3n3,  &
                                    zxsa,   zxocnv, ppbl            )
    zdvolsa=0._dp
    zn_vs_c=0._dp
    DO jj = 1,klev
       DO ii = 1,kproma

          zdvoloc = 0._dp

          !-- 1) Properties of air and condensing gases --------------------
          zvisc  = (7.44523e-3_dp*ptemp(ii,jj)**1.5_dp)/(5093._dp*(ptemp(ii,jj)+110.4_dp))! viscosity of air [kg/(m s)] 
          zdfvap = 5.1111e-10_dp*ptemp(ii,jj)**1.75_dp*pstand/ppres(ii,jj)                ! diffusion coefficient [m2/s]
          zmfp   = 3._dp*zdfvap*sqrt(pi*msu/(8._dp*rg*ptemp(ii,jj)))                      ! mean free path [m]
         
          !zmfph2o = 3._dp*zdfh2o*sqrt(pi*mwa/(8._dp*rg*ptemp(ii,jj)))                                ! Mean free path fro h2o 

          !-- 2) Transition regime correction factor for particles ---------
          !  
          !  Fuchs and Sutugin (1971), In: Hidy et al. (ed.)
          !  Topics in current aerosol research, Pergamon.  
          !
          !  Size of condensing molecule considered only for 
          !  nucleation mode (3 - 20 nm) 
          !

          !-- particle Knudsen numbers
          zknud(in1a:in1a+1) = 2._dp*zmfp/(paero(ii,jj,in1a:in1a+1)%dwet+d_sa)              ! Gases on aerosols
          zknud(in1a+2:fn2b) = 2._dp*zmfp/paero(ii,jj,in1a+2:fn2b)%dwet

          zknca(1:ncld) = 2._dp*zmfp/pcloud(ii,jj,1:ncld)%dwet          ! Knudsen number for gases on cloud drplets

          zknpa(1:nprc) = 2._dp*zmfp/pprecp(ii,jj,1:nprc)%dwet          ! Knudsen number for gases on rain drops

          !-- transitional correction factor
          zbeta = (zknud + 1.)/(0.377_dp*zknud+1._dp+4._dp/ &     ! Aerosol + gas
                  (3._dp*massacc)*(zknud+zknud**2))  
          
          zbetaca = 1._dp + zknca*( 1.33_dp + (0.71_dp/zknca) )/( 1._dp + (1._dp/zknca) ) ! Hydrometeor + gas
          zbetaca = 1._dp/zbetaca

          zbetapa = 1._dp + zknpa*( 1.33_dp + (0.71_dp/zknpa) )/( 1._dp + (1._dp/zknpa) ) ! Rain drop + gas
          zbetapa = 1._dp/zbetapa

          !-- 3) Collision rate of molecules to particles -------------------
          !
          !  Particle diffusion coefficient considered only for 
          !  nucleation mode (3 - 20 nm) 
          !

          !-- particle diffusion coefficient [m2/s]
          zdfpart = boltz*ptemp(ii,jj)*zbeta(in1a:in1a+1)/ &    
                    (3._dp*pi*zvisc*paero(ii,jj,in1a:in1a+1)%dwet)  

          !-- collision rate (gases on aerosols) [1/s]
          zcolrate = 0._dp
          IF (lscndgas) &
               zcolrate(in1a:in1a+1) = MERGE( 2._dp*pi*(paero(ii,jj,in1a:in1a+1)%dwet+d_sa)*    & 
                                              (zdfvap+zdfpart)*zbeta(in1a:in1a+1)*              &
                                              paero(ii,jj,in1a:in1a+1)%numc,                    &
                                              0._dp,                                            &
                                              paero(ii,jj,in1a:in1a+1)%numc > nlim              )
         
          IF (lscndgas) &
               zcolrate(in1a+2:fn2b) = MERGE( 2._dp*pi*paero(ii,jj,in1a+2:fn2b)%dwet*zdfvap*    &
                                              zbeta(in1a+2:fn2b)*paero(ii,jj,in1a+2:fn2b)%numc, &
                                              0._dp,                                            &
                                              paero(ii,jj,in1a+2:fn2b)%numc > nlim              )

          !-- gases on hydrometeors
          zcolrateca = 0._dp
          IF (lscndgas) &
               zcolrateca(1:ncld) = MERGE( 2._dp*pi*pcloud(ii,jj,1:ncld)%dwet*zdfvap*    &
                                           zbetaca(1:ncld)*pcloud(ii,jj,1:ncld)%numc,    &
                                           0._dp,                                        &
                                           pcloud(ii,jj,1:ncld)%numc > nlim              )
          
          ! Gases on rain drops
          zcolratepa = 0._dp
          IF (lscndgas) &
               zcolratepa(1:nprc) = MERGE( 2._dp*pi*pprecp(ii,jj,1:nprc)%dwet*zdfvap*    &
                                           zbetapa(1:nprc)*pprecp(ii,jj,1:nprc)%numc,    &
                                           0._dp,                                        &
                                           pprecp(ii,jj,1:nprc)%numc > prlim             )

          !-- 4) Condensation sink [1/s] -------------------------------------

          zcs_tot = sum(zcolrate) + sum(zcolrateca) + sum(zcolratepa)  ! total sink


          !-- 5) Changes in gas-phase concentrations and particle volume -----
          !
          !--- 5.1) Organic vapours ------------------------

          !---- 5.1.1) Non-volatile organic compound: condenses onto all bins 
          IF(pcocnv(ii,jj) > 1.e-10_dp .and. zcs_tot > 1.e-30_dp .AND. IsUsed(prtcl,'OC')) THEN

             zn_vs_c = 0._dp

             IF(zj3n3(ii,jj,2) > 1._dp) zn_vs_c = (zj3n3(ii,jj,2))/(zj3n3(ii,jj,2) + &
                                                  pcocnv(ii,jj) * zcolrate(in1a))

             !   collision rate in the smallest bin, including nucleation and condensation
             !   see Mark Z. Jacobson, Fundamentals of Atmospheric Modeling, Second Edition (2005)
             !   equation (16.73)
             zcolrate_ocnv = zcolrate
             zcolrate_ocnv(in1a) = zcolrate_ocnv(in1a) + zj3n3(ii,jj,2)/pcocnv(ii,jj)

             zcs_ocnv = zcs_tot + zj3n3(ii,jj,2)/pcocnv(ii,jj)                   ! total sink for organic vapor

             zcvap_new2 = pcocnv(ii,jj)/(1._dp+ptstep*zcs_ocnv)                  ! new gas phase concentration [#/m3]
             zdvap2 = pcocnv(ii,jj) - zcvap_new2                                 ! change in gas concentration [#/m3]
             pcocnv(ii,jj) = zcvap_new2                                          ! updating vapour concentration [#/m3]
             
             zdvoloc = zcolrate_ocnv(in1a:fn2b)/zcs_ocnv*mvoc*zdvap2             ! volume change of particles 
                                                                                 !  [m3(OC)/m3(air)]

             paero(ii,jj,in1a:fn2b)%volc(2) = paero(ii,jj,in1a:fn2b)%volc(2) + & !-- change of volume
                                                    zdvoloc                      !   due to condensation in 1a-2b

             ! Condensation on hydromets 
             pcloud(ii,jj,1:ncld)%volc(2) = pcloud(ii,jj,1:ncld)%volc(2) +  &
                  zcolrateca(1:ncld)/zcs_ocnv*mvoc*zdvap2

             ! Condensation on rain drops
             pprecp(ii,jj,1:nprc)%volc(2) = pprecp(ii,jj,1:nprc)%volc(2) +  &
                  zcolratepa(1:nprc)/zcs_ocnv*mvoc*zdvap2

             !-- Change of number concentration in the smallest bin caused by nucleation
             !   Jacobson (2005), equation (16.75)
             ! If zxocnv = 0, then the chosen nucleation mechanism does not take into account
             ! the nonvolatile organic vapors and thus the pnaero does not have to be updated.
             IF (zxocnv(ii,jj) > 0._dp) THEN 
                paero(ii,jj,in1a)%numc = paero(ii,jj,in1a)%numc + &
                     zn_vs_c * zdvoloc(in1a)/mvoc/(n3*zxocnv(ii,jj))
             END IF
             
          END IF


          !---- 5.1.2) Semivolatile organic compound: regimes 1, 2 and 3
          zcs_ocsv = sum(zcolrate(in2a:fn2b)) +  &       ! sink for semivolatile organics
                     sum(zcolrateca(1:ncld))  +  &       ! ... including condensation on cloud droplets
                     sum(zcolratepa(1:nprc))             ! and rain drops 

          IF(pcocsv(ii,jj) > 1.e-10_dp .and. zcs_ocsv > 1.e-30_dp .AND. IsUsed(prtcl,'OC')) THEN

 
             zcvap_new3 = pcocsv(ii,jj)/(1._dp+ptstep*zcs_ocsv)   ! new gas phase concentration [#/m3]
             zdvap3 = pcocsv(ii,jj) - zcvap_new3                  ! change in gas concentration [#/m3]
             pcocsv(ii,jj) = zcvap_new3                           ! updating gas concentration [#/m3]
             
             zdvoloc(in2a:fn2b) = zdvoloc(in2a:fn2b) + &          ! volume change of particles 
                  zcolrate(in2a:fn2b)/zcs_ocsv*mvoc*zdvap3        !  [m3(OC)/m3(air)]

             paero(ii,jj,in1a:fn2b)%volc(2) = &                   !-- change of volume due
                  paero(ii,jj,in1a:fn2b)%volc(2) + zdvoloc        !   due to condensation in 1a-2b

             ! Condensation on hydromets
             pcloud(ii,jj,1:ncld)%volc(2) = pcloud(ii,jj,1:ncld)%volc(2)  +  &
                  zcolrateca(1:ncld)/zcs_ocsv*mvoc*zdvap3

             ! Condensation on rain drops
             pprecp(ii,jj,1:nprc)%volc(2) = pprecp(ii,jj,1:nprc)%volc(2)  +  &
                  zcolratepa(1:nprc)/zcs_ocsv*mvoc*zdvap3
             
          END IF


          ! ---- 5.2) Sulphate -------------------------------------------
          IF(pcsa(ii,jj) > 1.e-10_dp .and. zcs_tot > 1.e-30_dp .AND. IsUsed(prtcl,'SO4')) THEN 

             !-- Ratio of mass transfer between nucleation and condensation

             zn_vs_c = 0._dp

             IF(zj3n3(ii,jj,1) > 1._dp) zn_vs_c = (zj3n3(ii,jj,1)) / &
                                              (zj3n3(ii,jj,1) +  &
                                              pcsa(ii,jj) * zcolrate(in1a))

             !   collision rate in the smallest bin, including nucleation and condensation
             !   see Mark Z. Jacobson, Fundamentals of Atmospheric Modeling, Second Edition (2005)
             !   equation (16.73)
             zcolrate(in1a) = zcolrate(in1a) + zj3n3(ii,jj,1) / pcsa(ii,jj)

             zcs_su = zcs_tot + zj3n3(ii,jj,1) / pcsa(ii,jj)      ! total sink for sulfate

             !--- Sulphuric acid -------------------------
             !
             zcvap_new1 = pcsa(ii,jj) /(1.+ptstep*zcs_su)         ! new gas phase concentration [#/m3]
             zdvap1 = pcsa(ii,jj) - zcvap_new1                    ! change in gas concentration [#/m3]
             pcsa(ii,jj) = zcvap_new1                             ! updating vapour concentration [#/m3]
             
             zdvolsa = zcolrate(in1a:fn2b)/zcs_su*mvsu*zdvap1     ! volume change of particles
             ! [m3(SO4)/m3(air)] by condensation

             !-- Change of volume concentration of sulphate in aerosol [fxm]
             paero(ii,jj,in1a:fn2b)%volc(1) = paero(ii,jj,in1a:fn2b)%volc(1) + zdvolsa

             !-- Clouds
             pcloud(ii,jj,1:ncld)%volc(1) = pcloud(ii,jj,1:ncld)%volc(1)  +  &
                  zcolrateca(1:ncld)/zcs_su*mvsu*zdvap1

             ! Rain drops
             pprecp(ii,jj,1:nprc)%volc(1) = pprecp(ii,jj,1:nprc)%volc(1)  +  &
                  zcolratepa(1:nprc)/zcs_su*mvsu*zdvap1

             !-- Change of number concentration in the smallest bin caused by nucleation
             !   Jacobson (2005), equation (16.75)
             IF (zxsa(ii,jj) > 0._dp) THEN
                paero(ii,jj,in1a)%numc = paero(ii,jj,in1a)%numc +          &
                     zn_vs_c * zdvolsa(in1a)/mvsu/(n3*zxsa(ii,jj))
             END IF

          END IF
           
       END DO ! kproma

    END DO ! klev

    ! -- 5.3) Water vapour 
    CALL gpparth2o(kproma,kbdim,klev,krow,  &
                   paero, pcloud, pprecp,   &
                   ptemp,ppres,prs,prv,     &
                   ptstep)

    ! -- 5.4) HNO3/NH3
    !CALL gpparthno3(kproma,kbdim,klev,krow,ppres,ptemp,paero,pcloud,   &
    !                pprecp,pchno3,pcnh3,prv,prs,zbeta,ptstep           )



  END SUBROUTINE condensation

!
! ----------------------------------------------------------------------------------------------------------
!

  SUBROUTINE gpparth2o(kproma, kbdim,  klev, krow,  &
                       paero,  pcloud, pprecp,      &
                       ptemp,  ppres,  prs, prv,    &
                       ptstep)
    USE mo_kind, ONLY : dp
    USE mo_submctl, ONLY : t_section,            &
                               nbins, ncld, nprc,    &
                               rhowa, mwa, mair,     &
                               surfw0, rg,           &
                               pi, prlim, nlim,      &
                               massacc,avog,pstand,  &
                               in1a,fn1a,in2a,fn2a,  &
                               in2b,fn2b,            &
                               lscndh2oae, lscndh2ocl
    USE mo_constants, ONLY : alv
    USE mo_salsa_properties, ONLY : equilibration
    IMPLICIT NONE

    INTEGER, INTENT(in) :: kproma,kbdim,klev,krow
    REAL(dp), INTENT(in) :: ptstep
    REAL(dp), INTENT(in) :: ptemp(kbdim,klev), ppres(kbdim,klev), prs(kbdim,klev)
    TYPE(t_section), INTENT(inout) :: paero(kbdim,klev,nbins),  &
                                      pcloud(kbdim,klev,ncld),  &
                                      pprecp(kbdim,klev,nprc)
    REAL(dp), INTENT(inout) :: prv(kbdim,klev)
    
    REAL(dp) :: zkelvin(nbins), zkelvincd(ncld), zkelvinpd(nprc)  ! Kelvin effects
    REAL(dp) :: zka(nbins), zkacd(ncld), zkapd(nprc)            ! Activity coefficients
    REAL(dp) :: zcwsurfae(nbins), zcwsurfcd(ncld), zcwsurfpd(nprc) ! Surface mole concentrations
    REAL(dp) :: zmtae(nbins), zmtcd(ncld), zmtpd(nprc)        ! Mass transfer coefficients
    REAL(dp) :: zwsatae(nbins), zwsatcd(ncld), zwsatpd(nprc)  ! Water saturation ratios above
    REAL(dp) :: zmttot                                        ! Total condensation rate
    REAL(dp) :: zcwtot                                        ! Total water mole concentration
    REAL(dp) :: zcwc, zcwn, zcwint                            ! Current and new water vapour mole concentrations
    REAL(dp) :: zcwcae(nbins), zcwnae(nbins), zcwintae(nbins) ! Current and new water mole concentrations in aerosols
    REAL(dp) :: zcwccd(ncld), zcwncd(ncld), zcwintcd(ncld)    !     -  ''  -     in cloud droplets
    REAL(dp) :: zcwcpd(nprc), zcwnpd(nprc), zcwintpd(nprc)    !     -  ''  -     in rain drops
    REAL(dp) :: zdcwae(nbins), zdcwcd(ncld), zdcwpd(nprc)
    REAL(dp) :: zdfh2o, zthcond,rhoair
    REAL(dp) :: zbeta,zknud,zmfph2o
    REAL(dp) :: zact, zhlp1,zhlp2,zhlp3
    REAL(dp) :: adt,adtc(nbins),ttot
    REAL(dp) :: testi(nbins)
    REAL(dp) :: hdasha(nbins)                                 ! H' in Eq (17.104) in Jacobson (2005) for aerosol
    REAL(dp) :: hdashc(nbins)                                 !           "                       cloud droplets
    REAL(dp) :: hdashp(nbins)                                 !           "           precipitation hydrometeors

    
    REAL(dp) :: zrh(kbdim,klev)

    REAL(dp) :: zaelwc1(kbdim,klev), zaelwc2(kbdim,klev)

    INTEGER :: nstr
    INTEGER :: ii,jj,cc

    INTEGER :: poista

    zrh(:,:) = prv(:,:)/prs(:,:)
    
    ! Calculate the condensation only for 2a/2b aerosol bins
    nstr = in2a

    ! Save the current aerosol water content 
    zaelwc1(:,:) = SUM(paero(:,:,in1a:fn2b)%volc(8),DIM=3)*rhowa

    ! For 1a bins do the equilibrium calculation
    CALL equilibration(kproma,kbdim,klev,      &
                       zrh,ptemp,paero,.FALSE. )

    ! If RH < 98 % OR dynamic condensation for aerosols switched off, do equilibrium for all bins
    IF (zrh(1,1) < 0.98_dp .OR. .NOT. lscndh2oae)  CALL equilibration(kproma,kbdim,klev,      &
                                                                      zrh,ptemp,paero,.TRUE. )

    ! The new aerosol water content after equilibrium calculation
    zaelwc2(:,:) = SUM(paero(:,:,in1a:fn2b)%volc(8),DIM=3)*rhowa

    prv(:,:) = prv(:,:) - ( zaelwc2(:,:) - zaelwc1(:,:) )*ppres(:,:)*mair/(rg*ptemp(:,:))


    adtc(:) = 0._dp
    zcwc = 0._dp; zcwint = 0._dp; zcwn = 0._dp
    zcwcae = 0._dp; zcwccd = 0._dp; zcwcpd = 0._dp
    zcwintae = 0._dp; zcwintcd = 0._dp; zcwintpd = 0._dp
    zcwnae = 0._dp; zcwncd = 0._dp; zcwnpd = 0._dp

    DO jj = 1,klev
       DO ii = 1,kproma

          rhoair = mair*ppres(ii,jj)/(rg*ptemp(ii,jj))

          zdfh2o = ( 5._dp/(16._dp*avog*rhoair*1.e-3_dp*(3.11e-8_dp)**2) ) * &
                   SQRT( rg*1e7_dp*ptemp(ii,jj)*mair*1.e3_dp*(mwa+mair)*1.e3_dp/( 2._dp*pi*mwa*1.e3_dp ) )
          zdfh2o = zdfh2o*1.e-4

          zmfph2o = 3._dp*zdfh2o*sqrt(pi*mwa/(8._dp*rg*ptemp(ii,jj)))
          zthcond = 0.023807_dp + 7.1128e-5_dp*(ptemp(ii,jj) - 273.16_dp) ! Thermal conductivity of air 

          ! -- Water vapour (Follows the analytical predictor method by Jacobson 2005)
          zkelvinpd = 1._dp; zkelvincd = 1._dp; zkelvin = 1._dp
          zka = 1._dp; zkacd = 1._dp; zkapd = 1._dp ! Assume activity coefficients as 1 for now.

          ! Kelvin effects
          zkelvin(1:nbins) = exp( 4._dp*surfw0*mwa /  &
               (rg*ptemp(ii,jj)*rhowa*paero(ii,jj,1:nbins)%dwet) )

          zkelvincd(1:ncld) = exp( 4._dp*surfw0*mwa /  &
               (rg*ptemp(ii,jj)*rhowa*pcloud(ii,jj,1:ncld)%dwet) )
             
          zkelvinpd(1:nprc) = exp( 4._dp*surfw0*mwa /  & 
               (rg*ptemp(ii,jj)*rhowa*MIN(pprecp(ii,jj,1:nprc)%dwet,2.e-3_dp)) )
 
          ! Cloud droplets --------------------------------------------------------------------------------
          zmtcd(:) = 0._dp
          zcwsurfcd(:) = 0._dp
          DO cc = 1,ncld
             IF (pcloud(ii,jj,cc)%numc > prlim .AND. lscndh2ocl) THEN
          
                ! Activity + Kelvin effect
                zact = acth2o(pcloud(ii,jj,cc))

                ! Saturation mole concentration over flat surface
                zcwsurfcd(cc) = prs(ii,jj)*rhoair/mwa
                
                ! Equilibrium saturation ratio
                zwsatcd(cc) = zact*zkelvincd(cc)

                !-- transitional correction factor
                zknud = 2._dp*zmfph2o/pcloud(ii,jj,cc)%dwet   
                zbeta = (zknud + 1._dp)/(0.377_dp*zknud+1._dp+4._dp/ &    
                     (3._dp)*(zknud+zknud**2))  

                ! Mass transfer according to Jacobson
                zhlp1 = pcloud(ii,jj,cc)%numc*2._dp*pi*pcloud(ii,jj,cc)%dwet*zdfh2o*zbeta       ! Jacobson Eq (16.64) without D^eff 
                                                                                                !                     fully calculated
                zhlp2 = mwa*zbeta*zdfh2o*alv*zwsatcd(cc)*zcwsurfcd(cc)/(zthcond*ptemp(ii,jj))   !     "    Eq (16.55) 1st term
                                                                                                ! in the left side of the denominator   
                zhlp3 = ( (alv*mwa)/(rg*ptemp(ii,jj)) ) - 1._dp                                 !     "    Eq (16.55) 2nd term
                                                                                                ! in the left side of the denominator           
                zmtcd(cc) = zhlp1/( zhlp2*zhlp3 + 1._dp )!*(min( zrh(ii,jj)*pcloud(ii,jj,cc)%dwet/1.e-5,1.0))**2.0 ! " Eq (16.64) 
                                                                                                !mass transfer coefficient **

             END IF
          END DO

          ! Rain drops --------------------------------------------------------------------------------
          zmtpd(:) = 0._dp
          zcwsurfpd(:) = 0._dp
          DO cc = 1,nprc
             IF (pprecp(ii,jj,cc)%numc > prlim .AND. lscndh2ocl) THEN
                
                ! Activity + Kelvin effect
                zact = acth2o(pprecp(ii,jj,cc))
                
                ! Saturation mole concentrations over flat surface
                zcwsurfpd(cc) = prs(ii,jj)*rhoair/mwa
                   
                ! Equilibrium saturation ratio
                zwsatpd(cc) = zact*zkelvinpd(cc)

                !-- transitional correction factor
                zknud = 2._dp*zmfph2o/pprecp(ii,jj,cc)%dwet   
                zbeta = (zknud + 1._dp)/(0.377_dp*zknud+1._dp+4._dp/ &    
                     (3._dp)*(zknud+zknud**2))  

                ! Mass transfer according to Jacobson
                zhlp1 = pprecp(ii,jj,cc)%numc*2._dp*pi*pprecp(ii,jj,cc)%dwet*zdfh2o*zbeta 
                zhlp2 = mwa*zdfh2o*alv*zwsatpd(cc)*zcwsurfpd(cc)/(zthcond*ptemp(ii,jj))
                zhlp3 = ( (alv*mwa)/(rg*ptemp(ii,jj)) ) - 1._dp
          
                zmtpd(cc) = zhlp1/( zhlp2*zhlp3 + 1._dp )                                       ! see above (**)

             END IF
          END DO

          ! -- Aerosols: ------------------------------------------------------------------------------------
          zmtae(:) = 0._dp
          zcwsurfae(:) = 0._dp
          DO cc = 1,nbins
             IF (paero(ii,jj,cc)%numc > nlim .AND. zrh(ii,jj) > 0.98_dp .AND. lscndh2oae) THEN

                ! Water activity
                zact = acth2o(paero(ii,jj,cc))
                testi(cc) = zact
!                write(*,*) zact 
                ! Ssaturation mole concentration over flat surface
                ! Limit the supersaturation to max 1.01 for the mass transfer
                ! EXPERIMENTAL
                zcwsurfae(cc) = MAX(prs(ii,jj),prv(ii,jj)/1.01_dp)*rhoair/mwa

                ! Equilibrium saturation ratio
                zwsatae(cc) = zact*zkelvin(cc)

                !-- transitional correction factor
                zknud = 2._dp*zmfph2o/paero(ii,jj,cc)%dwet   
                zbeta = (zknud + 1._dp)/(0.377_dp*zknud+1._dp+4._dp/ &    
                     (3._dp*massacc(cc))*(zknud+zknud**2))  

                ! Mass transfer
                zhlp1 = paero(ii,jj,cc)%numc*2._dp*pi*paero(ii,jj,cc)%dwet*zdfh2o*zbeta
                zhlp2 = mwa*zdfh2o*alv*zwsatae(cc)*zcwsurfae(cc)/(zthcond*ptemp(ii,jj))
                zhlp3 = ( (alv*mwa)/(rg*ptemp(ii,jj)) ) - 1._dp
                
                zmtae(cc) = zhlp1/( zhlp2*zhlp3 + 1._dp )                                       ! see above (**) 
          
             END IF
          END DO

          ! UGLY FIX
          ! See above for possible replacement
          !IF ( zrh(ii,jj) > 1.01_dp) zmtae(:) = 0._dp

          ! Current mole concentrations
          zcwc = prv(ii,jj)*rhoair/mwa
          zcwcae(1:nbins) = paero(ii,jj,1:nbins)%volc(8)*rhowa/mwa                             
          zcwccd(1:ncld) = pcloud(ii,jj,1:ncld)%volc(8)*rhowa/mwa
          zcwcpd(1:nprc) = pprecp(ii,jj,1:nprc)%volc(8)*rhowa/mwa

          zcwtot = zcwc + SUM(zcwcae) + &                                                       ! Jacobson Eq (16.70)
                          SUM(zcwccd) + &
                          SUM(zcwcpd)
          ttot = 0._dp
          adtc = 0._dp

          zcwintae = zcwcae; zcwintcd = zcwccd; zcwintpd = zcwcpd
          zcwint = 0._dp
          poista = 0
          DO WHILE (ttot < ptstep)

             poista = poista + 1

             adt=2.e-2
             ! New vapor concentration
             zhlp1 = zcwc + adt * ( SUM(zmtae(nstr:nbins)*zwsatae(nstr:nbins)*zcwsurfae(nstr:nbins))  + & ! Jacobson Eq (16.71)
                                    SUM(zmtcd(1:ncld)*zwsatcd(1:ncld)*zcwsurfcd(1:ncld))              + & ! numerator
                                    SUM(zmtpd(1:nprc)*zwsatpd(1:nprc)*zcwsurfpd(1:nprc))                )
             zhlp2 = 1._dp + adt * ( SUM(zmtae(nstr:nbins)) + SUM(zmtcd(1:ncld)) + SUM(zmtpd(1:nprc)) )   ! denominator
             zcwint = zhlp1/zhlp2
             zcwint = MIN(zcwint,zcwtot)
             
             IF ( ANY(paero(ii,jj,:)%numc > nlim) .AND. zrh(ii,jj) > 0.98_dp ) THEN
                DO cc = nstr,nbins
                   zcwintae(cc) = zcwcae(cc) + min(max(adt*zmtae(cc)*(zcwint - zwsatae(cc)*zcwsurfae(cc)), &
                        -0.02*zcwcae(cc)),0.05*zcwcae(cc))  
                   zwsatae(cc) = acth2o(paero(ii,jj,cc),zcwintae(cc))*zkelvin(cc)
                END DO
             END IF
             IF ( ANY(pcloud(ii,jj,:)%numc > nlim) ) THEN
                DO cc = 1,ncld
                   zcwintcd(cc) = zcwccd(cc) + min(max(adt*zmtcd(cc)*(zcwint - zwsatcd(cc)*zcwsurfcd(cc)), &
                        -0.02*zcwccd(cc)),0.05*zcwccd(cc))
                   zwsatcd(cc) = acth2o(pcloud(ii,jj,cc),zcwintcd(cc))*zkelvincd(cc)
                END DO
             END IF
             IF ( ANY(pprecp(ii,jj,:)%numc > prlim) ) THEN
                DO cc = 1,nprc
                   zcwintpd(cc) = zcwcpd(cc) + min(max(adt*zmtpd(cc)*(zcwint - zwsatpd(cc)*zcwsurfpd(cc)), &
                        -0.02*zcwcpd(cc)),0.05*zcwcpd(cc))
                   zwsatpd(cc) = acth2o(pprecp(ii,jj,cc),zcwintpd(cc))*zkelvinpd(cc)
                END DO
             END IF

             zcwintae(nstr:nbins) = MAX(zcwintae(nstr:nbins),0._dp)
             zcwintcd(1:ncld) = MAX(zcwintcd(1:ncld),0._dp)
             zcwintpd(1:nprc) = MAX(zcwintpd(1:nprc),0._dp)

             ! Update saturation ratios
             !DO cc = nstr,nbins
             !   zwsatae(cc) = acth2o(paero(ii,jj,cc),zcwintae(cc))*zkelvin(cc)
             !END DO
             !DO cc = 1,ncld
             !   zwsatcd(cc) = acth2o(pcloud(ii,jj,cc),zcwintcd(cc))*zkelvincd(cc)
             !END DO
             !DO cc = 1,nprc
             !   zwsatpd(cc) = acth2o(pprecp(ii,jj,cc),zcwintpd(cc))*zkelvinpd(cc)
             !END DO

             ! Updae vapor concentration for consistency
             zcwint = zcwtot - SUM(zcwintae(1:nbins)) - &
                               SUM(zcwintcd(1:ncld))     - &
                               SUM(zcwintpd(1:nprc))
          
             ! Update "old" values for next cycle
             zcwcae = zcwintae; zcwccd = zcwintcd; zcwcpd = zcwintpd
             zcwc = zcwint

             ttot = ttot + adt
          
          END DO ! ADT

          zcwn = zcwint
          zcwnae = zcwintae
          zcwncd = zcwintcd
          zcwnpd = zcwintpd

          prv(ii,jj) = zcwn*mwa/rhoair
          
          paero(ii,jj,1:nbins)%volc(8) = zcwnae(1:nbins)*mwa/rhowa
          pcloud(ii,jj,1:ncld)%volc(8) = zcwncd(1:ncld)*mwa/rhowa
          pprecp(ii,jj,1:nprc)%volc(8) = zcwnpd(1:nprc)*mwa/rhowa

       END DO
       
    END DO




    
  END SUBROUTINE gpparth2o
  !-------------------------------------------------------
  REAL(dp) FUNCTION acth2o(ppart,pcw)
    USE mo_kind, ONLY : dp
    USE mo_submctl, ONLY : t_section,  &
                               rhosu, msu,   &
                               rhooc, moc,   &
                               rhoss, mss,   &
                               rhowa, mwa,   &
                               rhono, mno,   &
                               rhonh, mnh
    IMPLICIT NONE

    TYPE(t_section), INTENT(in) :: ppart
    REAL(dp), INTENT(in), OPTIONAL :: pcw
    !REAL(dp), INTENT(out) :: pact

    REAL(dp) :: zns, znw
    
    zns =  ( 3._dp*(ppart%volc(1)*rhosu/msu) +  &
                   (ppart%volc(2)*rhooc/moc) +  &
             2._dp*(ppart%volc(5)*rhoss/mss) +  &
                   (ppart%volc(6)*rhono/mno) +  &
                   (ppart%volc(7)*rhonh/mnh) ) 
    
    IF (PRESENT(pcw)) THEN
       znw = pcw
    ELSE
       znw = ppart%volc(8)*rhowa/mwa
    END IF
                
    ! Assume activity coefficient of 1 for water...
    acth2o = MAX(znw/(znw+zns),0.1_dp)
  END FUNCTION acth2o

!
! ----------------------------------------------------------------------------------------------------------
!

  SUBROUTINE gpparthno3(kproma,kbdim,klev,krow,ppres,ptemp,paero,pcloud,    &
                        pprecp,pghno3,pgnh3,prv,prs,pbeta,ptstep)
    USE mo_kind, ONLY : dp
    USE mo_submctl, ONLY : t_section,           &
                               nbins, ncld, nprc,   &
                               in1a, fn2b,          &
                               surfw0, mvno, mvnh, boltz, &
                               rhono, mno,          &
                               rhonh, mnh,          &
                               rhosu, msu,          &
                               avog, pi,            &
                               pstand,              &
                               nlim, prlim
                               
    IMPLICIT NONE

    ! OTA MOLAALISUUDET PDFITESTA JA LASKE BETAT MYÖS PISAROILLE!!!!

    
    INTEGER, INTENT(in) :: kproma,kbdim,klev,krow
    REAL(dp), INTENT(in) :: ptstep
    REAL(dp), INTENT(in) :: ptemp(kbdim,klev), ppres(kbdim,klev)
    REAL(dp), INTENT(in) :: prv(kbdim,klev),prs(kbdim,klev)
    REAL(dp), INTENT(in) :: pbeta(kbdim,klev,nbins)

    TYPE(t_section), INTENT(inout) :: paero(kbdim,klev,nbins),   &
                                      pcloud(kbdim,klev,ncld),   &
                                      pprecp(kbdim,klev,nprc)
    REAL(dp), INTENT(inout) :: pghno3(kbdim,klev),   &
                               pgnh3(kbdim,klev)

    REAL(dp) :: zkelno3ae(nbins), zkelno3cd(ncld), zkelno3pd(nprc)     ! Kelvin effects for HNO3
    REAL(dp) :: zkelnh3ae(nbins), zkelnh3cd(ncld), zkelnh3pd(nprc)     ! Kelvin effects for NH3
    REAL(dp) :: zkelw(nbins), zkelwcd(ncld), zkelwpd(nprc)           ! Kelvin effects for H2O

    REAL(dp) :: zcno3cae(nbins), zcno3intae(nbins), zcno3nae(nbins),  & ! Current, intermediate and new HNO3 in aerosols
                zcnh3cae(nbins), zcnh3intae(nbins), zcnh3nae(nbins),  & !  -  ''  - NH3

                zcno3ccd(ncld), zcno3intcd(ncld), zcno3ncd(ncld),  & ! -  ''  - HNO3 in cloud droplets
                zcnh3ccd(ncld), zcnh3intcd(ncld), zcnh3ncd(ncld),  & ! -  ''  - NH3 

                zcno3cpd(nprc), zcno3intpd(nprc), zcno3npd(nprc),  & ! -  ''  - HNO3 in precipitation
                zcnh3cpd(nprc), zcnh3intpd(nprc), zcnh3npd(nprc)     ! -  ''  - NH3

    REAL(dp) :: zcno3c, zcno3int, zcno3n                             ! Current, intermediate and new HNO3 gas concentration
    REAL(dp) :: zcnh3c, zcnh3int, zcnh3n                             ! -  ''  - NH3

    REAL(dp) :: zcno3eqae(nbins), zcno3eqcd(ncld), zcno3eqpd(nprc)   ! Equilibrium particle concentrations of HNO3 (isorropia)
    REAL(dp) :: zcnh3eqae(nbins), zcnh3eqcd(ncld), zcnh3eqpd(nprc)   ! Equilibrium particle concentrations of NH3 (isorropia)

    REAL(dp) :: zcgnh3eqae(nbins), zcgno3eqae(nbins), &              ! Equilibrium gas concentrations 
                zcgnh3eqcd(ncld), zcgno3eqcd(ncld), &
                zcgnh3eqpd(nprc), zcgno3eqpd(nprc)

    REAL(dp) :: zacno3ae(nbins), zacno3cd(ncld), zacno3pd(nprc)       ! Activity coefficients for HNO3
    REAL(dp) :: zacnh3ae(nbins), zacnh3cd(ncld), zacnh3pd(nprc)       ! Activity coefficients for NH3
    REAL(dp) :: zacnh4hso2ae(nbins), zacnh4hso2cd(ncld), zacnh4hso2pd(nprc)
    REAL(dp) :: zachhso4ae(nbins), zachhso4cd(ncld), zachhso4pd(nprc)
    REAL(dp) :: zmolsae(nbins,7),zmolscd(ncld,7),zmolspd(nprc,7)  ! Ion molalities from pdfite
     
    REAL(dp) :: zcnh3tot, zcno3tot                                                 ! Total mole concentrations

    REAL(dp) :: zmtno3ae(nbins), zmtno3cd(ncld), zmtno3pd(nprc) ! Mass transfer coefficients for HNO3
    REAL(dp) :: zmtnh3ae(nbins), zmtnh3cd(ncld), zmtnh3pd(nprc) ! Mass transfer coefficients for NH3

    REAL(dp) :: zsathno3ae(nbins), zsathno3cd(ncld), zsathno3pd(nprc)
    REAL(dp) :: zsatnh3ae(nbins), zsatnh3cd(ncld), zsatnh3pd(nprc)

    REAL(dp) :: Hhno3ae(nbins),Hnh3ae(nbins)

    REAL(dp) :: zbeta ! LASKE
    REAL(dp) :: zdfvap ! Diffusion coefficient for vapors

    REAL(dp) :: zrh

    REAL(dp) :: zhlp1,zhlp2,zhlp3,zhlp4,zhlp5,zhlp6,zhlp7,zhlp8

    REAL(dp) :: adt,adtcae(2,nbins),adtcae2(2,nbins),adtccd(2,ncld),adtcpd(2,nprc) ! Adaptive timestep
    REAL(dp) :: adtc2(2,nbins)
    REAL(dp) :: telp,ttot ! Elapsed time

    INTEGER :: nstr
    INTEGER :: ii,jj,cc,tt

    nstr = 1
    ! ALUSTA KRIITTISET NOLLIKSI
    adtcae(:,:) = 0._dp
    adtccd(:,:) = 0._dp
    adtcpd(:,:) = 0._dp
    adtc2(:,:) = 0._dp

    DO jj = 1,klev
       DO ii = 1,kproma

          zkelno3pd = 1._dp; zkelno3cd = 1._dp; zkelno3ae = 1._dp
          zkelnh3pd = 1._dp; zkelnh3cd = 1._dp; zkelnh3ae = 1._dp
          zacno3ae = 1._dp; zacno3cd = 1._dp; zacno3pd = 1._dp
          zacnh3ae = 1._dp; zacnh3cd = 1._dp; zacnh3pd = 1._dp
          zsathno3ae = 1._dp; zsathno3cd = 1._dp; zsathno3pd = 1._dp
          zsatnh3ae = 1._dp; zsatnh3cd = 1._dp; zsatnh3pd = 1._dp

          zdfvap = 5.1111e-10_dp*ptemp(ii,jj)**1.75_dp*pstand/ppres(ii,jj)                ! diffusion coefficient [m2/s]

          ! Kelvin effects
          zkelno3ae(1:nbins) = exp( 4._dp*surfw0*mvno /  &
                                    (boltz*ptemp(ii,jj)*paero(ii,jj,1:nbins)%dwet) ) 
          zkelnh3ae(1:nbins) = exp( 4._dp*surfw0*mvnh /  &
                                    (boltz*ptemp(ii,jj)*paero(ii,jj,1:nbins)%dwet) )

          zkelno3cd(1:ncld) = exp( 4._dp*surfw0*mvno /  & 
                                   (boltz*ptemp(ii,jj)*pcloud(ii,jj,1:ncld)%dwet) )
          zkelnh3cd(1:ncld) = exp( 4._dp*surfw0*mvnh /  &
                                   (boltz*ptemp(ii,jj)*pcloud(ii,jj,1:ncld)%dwet) )

          zkelno3pd(1:nprc) = exp( 4._dp*surfw0*mvno /  &
                                   (boltz*ptemp(ii,jj)*pprecp(ii,jj,1:nprc)%dwet) )
          zkelnh3pd(1:nprc) = exp( 4._dp*surfw0*mvnh /  &
                                   (boltz*ptemp(ii,jj)*pprecp(ii,jj,1:nprc)%dwet) )

          ! Current gas concentrations
          zcno3c = pghno3(ii,jj)/avog
          zcnh3c = pgnh3(ii,jj)/avog

          ! Current particle concentrations
          zcno3cae(1:nbins) = paero(ii,jj,1:nbins)%volc(6)*rhono/mno
          zcnh3cae(1:nbins) = paero(ii,jj,1:nbins)%volc(7)*rhonh/mnh

          zcno3ccd(1:ncld) = pcloud(ii,jj,1:ncld)%volc(6)*rhono/mno
          zcnh3ccd(1:ncld) = pcloud(ii,jj,1:ncld)%volc(7)*rhonh/mnh

          zcno3cpd(1:nprc) = pprecp(ii,jj,1:nprc)%volc(6)*rhono/mno
          zcnh3cpd(1:nprc) = pprecp(ii,jj,1:nprc)%volc(7)*rhonh/mnh

          ! Total concentrations
          zcno3tot = zcno3c + SUM(zcno3cae(1:nbins)) +   &
                              SUM(zcno3ccd(1:ncld))       +   &
                              SUM(zcno3cpd(1:nprc))

          zcnh3tot = zcnh3c + SUM(zcnh3cae(1:nbins)) +   &
                              SUM(zcnh3ccd(1:ncld))       +   &
                              SUM(zcnh3cpd(1:nprc))

          ! BETA PUUTTUU PILVILTÄ
          zbeta = 1._dp

          ! Mass transfer coefficients 
          zmtno3ae = 0._dp; zmtnh3ae = 0._dp
          zmtno3cd = 0._dp; zmtnh3cd = 0._dp
          zmtno3pd = 0._dp; zmtnh3pd = 0._dp
          zmtno3ae(1:nbins) = 2._dp*pi*paero(ii,jj,1:nbins)%dwet *  &
                              zdfvap*paero(ii,jj,1:nbins)%numc*pbeta(ii,jj,1:nbins)
          zmtnh3ae(1:nbins) = 2._dp*pi*paero(ii,jj,1:nbins)%dwet *  &
                              zdfvap*paero(ii,jj,1:nbins)%numc*pbeta(ii,jj,1:nbins)

          zmtno3cd(1:ncld) = 2._dp*pi*pcloud(ii,jj,1:ncld)%dwet *  &
                             zdfvap*pcloud(ii,jj,1:ncld)%numc*zbeta
          zmtnh3cd(1:ncld) = 2._dp*pi*pcloud(ii,jj,1:ncld)%dwet *  &
                             zdfvap*pcloud(ii,jj,1:ncld)%numc*zbeta
          
          zmtno3pd(1:nprc) = 2._dp*pi*pprecp(ii,jj,1:nprc)%dwet *  &
                             zdfvap*pprecp(ii,jj,1:nprc)%numc*zbeta
          zmtnh3pd(1:nprc) = 2._dp*pi*pprecp(ii,jj,1:nprc)%dwet *  &
                             zdfvap*pprecp(ii,jj,1:nprc)%numc*zbeta

          zrh = prv(ii,jj)/prs(ii,jj)

          zmolsae = 0._dp
          zmolscd = 0._dp
          zmolspd = 0._dp

          ! Get the equilibrium concentrations
          ! Aerosols
          CALL NONHEquil(nbins,zrh,ptemp(ii,jj),paero(ii,jj,:),    &
                         zcgno3eqae,zcgnh3eqae,zacno3ae,zacnh3ae,  &
                         zacnh4hso2ae,zachhso4ae,zmolsae           )

          ! Cloud droplets
          CALL NONHEquil(ncld,zrh,ptemp(ii,jj),pcloud(ii,jj,:),    &
                         zcgno3eqcd,zcgnh3eqcd,zacno3cd,zacnh3cd,  &
                         zacnh4hso2cd,zachhso4cd,zmolscd           )

          ! Precipitation
          CALL NONHEquil(nprc,zrh,ptemp(ii,jj),pprecp(ii,jj,:),    &
                         zcgno3eqpd,zcgnh3eqpd,zacno3pd,zacnh3pd,  &
                         zacnh4hso2pd,zachhso4pd,zmolspd           )


          !WRITE(*,*) zmolsae(1:7,1)

          zsatnh3ae = 1._dp; zsathno3ae = 1._dp
          zsatnh3cd = 1._dp; zsathno3cd = 1._dp
          zsatnh3pd = 1._dp; zsathno3pd = 1._dp
          ! NH4/HNO3 saturation ratios for
          ! aerosols
          CALL SVsat(nbins,ptemp(ii,jj),paero(ii,jj,:),zacno3ae,zacnh3ae,zacnh4hso2ae,  &
               zachhso4ae,zcgno3eqae,zcno3cae,zcnh3cae,zkelno3ae,zkelnh3ae,       &
               zsathno3ae,zsatnh3ae,zmolsae,nlim                                          )
          ! clouds
          CALL SVsat(ncld,ptemp(ii,jj),pcloud(ii,jj,:),zacno3cd,zacnh3cd,zacnh4hso2cd,  &
               zachhso4cd,zcgno3eqcd,zcno3ccd,zcnh3ccd,zkelno3cd,zkelnh3cd,       &
               zsathno3cd,zsatnh3cd,zmolscd,nlim                                          )
          ! precipitation
          CALL SVsat(nprc,ptemp(ii,jj),pprecp(ii,jj,:),zacno3pd,zacnh3pd,zacnh4hso2pd,  &
               zachhso4pd,zcgno3eqpd,zcno3cpd,zcnh3cpd,zkelno3pd,zkelnh3pd,       &
               zsathno3pd,zsatnh3pd,zmolspd,prlim                                         )

          adt = ptstep

          !WRITE(*,*) adt

          zhlp1 = SUM( zcno3cae(nstr:nbins) /  &
               (1._dp + adt*zmtno3ae(nstr:nbins)*zsathno3ae(nstr:nbins)) )
          zhlp2 = SUM( zcno3ccd(1:ncld)     /  &
               (1._dp + adt*zmtno3cd(1:ncld)*zsathno3cd(1:ncld)) )
          zhlp3 = SUM( zcno3cpd(1:nprc)     /  &
               (1._dp + adt*zmtno3pd(1:nprc)*zsathno3pd(1:nprc)) )
          zhlp4 = SUM( zmtno3ae(nstr:nbins)/(1._dp + adt*zmtno3ae(nstr:nbins)*zsathno3ae(nstr:nbins)) )
          zhlp5 = SUM( zmtno3cd(1:ncld)/(1._dp + adt*zmtno3cd(1:ncld)*zsathno3cd(1:ncld)) )
          zhlp6 = SUM( zmtno3pd(1:nprc)/(1._dp + adt*zmtno3pd(1:nprc)*zsathno3pd(1:nprc)) )
          zcno3int = ( zcno3tot - (zhlp1 + zhlp2 + zhlp3) )/( 1._dp + adt*(zhlp4 + zhlp5 + zhlp6) )
            
          zhlp1 = SUM( zcnh3cae(nstr:nbins) / &
               (1._dp + adt*zmtnh3ae(nstr:nbins)*zsatnh3ae(nstr:nbins)) )
          zhlp2 = SUM( zcnh3ccd(1:ncld)     / &
               (1._dp + adt*zmtnh3cd(1:ncld)*zsatnh3cd(1:ncld)) )
          zhlp3 = SUM( zcnh3cpd(1:nprc)     / &
               (1._dp + adt*zmtnh3pd(1:nprc)*zsatnh3pd(1:nprc)) )
          zhlp4 = SUM( zmtnh3ae(nstr:nbins)/(1._dp + adt*zmtnh3ae(nstr:nbins)*zsatnh3ae(nstr:nbins)) )
          zhlp5 = SUM( zmtnh3cd(1:ncld)/(1._dp + adt*zmtnh3cd(1:ncld)*zsatnh3cd(1:ncld))  )
          zhlp6 = SUM( zmtnh3pd(1:nprc)/(1._dp + adt*zmtnh3pd(1:nprc)*zsatnh3pd(1:nprc))  )
          zcnh3int = ( zcnh3tot - (zhlp1 + zhlp2 + zhlp3) )/( 1._dp + adt*(zhlp4 + zhlp5 + zhlp6) )
                
          zcno3int = MIN(zcno3int, zcno3tot)
          zcnh3int = MIN(zcnh3int, zcnh3tot)
                
          ! Calculate the new particle concentrations
          zcno3intae(:) = zcno3cae(:)
          zcno3intcd(:) = zcno3ccd(:)
          zcno3intpd(:) = zcno3cpd(:)
          zcnh3intae(:) = zcnh3cae(:)
          zcnh3intcd(:) = zcnh3ccd(:)
          zcnh3intpd(:) = zcnh3cpd(:)
          DO cc = nstr,nbins
             zcno3intae(cc) = ( zcno3cae(cc) + adt*zmtno3ae(cc)*zcno3int ) /  &
                  ( 1._dp + adt*zmtno3ae(cc)*zsathno3ae(cc) )
             zcnh3intae(cc) = ( zcnh3cae(cc) + adt*zmtnh3ae(cc)*zcnh3int ) /  &
                  ( 1._dp + adt*zmtnh3ae(cc)*zsatnh3ae(cc) )
          END DO
          DO cc = 1,ncld
             zcno3intcd(cc) = ( zcno3ccd(cc) + adt*zmtno3cd(cc)*zcno3int ) /  &
                  ( 1._dp + adt*zmtno3cd(cc)*zsathno3cd(cc) )
             zcnh3intcd(cc) = ( zcnh3ccd(cc) + adt*zmtnh3cd(cc)*zcnh3int ) /  &
                  ( 1._dp + adt*zmtnh3cd(cc)*zsatnh3cd(cc) )
          END DO
          DO cc = 1,nprc
             zcno3intpd(cc) = ( zcno3cpd(cc) + adt*zmtno3pd(cc)*zcno3int ) /  &
                  ( 1._dp + adt*zmtno3pd(cc)*zsathno3pd(cc) )
             zcnh3intpd(cc) = ( zcnh3cpd(cc) + adt*zmtnh3pd(cc)*zcnh3int ) /  &
                  ( 1._dp + adt*zmtnh3pd(cc)*zsatnh3pd(cc) )
          END DO
          zcno3intae(1:nbins) = MAX(zcno3intae(1:nbins),0._dp)
          zcno3intcd(1:ncld) = MAX(zcno3intcd(1:ncld),0._dp)
          zcno3intpd(1:nprc) = MAX(zcno3intpd(1:nprc),0._dp)
          zcnh3intae(1:nbins) = MAX(zcnh3intae(1:nbins),0._dp)
          zcnh3intcd(1:ncld) = MAX(zcnh3intcd(1:ncld),0._dp)
          zcnh3intpd(1:nprc) = MAX(zcnh3intpd(1:nprc),0._dp)
                
          ! JACOBSONIN RAJOTUS TÄHÄ?

          !zcno3int = zcno3tot - SUM(zcno3intae(1:nbins)) - &
          !     SUM(zcno3intcd(1:ncld))  - &
          !     SUM(zcno3intpd(1:nprc))
                
          !zcnh3int = zcnh3tot - SUM(zcnh3intae(1:nbins)) - &
          !     SUM(zcnh3intcd(1:ncld))  - &
          !     SUM(zcnh3intpd(1:nprc))
          
          !WRITE(*,*) zcnh3int, SUM(zcnh3intae(1:nbins))
          !WRITE(*,*) zsatnh3ae(1:7)

          zcno3n = zcno3int
          zcno3nae(:) = zcno3intae(:); zcno3ncd(:) = zcno3intcd(:); zcno3npd(:) = zcno3intpd(:)
          zcnh3n = zcnh3int
          zcnh3nae(:) = zcnh3intae(:); zcnh3ncd(:) = zcnh3intcd(:); zcnh3npd(:) = zcnh3intpd(:)

          ! Model timestep reached - update the new arrays
          pghno3(ii,jj) = zcno3n*avog
          pgnh3(ii,jj) = zcnh3n*avog
          
          paero(ii,jj,1:nbins)%volc(6) = zcno3nae(1:nbins)*mno/rhono
          pcloud(ii,jj,1:ncld)%volc(6) = zcno3ncd(1:ncld)*mno/rhono
          pprecp(ii,jj,1:nprc)%volc(6) = zcno3npd(1:nprc)*mno/rhono

          paero(ii,jj,1:nbins)%volc(7) = zcnh3nae(1:nbins)*mnh/rhonh
          pcloud(ii,jj,1:ncld)%volc(7) = zcnh3ncd(1:ncld)*mnh/rhonh
          pprecp(ii,jj,1:nprc)%volc(7) = zcnh3npd(1:nprc)*mnh/rhonh

       END DO

    END DO

  END SUBROUTINE gpparthno3

  ! ---------------------------------------------------------------

  REAL(dp) FUNCTION acthno3(ppart,pgamma,pchno3p)
    USE mo_kind, ONLY : dp
    USE mo_submctl, ONLY : t_section,  &
                               rhosu, msu,   &
                               rhooc, moc,   &
                               rhobc, mbc,   &
                               rhodu, mdu,   &
                               rhoss, mss,   &
                               rhono, mno,   &
                               rhonh, mnh,   &
                               rhowa, mwa

    IMPLICIT NONE 
    
    TYPE(t_section), INTENT(in) :: ppart
    REAL(dp), INTENT(in), OPTIONAL :: pchno3p  ! Current particle HNO3 mole concentration
    REAL(dp), INTENT(in) :: pgamma

    REAL(dp) :: zns,znhno3

    ! Solute ~everything else (soluble)?
    zns = ( 3._dp*(ppart%volc(1)*rhosu/msu)  + &
                  (ppart%volc(2)*rhooc/moc)  + &
            2._dp*(ppart%volc(5)*rhoss/mss)  + &
                  (ppart%volc(7)*rhonh/mnh)  + &
                  (ppart%volc(8)*rhowa/mwa))

    IF (PRESENT(pchno3p)) THEN
       znhno3 = pchno3p
    ELSE
       znhno3 = ppart%volc(6)*rhono/mno
    END IF
    
    acthno3 = MAX(znhno3/(zns+znhno3),1.e-3_dp)*pgamma
  END FUNCTION acthno3
  ! -------------------------------------------------------
  REAL(dp) FUNCTION actnh3(ppart,pgamma,pcnh3p)
    USE mo_kind, ONLY : dp
    USE mo_submctl, ONLY : t_section,  &
                               rhosu, msu,   &
                               rhooc, moc,   &
                               rhobc, mbc,   &
                               rhodu, mdu,   &
                               rhoss, mss,   &
                               rhono, mno,   &
                               rhonh, mnh,   &
                               rhowa, mwa

    IMPLICIT NONE 
    
    TYPE(t_section), INTENT(in) :: ppart
    REAL(dp), INTENT(in), OPTIONAL :: pcnh3p  ! Current particle HNO3 mole concentration
    REAL(dp), INTENT(in) :: pgamma

    REAL(dp) :: zns,znnh3

    ! Solute ~everything else (soluble)?
    zns = ( 3._dp*(ppart%volc(1)*rhosu/msu)  + &
                  (ppart%volc(2)*rhooc/moc)  + &
            2._dp*(ppart%volc(5)*rhoss/mss)  + &
                  (ppart%volc(6)*rhono/mno)  + &
                  (ppart%volc(8)*rhowa/mwa))

    IF (PRESENT(pcnh3p)) THEN
       znnh3 = pcnh3p
    ELSE
       znnh3 = ppart%volc(7)*rhonh/mnh
    END IF
    
    actnh3 = MAX(znnh3/(zns+znnh3),1.e-3_dp)*pgamma
  END FUNCTION actnh3


  ! ------------------------------------------------------------------
  SUBROUTINE NONHEquil(nb,prh,ptemp,ppart,pcgno3eq,pcgnh3eq,      &
                       pgammano,pgammanh,pgammanh4hso2,pgammahhso4,pmols)
    USE mo_kind, ONLY : dp
    USE mo_submctl, ONLY : t_section,    &
                               rhosu,msu,    &
                               rhoss,mss,    &
                               rhono,mno,    &
                               rhonh,mnh,    &
                               rhowa,mwa,    &
                               rg,nlim
    USE aerosol_thermodynamics, ONLY : inorganic_pdfite
    IMPLICIT NONE
    
    INTEGER, INTENT(in) :: nb
    TYPE(t_section), INTENT(in) :: ppart(nb)
    REAL(dp), INTENT(in) :: prh,ptemp

    REAL(dp), INTENT(out) :: pgammano(nb),pgammanh(nb),pgammanh4hso2(nb),pgammahhso4(nb)
    REAL(dp), INTENT(out) :: pcgno3eq(nb),pcgnh3eq(nb)
    REAL(dp), INTENT(out) :: pmols(nb,7)

    REAL(dp) :: zions(7)                ! mol/m3

    REAL(dp) :: zwatertotal,        &   ! Total water in particles (mol/m3) ???
                zphno3,zphcl,zpnh3, &   ! Equilibrium vapor pressures (Pa??) 
                zgammas(7)              ! Activity coefficients
    
    REAL(dp) :: zhlp

    INTEGER :: cc

    pmols = 0._dp

    DO cc = 1,nb

       ! 2*H2SO4 + CL + NO3 - Na - NH4

       zhlp = 2._dp*ppart(cc)%volc(1)*rhosu/msu + ppart(cc)%volc(5)*rhoss/mss + &
              ppart(cc)%volc(6)*rhono/mno - ppart(cc)%volc(5)*rhoss/mss -  &
              ppart(cc)%volc(7)*rhonh/mnh

       zhlp = MAX(zhlp, 1.e-30_dp)

       zions(1) = zhlp !0._dp ! H+
       zions(2) = ppart(cc)%volc(7)*rhonh/mnh ! NH4
       zions(3) = ppart(cc)%volc(5)*rhoss/mss ! Na
       zions(4) = ppart(cc)%volc(1)*rhosu/msu ! SO4
       zions(5) = 0._dp ! HSO4
       zions(6) = ppart(cc)%volc(6)*rhono/mno ! NO3
       zions(7) = ppart(cc)%volc(5)*rhoss/mss ! Cl
       
       zwatertotal = ppart(cc)%volc(8)*rhowa/mwa

       CALL inorganic_pdfite(prh,ptemp,zions,zwatertotal,zphno3,zphcl,zpnh3,zgammas,pmols(cc,:))

       pgammano(cc) = zgammas(1)
       pgammanh(cc) = zgammas(3)
       pgammanh4hso2(cc) = zgammas(6)
       pgammahhso4(cc) = zgammas(7)

       !IF (ppart(cc)%numc < nlim) pmols(cc,:) = 0._dp

       !WRITE(*,*) pgammano(cc), pgammanh(cc), pgammanh4hso2(cc), pgammahhso4(cc)

       pcgno3eq(cc) = zphno3/(rg*ptemp)
       pcgnh3eq(cc) = zpnh3/(rg*ptemp)

    END DO

  END SUBROUTINE NONHEquil

  ! ------------------------------------------------------------------

  SUBROUTINE SVsat(nb,ptemp,ppart,pachno3,pacnh3,pacnh4hso2,   &
                   pachhso4,pchno3eq,pchno3,pcnh3,pkelhno3,    &
                   pkelnh3,psathno3,psatnh3,pmols,plim         )
    USE mo_kind, ONLY : dp
    USE mo_submctl, ONLY : t_section,   &
                               rhosu,msu,   &
                               rhooc,moc,   &
                               rhoss,mss,   &
                               rhono,mno,   &
                               rhonh,mnh,   &
                               rhowa,mwa,   &
                               rg, nlim
    IMPLICIT NONE
    
    ! Calculates the saturation ratio for semivolatile species
    
    INTEGER, INTENT(in) :: nb
    REAL(dp), INTENT(in) :: ptemp
    TYPE(t_section), INTENT(in) :: ppart(nb)
    REAL(dp), INTENT(in) :: pachno3(nb),pacnh3(nb),pacnh4hso2(nb),pachhso4(nb) ! Activity coefficients
    REAL(dp), INTENT(in) :: pchno3eq(nb) ! Equolibrium surface concentration of HNO3
    REAL(dp), INTENT(in) :: pchno3(nb)   ! Current particle concentration of HNO3
    REAL(dp), INTENT(in) :: pcnh3(nb)    ! Current particle concentration of NH3
    REAL(dp), INTENT(in) :: pkelhno3(nb), pkelnh3(nb)  ! Kelvin effects
    REAL(dp), INTENT(in) :: pmols(nb,7)
    REAL(dp), INTENT(in) :: plim
    REAL(dp), INTENT(out) :: psathno3(nb), psatnh3(nb)

    REAL(dp) :: KllH2O, KllNH3, KglNH3, KglHNO3

    REAL(dp), PARAMETER :: zt0 = 298.15_dp    ! Reference temp
    REAL(dp), PARAMETER :: zatm = 101325._dp  ! Unit atm in Pa
    
    REAL(dp), PARAMETER :: K01 = 1.01e-14_dp,   &
                           K02 = 1.81e-5_dp,    &
                           K03 = 57.64_dp,      &
                           K04 = 2.51e6_dp

    REAL(dp), PARAMETER :: a1 = -22.52_dp,      &
                           a2 = 1.50_dp,        &
                           a3 = 13.79_dp,       &
                           a4 = 29.47_dp,       &
                           b1 = 26.92_dp,       &
                           b2 = 26.92_dp,       &
                           b3 = -5.39_dp,       &
                           b4 = 16.84_dp
    
    REAL(dp) :: zmolno3     ! molality of NO3-
    REAL(dp) :: zmolhp      ! molality of H+
    REAL(dp) :: zmolso4,  & ! molality of SO4
                zmolcl,   & ! molality of Cl
                zmolnh4,  & ! Molality of NH4
                zmolna      ! Molality of Na
    REAL(dp) :: zhlp1,zhlp2,zhlp3, zxi

    INTEGER :: cc

    zhlp1 = 0._dp; zhlp2 = 0._dp; zhlp3 = 0._dp

    zmolno3 = 0._dp
    zmolhp = 0._dp

    ! Calculates K^ll_h20, K^ll_NH3, K^gl_NH3, K^gl_HNO3
    zhlp1 = zt0/ptemp
    zhlp2 = zhlp1-1._dp
    zhlp3 = 1._dp + LOG(zhlp1) - zhlp1

    KllH2O = K01*EXP( a1*zhlp2 + b1*zhlp3 )
    KllNH3 = K02*EXP( a2*zhlp2 + b2*zhlp3 )
    KglNH3 = k03*EXP( a3*zhlp2 + b3*zhlp3 )
    KglHNO3 = k04*EXP( a4*zhlp2 + b4*zhlp3 )

    !WRITE(*,*) 'testi', KllH2O, KllNH3, KglNH3, KglHNO3

    ! Get NO3- and H+ molality
    DO cc = 1,nb

       IF ( ppart(cc)%numc > 1.e-40_dp) THEN 

          zhlp1 = pcnh3(cc)*mnh + ppart(cc)%volc(1)*rhosu + ppart(cc)%volc(2)*rhooc +  &
               ppart(cc)%volc(5)*rhoss + ppart(cc)%volc(8)*rhowa
          zmolno3 = pchno3(cc) / zhlp1   

          zxi = ( pcnh3(cc) + ppart(cc)%volc(5)*rhoss/mss ) / &
                ( ppart(cc)%volc(1)*rhosu/msu )
          IF ( zxi <= 2._dp ) THEN
             !psathno3(cc) = acthno3(ppart(cc),pachno3(cc),pchno3(cc)) * &
             !        pkelhno3(cc) / KglHNO3
             !psatnh3(cc) = actnh3(ppart(cc),pacnh3(cc),pcnh3(cc)) *   &
             !     pkelnh3(cc) / KglNH3
             ! Molality of SO4
             zhlp1 = pcnh3(cc)*mnh + pchno3(cc)*mno + ppart(cc)%volc(2)*rhooc + &
                  ppart(cc)%volc(5)*rhoss + ppart(cc)%volc(8)*rhowa
             zmolso4 = (ppart(cc)%volc(1)*rhosu/msu)/zhlp1
             ! Molality of Cl
             zhlp1 = pcnh3(cc)*mnh + pchno3(cc)*mno + ppart(cc)%volc(2)*rhooc + &
                  ppart(cc)%volc(1)*rhosu + ppart(cc)%volc(8)*rhowa
             zmolcl = (ppart(cc)%volc(5)*rhoss/mss)/zhlp1
             ! Molality of NH4
             zhlp1 =  pchno3(cc)*mno + ppart(cc)%volc(1)*rhosu + ppart(cc)%volc(2)*rhooc + &
                  ppart(cc)%volc(5)*rhoss + ppart(cc)%volc(8)*rhowa
             zmolnh4 = pcnh3(cc)/zhlp1
             ! Molality of Na
             zhlp1 = pcnh3(cc)*mnh + pchno3(cc)*mno + ppart(cc)%volc(2)*rhooc + &
                  ppart(cc)%volc(1)*rhosu + ppart(cc)%volc(8)*rhowa
             zmolna = (ppart(cc)%volc(5)*rhoss/mss)/zhlp1

             !zmolhp = 2._dp*pmols(cc,4) + pmols(cc,5) + pmols(cc,6) + pmols(cc,7) - &
             !     (pmols(cc,2) + pmols(cc,3))!2._dp*zmolso4 + zmolno3 + zmolcl - (zmolnh4 + zmolna)

             zmolhp = 2._dp*zmolso4 + zmolno3 + zmolcl - (zmolnh4 + zmolna)
             
             !zmolhp = 2._dp*pmols(cc,4) + pmols(cc,5) + zmolno3 + zmolcl - (zmolnh4 + zmolna)
             
             !WRITE(*,*) zmolso4, pmols(cc,4), pmols(cc,5),cc

          ELSE

             !WRITE(*,*) 'HEP'
             zhlp2 = pkelhno3(cc)*zmolno3*pachno3(cc)**2._dp
             zmolhp = KglHNO3*pchno3eq(cc)/zhlp2

          END IF
          
          zhlp1 = ppart(cc)%volc(8)*rhowa*rg*ptemp*KglHNO3
          ! Saturation ratio for NH3
          zhlp2 = pkelnh3(cc)/(zhlp1*zmolhp)
          zhlp3 = KllH2O/(KllNH3+KglNH3)
          psatnh3(cc) = zhlp2 * ( (pacnh4hso2(cc)/pachhso4(cc))**2._dp ) * zhlp3
          ! Saturation ratio for HNO3
          psathno3(cc) = ( pkelhno3(cc)*zmolhp*pachno3(cc)**2 ) / zhlp1

       ELSE
          
          ! Ei hyvä näin!!!!
          psatnh3(cc) = 0._dp
          psathno3(cc) = 0._dp
       END IF

    END DO

  END SUBROUTINE SVsat


  ! ------------------------------------------------------------------

  FUNCTION GetTstep(nb,zcg,zcs,zmt,zconst) RESULT(tscale)
    USE mo_kind, ONLY : dp
    IMPLICIT NONE
    
    INTEGER, INTENT(in) :: nb
    REAL(dp), INTENT(in) :: zcg
    REAL(dp), INTENT(in) :: zcs(nb), zmt(nb)
    REAL(dp), INTENT(in) :: zconst
    
    REAL(dp) :: th(nb)
    REAL(dp) :: tscale

    INTEGER :: cc
    
    DO cc = 1,nb
       th(cc) = (zcg - zcs(cc))/MAX(zcg,zcs(cc))   
    END DO

    tscale = zconst/SUM( ABS(zmt(:)*th(:)),MASK=(th(:)*zmt(:) /= 0._dp) )

  END FUNCTION GetTstep

!
! ----------------------------------------------------------------------------------------------------------
!

  FUNCTION satvaph2o(ptemp) RESULT(psat)
    !-----------------------------------------------------------------
    ! Saturation vapour pressure of water vapour
    ! This is a local function for the subroutine *cloud_condensation*
    !
    ! J. Tonttila, FMI, 03/2014
    !-----------------------------------------------------------------
    USE mo_kind, ONLY : dp
    IMPLICIT NONE
    
    REAL(dp), INTENT(in) :: ptemp
    
    REAL(dp), PARAMETER ::        & 
         za0 = 6.107799961_dp,    &
         za1 = 4.436518521e-1_dp, &
         za2 = 1.428945805e-2_dp, &
         za3 = 2.650648471e-4_dp, &
         za4 = 3.031240396e-6_dp, &
         za5 = 2.034080948e-8_dp, &
         za6 = 6.136820929e-11_dp 

    REAL(dp) :: zt
    
    REAL(dp) :: psat

    zt = ptemp - 273.15_dp

    psat = za0 + za1*zt + za2*zt**2 + za3*zt**3 +   &
           za4*zt**4 + za5*zt**5 + za6*zt**6

    ! To Pascals
    psat = psat*100._dp

  END FUNCTION satvaph2o
!
! --------------------------------------------------------------------------------------------------------------
!
  
  !------------------------------------------------
  !
  ! ***************
  ! FUNCTION coagc
  ! ***************
  !
  ! Calculation of coagulation coefficients.
  ! Extended version of the function originally 
  ! found in mo_salsa_init. This is now placed
  ! here to avoid cyclic dependencies between 
  ! modules upon coupling with UCLALES.
  !
  ! J. Tonttila, FMI, 05/2014 
  !
  !-------------------------------------------------
  FUNCTION coagc(diam1,diam2,mass1,mass2,temp,pres,kernel)

    USE mo_kind, ONLY : dp

    USE mo_submctl, ONLY : pi, pi6, boltz, pstand, grav
    USE mo_constants, ONLY : rd, amd

    IMPLICIT NONE

    !-- Input variables ----------
    REAL(dp), INTENT(IN) :: &
         diam1,  &   ! diameters of colliding particles [m]
         diam2,  &   !
         mass1,  &   ! masses -"- [kg]
         mass2,  &
         temp,   &   ! ambient temperature [K]
         pres        ! ambient pressure [fxm]

    INTEGER, INTENT(in) :: kernel ! Select the type of kernel: 1 - aerosol-aerosol coagulation (the original version)
                                  !                            2 - hydrometeor-aerosol or hydrometeor-hydrometeor coagulation 

    !-- Output variables ---------
    REAL(dp) ::  &
         coagc       ! coagulation coefficient of particles [m3/s]

    !-- Local variables ----------  
    REAL(dp) ::  &
         visc,   &   ! viscosity of air [kg/(m s)]
         vkin,   &   ! Kinematic viscosity of air [m2 s-1]
         zrhoa,  &   ! Density of air [kg m-3]
         mfp,    &   ! mean free path of air molecules [m]
         mdiam,  &   ! mean diameter of colliding particles [m]
         fmdist, &   ! distance of flux matching [m]
         zecoll, &   ! Collition efficiency for graviational collection
         zev,    &   ! 
         zea,    &
         zbrown, &   ! Components for coagulation kernel; Brownian
         zbrconv,&   !                                    Convective diffusion enhancement
         zgrav       !                                    Gravitational collection
    

    REAL(dp), DIMENSION (2) :: &
         diam,   &   ! diameters of particles [m]
         mpart,  &   ! masses of particles [kg]
         knud,   &   ! particle knudsen number [1]
         beta,   &   ! Cunningham correction factor [1]
         zrhop,  &   ! Particle density [kg m-3]
         dfpart, &   ! particle diffusion coefficient [m2/s]
         mtvel,  &   ! particle mean thermal velocity [m/s]
         termv,  &   ! Particle terminal velocity
         omega,  &   !
         tva,    &   ! temporary variable [m]
         flux        ! flux in continuum and free molec. regime [m/s]


    REAL(dp), DIMENSION(2) ::  &
         schm,   &   ! Schmidt nubmer 
         reyn        ! Reynolds number
    REAL(dp) :: &
         stok        ! Stokes number
    INTEGER :: lrg,sml 

    !------------------------------------------------------------------------------- 

    !-- 0) Initializing particle and ambient air variables --------------------
    diam = (/ diam1, diam2 /)       ! particle diameters [m]
    mpart = (/ mass1, mass2 /)       ! particle masses [kg]

    visc = (7.44523e-3_dp*temp**1.5_dp)/ &
         (5093._dp*(temp+110.4_dp))                   ! viscosity of air [kg/(m s)]
 
    mfp = (1.656e-10_dp*temp+1.828e-8_dp)*pstand/pres ! mean free path of air [m]


    !-- 2) Slip correction factor for small particles -------------------------

    knud = 2._dp*mfp/diam                                    ! Knudsen number
    beta = 1._dp+knud*(1.142_dp+0.558_dp*exp(-0.999_dp/knud))! Cunningham correction factor
    ! (Allen and Raabe, Aerosol Sci. Tech. 4, 269)

    !-- 3) Particle properties ------------------------------------------------

    dfpart = beta*boltz*temp/(3._dp*pi*visc*diam)  ! diffusion coefficient [m2/s]
    mtvel = sqrt((8._dp*boltz*temp)/(pi*mpart))    ! mean thermal velocity [m/s]
    omega = 8._dp*dfpart/(pi*mtvel)

    mdiam = 0.5_dp*(diam(1)+diam(2))               ! mean diameter [m]

    !-- 4) Calculation of fluxes and flux matching ----------------------------

    flux(1) = 4._dp*pi*mdiam*( dfpart(1)+dfpart(2) )    ! flux in continuum regime [m3/s]
    flux(2) = pi*sqrt((mtvel(1)**2)+(mtvel(2)**2))*(mdiam**2) !  -"- in free molec. regime [m3/s]

    tva(1) = ((mdiam+omega(1))**3 - &              ! temporary variable [m]
         (mdiam**2+omega(1)**2)* &
         sqrt((mdiam**2+omega(1)**2)))/ &
         (3._dp*mdiam*omega(1)) - mdiam

    tva(2) = ((mdiam+omega(2))**3 - &              ! temporary variable [m]
         (mdiam**2+omega(2)**2)* &
         sqrt((mdiam**2+omega(2)**2)))/ &
         (3._dp*mdiam*omega(2)) - mdiam

    fmdist = sqrt(tva(1)**2+tva(2)**2)             ! flux matching distance [m]
    
    SELECT CASE(kernel)
       CASE(1)

          ! Aerosol-Aerosol coagulation - like the original version
          !-- 5) Coagulation coefficient [m3/s] -------------------------------------
          coagc = flux(1) / (mdiam/(mdiam+fmdist) + flux(1)/flux(2)) 

       CASE(2)
          
          ! Which particle is larger?
          sml = 1; lrg = 2
          IF (diam(1) >= diam(2)) THEN
             lrg = 1; sml = 2
          END IF 

          zrhoa = pres/(rd*temp)   ! Density of air
          zrhop = mpart/(pi6*diam**3)             ! Density of particles
          vkin = visc/zrhoa   ! Kinematic viscosity of air [m2 s-1]
          termv = ( (diam**2) * (zrhop - zrhoa) * grav * beta )/( 18._dp*visc  ) ![m s-1]          

          ! Reynolds number
          reyn = diam*termv/vkin
          ! Schmidt number for the smaller particle
          schm = vkin/dfpart
          ! Stokes number
          stok = 2._dp*termv(sml)*ABS(termv(1) - termv(2))/( diam(lrg)*grav )
          
          !Brownian component
          zbrown = flux(1) / (mdiam/(mdiam+fmdist) + flux(1)/flux(2)) 

          ! Convective enhancement
          IF (reyn(lrg) <= 1._dp) THEN
             zbrconv = 0.45_dp*zbrown*( reyn(lrg)**(1._dp/3._dp) )*( schm(sml)**(1._dp/3._dp) )
          ELSE IF (reyn(lrg) > 1._dp) THEN
             zbrconv = 0.45_dp*zbrown*SQRT(reyn(lrg))*( schm(sml)**(1._dp/3._dp) )
          END IF
          
          ! gravitational collection
          zea = stok**2/( stok + 0.5_dp )**2
          IF (stok > 1.214_dp) THEN
             zev = 0.75_dp*LOG(2._dp*stok)/(stok - 1.214)
             zev = (1._dp + zev)**(-2._dp)
          ELSE IF (stok <= 1.214) THEN
             zev = 0._dp
          END IF
          
          zecoll = (60._dp*zev + zea*reyn(lrg))/(60._dp + reyn(lrg))
          zgrav = zecoll * pi * mdiam**2
          zgrav = zgrav * ABS(termv(1)-termv(2))

          ! Total coagulation kernel
          coagc = zbrown  + zbrconv + zgrav

    END SELECT

  END FUNCTION coagc



END MODULE mo_salsa_dynamics
