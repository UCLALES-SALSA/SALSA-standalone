C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C                        Program AEROTR                                +
C                   For WWW or file-input use                          +
C                                                                      +
C                                                                      +
C            Development:  22/2/98 - tested vs original aeromix.       +
C                          3/3/98  - improved error messages           +
C                          6/6/98  - inclusion of HBr. Results         +
C                                    validated against the original    +
C                                    program OK.                       +
C                          15/7/98 - minor change to main program.     +
C                                                                      +
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C ======================================================================
C ======================================================================
C
      SUBROUTINE THERMO(T,MOLAL,pHNO3,pNH3)
      IMPLICIT REAL*8(A-H,O-Z)
C
      PARAMETER(NCmax=3, NAmax=5, NNmax=1)
      PARAMETER(NCOEF=80)
C
      COMMON/PARAM/COEFF(NCOEF)
      COMMON/SPEC1ES/NC,NA,NN,nCAT(NCmax),nAN(NAmax),nNEUT(NNmax)
      COMMON/COEFF1/RHO,A1ca(NCmax,NAmax),A2ca(NCmax,NAmax),
     >       B1ca(NCmax,NAmax),B2ca(NCmax,NAmax),
     >       Wnca(NNmax,NCmax,NAmax),Unca(NNmax,NCmax,NAmax),
     >       Vnca(NNmax,NCmax,NAmax),Wcca(NCmax,NCmax,NAmax),
     >       Waac(NAmax,NAmax,NCmax),Qncca(NNmax,NCmax,NCmax,NAmax),
     >       Qnaac(NNmax,NAmax,NAmax,NCmax),Ucca(NCmax,NCmax,NAmax),
     >       Uaac(NAmax,NAmax,NCmax),Wnn(NNmax,NNmax),Unn(NNmax,NNmax)
      COMMON/COEFF2/Ynnca(NNmax,NNmax,NCmax,NAmax),
     >              Xaaac(NAmax,NAmax,NAmax,NCmax),
     >              Xccca(NCmax,NCmax,NCmax,NAmax),
     >              Zccaa(NCmax,NCmax,NAmax,NAmax),
     >              Cnnn(NNmax,NNmax,NNmax)
      COMMON/CHARGE/CATchrg(NCmax),ANchrg(NAmax),
     >              Vca(NCmax,NAmax),Vac(NAmax,NCmax)
      INTEGER SOLID_R(1:30,4),SOLID_S(1:30,4),SOLID_N(1:30)
      COMMON/SLDS/SOLID_R,SOLID_S,SOLID_N
C
C
      LOGICAL Test_Charge_Bal,WWW
      INTEGER iSLD(25),refSalt(10)
      REAL*8  NH3g,MOLAL(-NAmax:NCmax),MFRAC(-NAmax:NCmax),
     >        ACT(-NAmax:NCmax), SRATIO(25), mSalt(10), L1, Ls,
     >        J1, Js
C
C ===================================================================
C ================ EXECUTABLE STATEMENTS FOLLOW =====================
C
C    -------------------------- 
C   | Set the mode of use here |
C    -------------------------- 
      WWW = .FALSE.             
C    -------------------------- 
C   |--------------------------|
C    -------------------------- 
C
C      ww=100.-wt1-wt2-wt3-wt4
C      xm1=1000./98.*wt1/ww
C      xm2=1000./132.*wt2/ww
C      xm3=1000./63.*wt3/ww
C      xm4=1000./36.5*wt4/ww
C      x_leto = pleto/ww
C
      infault0=0
C       DO 2 I=-NAmax,NCmax
C         MOLAL(I)=0.d0
C 2    CONTINUE 
C     
C      MOLAL(-3)=xm3
C      MOLAL(-2)=xm1+xm2
C      MOLAL(-4)=0.d0
C      MOLAL(2)=2*xm2
C      MOLAL(1)=MOLAL(-3)+2*MOLAL(-2)+MOLAL(-4)-MOLAL(2)
c      IF(MOLAL(1).LT.0.d0) THEN
c         write(*,*) molal(1),molal(2),2.*molal(-2),molal(-3),xleto
c         write(*,*) wt1,xleto
C         MOLAL(1)=0.d0
C         MOLAL(2)=MOLAL(-3)+2*MOLAL(-2)+MOLAL(-4)
c      ENDIF
C
C               Ref.  Solid               Ref.  Solid              
C                 1   ice                  10   (NH4)2SO4          
C                (2   H2SO4) [a]           11   (NH4)3H(SO4)2      
C                 3   H2SO4.H2O            12   NH4HSO4            
C                 4   H2SO4.2H2O           13   NH4NO3             
C                 5   H2SO4.3H2O           14   (NH4)2SO4.2NH4NO3  
C                 6   H2SO4.4H2O           15   (NH4)2SO4.3NH4NO3  
C                 7   H2SO4.6.5H2O         16   NH4NO3.NH4HSO4     
C                 8   HNO3.H2O                                     
C                 9   HNO3.3H2O            28   HCl.3H2O           
C                                          29   HBr  [b]           
C      nSLD = 1
C      iSLD(1) = 11
C
C   ..set machine (double) precision to some reasonable value.
      EPS=2.0D-16
C      EPS=1.D-5
C      EPS=1.D-14
C
C   ..set temperature bounds
      T_Low_Bound=180.D0
      T_High_Bound=330.D0
C     
C      -----------------------------------
C     | (2) Check for negative molalities |
C      -----------------------------------
      DO 3 J=-NAmax,NCmax
        IF(MOLAL(J).LT.0.D0) THEN
          WRITE(*,'(1X,''Species input molality <0: STOP.'')')
          WRITE(*,*) MOLAL
          STOP
        ENDIF
3     CONTINUE
C
C      ---------------------------
C     | (3) Check for errors in T |
C      ---------------------------
      IF(T.LT.T_Low_Bound .OR. T.GT.T_High_Bound) THEN
        WRITE(*,105)
        STOP
      ENDIF
C
C      --------------------------------------
C     | (4) Charge balance test and revision |
C      --------------------------------------
      CALL Charge_Balance(Molal,Test_Charge_Bal)
      IF(.NOT.Test_Charge_Bal) THEN
         WRITE(*,*) MOLAL
        WRITE(*,125)
        STOP
      ENDIF
C
C      -------------------------------------
C     | (5a) test for Cl- and NH4+ together |
C      -------------------------------------
      IF(Molal(-4).GT.0.D0 .AND. Molal(2).GT.0.D0) THEN
        WRITE(*,127)
        STOP
      ENDIF
C
C      -------------------------------------
C     | (5a) test for Br- and NH4+ together |
C      -------------------------------------
      IF(Molal(-5).GT.0.D0 .AND. Molal(2).GT.0.D0) THEN
        WRITE(*,128)
        STOP
      ENDIF
C
C
C      ---------------------------------------------------------
C     | (6a) First calculate the speciation, returning the mole |
C     |      fractions and activities, and the water activity:  |
C      ---------------------------------------------------------
      CALL SPEC(T,MOLAL,EPS,MFRAC,ACT,AW) 
C
C      ---------------------------------------
C     | (6b) Calculate the equilibrium vapour |
C     |      pressure (atm) of H2O:           |
C      ---------------------------------------
      pH2O=H2Og(T,aw)*101325.
C
C      ---------------------------------------
C     | (6c) Calculate the equilibrium vapour |
C     |      pressure (atm) HNO3:             |
C      ---------------------------------------
      pHNO3=HNO3g(T,MFRAC(1),MFRAC(-3),ACT(1),ACT(-3))*101325.
      actHNO3=ACT(-3)
C
C      ----------------------------------------
C     | (6d)  Calculate the equilibrium vapour |
C     |       pressure (atm) HCl:              |
C      ----------------------------------------
C      pHCl=HClg(T,MFRAC(1),MFRAC(-4),ACT(1),ACT(-4))*1013.25
C
C      ------------------------------------------------
C     | (6e) Calculate the equilibrium vapour pressure |
C     |      (atm) of NH3.                             |
C      ------------------------------------------------
      pNH3=NH3g(T,MFRAC(1),MFRAC(2),ACT(1),ACT(2))*101325.
      actNH3=ACT(2)
C
C      ------------------------------------------------
C     | (6f) Calculate the equilibrium vapour pressure |
C     |      (atm) of H2SO4.                           |
C      ------------------------------------------------
C      pH2SO4=H2SO4g(T,MFRAC(1),MFRAC(-2),ACT(1),ACT(-2))
C
C      ----------------------------------------
C     | (6g)  Calculate the equilibrium vapour |
C     |       pressure (atm) HBr:              |
C      ----------------------------------------
C      pHBr=HBrg(T,MFRAC(1),MFRAC(-5),ACT(1),ACT(-5))
C
C      -------------------------------------------
C     | (6h) Calculate degrees of saturation with |
C     |      respect to specified solid phases:   |
C      -------------------------------------------
C      CALL CALCSAT(nSLD,iSLD,T,MFRAC,ACT,SRATIO)
C
C      ---------------------------------------------------
C     | (6i) Calculate partial molar enthalpies and heat  |
C     |      capacities of the water ('1') and a selected |
C     |      solute ion ('s').                            |
C      ---------------------------------------------------
C      IF(.NOT.WWW) THEN
C        CALL fL1J1(T,MOLAL,EPS,I1,L1,J1)
C        CALL fLsJs(T,MOLAL,EPS,Is,Ls,Js)
C      ENDIF
C
C     -----------------------------------
C     Read and end-of-file error handling
C     -----------------------------------
C   ..header error/eof
C200   IF(InFault0.GT.0) THEN
C        WRITE(*,2000) nProblem
C      ELSEIF(Infault0.LT.0) THEN
C        WRITE(*,2001) nProblem
C      ENDIF
C      STOP
CC
CC   ..in_choose error/eof
C205   IF(InFault5.GT.0) THEN
C        WRITE(*,2050) nProblem
C      ELSEIF(Infault5.LT.0) THEN
C        WRITE(*,2051) nProblem
C      ENDIF
C      STOP
CC
CC   ..input salts error/eof
C210   IF(InFault10.GT.0) THEN
C        WRITE(*,2100) nProblem
C      ELSEIF(Infault10.LT.0) THEN
C        WRITE(*,2101) nProblem
C      ENDIF
C      STOP
CC
CC   ..solid saturation error/eof
C215   IF(InFault15.GT.0) THEN
C        WRITE(*,2150) nProblem
C      ELSEIF(Infault15.LT.0) THEN
C        WRITE(*,2151) nProblem
C      ENDIF
C      STOP
CC
CC   ..input temp & electrolyte molalities error/eof
C220   IF(InFault20.GT.0) THEN
C        WRITE(*,2200) nProblem-1
C      ELSEIF(Infault20.LT.0) THEN
C        WRITE(*,2201) nProblem-1
C      ENDIF
C      STOP
CC
CC   ..input temp & ion molalities error/eof
C225   IF(InFault25.GT.0) THEN
C        WRITE(*,2250) nProblem-1
C      ELSEIF(Infault25.LT.0) THEN
C        WRITE(*,2251) nProblem-1
C      ENDIF
C      STOP
CC
C
90    FORMAT(/12X,
     >' ====================================================',/12X,
     >'||         Aerosol Inorganics Models I & II         ||',/12X,
     >'||                                                  ||',/12X,
     >'||  Systems: H-NH4-SO4-NO3-H2O,                     ||',/12X,
     >'||           H-SO4-NO3-Cl-Br-H2O                    ||',/12X,
     >'||  Temperatures: 180 to 330 K.                     ||',/12X,
     >'||                                                  ||',/12X,
     >'||  AEROTR, version 1.0d, compiled 12/6/98          ||',/12X,
     >'||  Contact: S. L. Clegg (s.clegg@uea.ac.uk)        ||',/12X,
     >' ====================================================',/1X)
C
105   FORMAT(/1X,'Error - temperature out of range. Stop.')
125   FORMAT(/1X,'Error - charges of cations and anions do not '/
     >           ' balance. Stop.'/
     >           ' Suggestion: check the input composition.')
127   FORMAT(/1X,'Error: Cl- not permitted in this model '/
     >  ' when NH4+ is also present. Stop.'/
     >  ' Suggestion: revise the input composition.')
128   FORMAT(/1X,'Error: Br- not permitted in this model '/
     >  ' when NH4+ is also present. Stop.'/
     >  ' Suggestion: revise the input composition.')
130   FORMAT(1X,'Error: input electrolyte reference number (',I3,')',
     >          /1X,'is invalid. Stop in main program.')

135   FORMAT(1X,'Error: solid phase reference number (',I3,')',
     >          /1X,'is out of range. Stop in main program.')
C
100   FORMAT(1X,F6.2,2X,F7.5,2X,25(E12.5,1X))
110   FORMAT(1X,F6.2,2X,6(E11.4,1X),1X,25(E12.5,1X))
101   FORMAT(1X,'    T      AW          xBr         xCl         xNO3',
     >'         xSO4         xHSO4        xH2O         xH',
     >'           xNH4         fBr          fCl          fNO3',
     >'         fSO4        fHSO4         fH2O         fH',
     >'           fNH4')
102   FORMAT(1X,'    T        pH2O       pHNO3        pNH3 ',
     >'       pHCl      pHBr      pH2SO4       mBr          mCl',
     >'          mNO3          mSO4         mHSO4        mH',
     >'           mNH4',8X,25(A6,I2,5X))      
103   FORMAT(1X,'    T         L1          J1          Ls  ',
     >'       Js')
C
2000  FORMAT(/1X,
     >'Error reading header in aerotr.dat. Stop before problem ',I7,'.')
2001  FORMAT(/1X,
     >'Error - end of file reached while reading header in ',/1X,
     >'aerotr.dat. Stop before problem ',I7,'.')
2050  FORMAT(/1X,
     >'Error - read error encountered while reading in_choose on',/1X,
     >'line 2 of aerotr.dat. Stop before problem ',I7,'.')
2051  FORMAT(/1X,
     >'Error - end of file reached while reading in_choose on',/1X,
     >'line 2 of aerotr.dat. Stop before problem ',I7,'.')
2100  FORMAT(/1X,
     >'Error - read error encountered while reading the reference',/1X,
     >'numbers of input electrolytes (for in_choose=1) in aerotr.dat.',
     >/1X,'Stop before problem ',I7,'.')
2101  FORMAT(/1X,
     >'Error - end of file reached while reading the reference',/1X,
     >'numbers of input electrolytes (for in_choose=1) in aerotr.dat.',
     >/1X,'Stop before problem ',I7,'.')
2150  FORMAT(/1X,
     >'Error - read error encountered while reading the reference',/1X,
     >'numbers of solids to be checked for saturation in aerotr.dat.',
     >/1X,'Stop before problem ',I7,'.')
2151  FORMAT(/1X,
     >'Error - end of file encountered while reading the reference',/1X,
     >'numbers of solids to be checked for saturation in aerotr.dat.',
     >/1X,'Stop before problem ',I7,'.')
2200  FORMAT(/1X,
     >'Read error encountered while reading temperature',/1X,
     >'and input salt molalities (for in_choose=1) in aerotr.dat.',
     >/1X,'Stop after problem ',I7,'.')
2201  FORMAT(/1X,
     >'End of file reached.  No. of problems solved = ',I7,'.')
2250  FORMAT(/1X,
     >'Read error encountered while reading temperature',/1X,
     >'and input ion molalities (for in_choose=2) in aerotr.dat.',
     >/1X,'Stop after problem ',I7,'.')
2251  FORMAT(/1X,
     >'End of file reached. No. of problems solved = ',I7,'.')
C
      END
C
C >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
C              beginning of aerotr code for web input/output
C <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
C >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
C <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
C
      SUBROUTINE Read_Web_Input(T,Molal,nSLD,iSolid)
      IMPLICIT REAL*8 (A-H,O-Z)
C
C  ========================================================
C  |  This routine takes input data from the AIM web pages |
C  ========================================================
C
      PARAMETER(NCmax=3, NAmax=5)
C
C     ---- local ---
      PARAMETER(nNamemax=38)
      CHARACTER Name(nNamemax)*20, LineIn*500
      INTEGER Name_length(nNamemax),endofnumber,startofnumber,type_size
C
C     ---- arguments ----
      INTEGER      iSolid(25)
      REAL*8       Molal(-NAmax:NCmax)
C
C
C      ------------------------------------------------------------- 
C     | Initialise Name for the variables present in the character  |
C     | string POSTED by the web page. Also initialise an integer   |
C     | variable giving the length of each name.                    |
C     | Make sure each Name is unambiguous in that it does not      |
C     | constitute a part of another one sent by the web page.      |
C      ------------------------------------------------------------- 
      DATA (Name(I),Name_length(I),I=1,nNamemax)/   
C    --- system and special conditions ---                      
     >   'temp=',5,
C
C    --- moles of ions (listed in order of reference no.) ---
     >   'hydrogen=',9,    'ammonium=',9,    'sodium=',7,
     >   'bisulphate=',11, 'sulphatex=',10,  'nitrate=',8,
     >   'chloride=',9, 'bromide=',8,
C
C    --- solids options (for saturation ratios) ---
     >   'ice=',4,                          'h2so4=',6, 
     >   'h2so4_h2o=',10,                   'h2so4_2h2o=',11,
     >   'h2so4_3h2o=',11,                  'h2so4_4h2o=',11,
     >   'h2so4_65h2o=',12,                 'hno3_h2o=',9,
     >   'hno3_3h2o=',10,                   'nh42so4x=',9,
     >   'nh43hso42=',10,                   'nh4hso4x=',9,
     >   'nh4no3=',7,                       '2nh4no3_nh42so4=',16,
     >   '3nh4no3_nh42so4=',16,             'nh4no3_nh4hso4=',15,
     >   'nh4cl=',6,                        'na2so4=',7,
     >   'na2so4_10h2o=',13,                'na3hso42=',9,
     >   'nahso4_h2o=',11,                  'nahso4=',7,
     >   'nah3so42_h2o=',13,                'na2so4_nh42so4_4h2o=',20,
     >   'nano3=',6,                        'nano3_na2so4_h2o=',17,
     >   'nacl=',5,                         'hcl_3h2o=',9,
C
C    --- other options ---
     >   'type_size=',10/        
C
C
C
C      ==========================================
C     | (1) Read the webpage data as a character |
C     | string, and determine the length of the  |
C     | string. Append a final '&' to mark the   |
C     | end of the last numerical value.         |
C      ==========================================
      READ(*,'(A)',ERR=999) LineIn
      linelength=INDEX(LineIn,' ')-1
C     ..zero length, stop on error.
      IF(linelength.EQ.0) Goto 999   
C
      linelength=linelength+1
      LineIn(linelength:linelength)='&'
C
C   ..print the input string to the standard output.
C     write(*,'(A)') linein
C                           
C
C   ..initialise number of solids for saturation ratio checking
      nSLD=1
C
C      ===================================
C     | (2) Read all the information in   |
C     | the character string 'LineIn',    |
C     | and assign to FORTRAN variables.  |
C      =================================== 
      DO 1 I=1,nNamemax
C     ..index of first character of Name in LineIn
        istart=INDEX(LineIn,Name(I)(1:Name_length(I)))
        IF(istart.GT.0) THEN
C         -------------------------------------------------------------
C        | Name(i) occurs in the input string, so find the indices     |
C        | of the initial and final characters of the numerical value  |
C        | associated with it. Note that we depend upon finding an '&' |
C        | after the final entry in the string (added in (1) above).   |
C         ------------------------------------------------------------- 
          startofnumber=istart+Name_length(I)  
          numberlength=INDEX(LineIn(startofnumber:linelength),'&')-1
          endofnumber=startofnumber + numberlength - 1
C
          IF(endofnumber.LT.startofnumber) THEN
C         ..this means a null value has been transmitted
C           (ie. Name=&), so we set var to zero:
            var=0.D0
          ELSE
C         ..we can read the number into a real variable
            READ(LineIn(startofnumber:endofnumber),*) var
          ENDIF
C
C       ..write the results of the 'reads' to standard output
C         write(*,*) I,Name(I)(1:Name_Length(I)),startofnumber,
C    >               endofnumber,var
C       
C
        ELSE
C         ------------------------------
C        | Name(i) does not occur, so   |
C        | the variable is set to zero. |
C         ------------------------------ 
          var=0.D0
        ENDIF
C
C       ------------------------------------------
C      | Assign values to T, the water variable,  |
C      | Moles_Species, Options_Species and pGas. |
C      | Note that any modification of the order  |
C      | in Names will require changes in the IF  |
C      | statements below.                        |
C       ------------------------------------------
        IF(I.EQ.1) THEN
C       ..temperature
          T=DBLE(var) 
        ELSEIF(I.GE.2 .AND. I.LE.4) THEN
C       ..Moles of cations
          iCat=I-1
          Molal(iCat)=DBLE(var)
        ELSEIF(I.GE.5 .AND. I.LE.9) THEN
C       ..Moles of anions
          iAn=-(I-4)
          Molal(iAn)=DBLE(var)
        ELSEIF(I.EQ.38) THEN
C       .. type_size (0=normal, 1=small)          
          type_size=INT(var)

        ENDIF
C
C     ..assign solids whose saturation ratios are to be checked. The
C       reference numbers are the same as in the gibbs minimisation
C       codes.
        IF(I.GE.10 .AND. I.LE.37 .AND. INT(var).EQ.1) nSLD=nSLD+1
C       ..H2O         
        IF(I.EQ.10 .AND. INT(var).EQ.1)  iSolid(nSLD)=1    
C       ..H2SO4.H2O   
        IF(I.EQ.12 .AND. INT(var).EQ.1) iSolid(nSLD)=3    
C       ..H2SO4.2H2O  
        IF(I.EQ.13 .AND. INT(var).EQ.1) iSolid(nSLD)=4    
C       ..H2SO4.3H2O  
        IF(I.EQ.14 .AND. INT(var).EQ.1) iSolid(nSLD)=5    
C       ..H2SO4.4H2O  
        IF(I.EQ.15 .AND. INT(var).EQ.1) iSolid(nSLD)=6    
C       ..H2SO4.6.5H2O
        IF(I.EQ.16 .AND. INT(var).EQ.1) iSolid(nSLD)=7    
C       ..HNO3.H2O   
        IF(I.EQ.17 .AND. INT(var).EQ.1) iSolid(nSLD)=8    
C       ..HNO3.3H2O   
        IF(I.EQ.18 .AND. INT(var).EQ.1) iSolid(nSLD)=9    
C       ..(NH4)2SO4           
        IF(I.EQ.19 .AND. INT(var).EQ.1) iSolid(nSLD)=10   
C       ..(NH4)3H(SO4)2        
        IF(I.EQ.20 .AND. INT(var).EQ.1) iSolid(nSLD)=11   
C       ..NH4HSO4          
        IF(I.EQ.21 .AND. INT(var).EQ.1) iSolid(nSLD)=12   
C       ..NH4NO3           
        IF(I.EQ.22 .AND. INT(var).EQ.1) iSolid(nSLD)=13   
C       ..2NH4NO3.(NH4)2SO4
        IF(I.EQ.23 .AND. INT(var).EQ.1) iSolid(nSLD)=14   
C       ..3NH4NO3.(NH4)2SO4
        IF(I.EQ.24 .AND. INT(var).EQ.1) iSolid(nSLD)=15   
C       ..NH4NO3.NH4HSO4   
        IF(I.EQ.25 .AND. INT(var).EQ.1) iSolid(nSLD)=16   
C       ..NH4Cl
        IF(I.EQ.26 .AND. INT(var).EQ.1) iSolid(nSLD)=17   
C       ..Na2SO4          
        IF(I.EQ.27 .AND. INT(var).EQ.1) iSolid(nSLD)=18   
C       ..Na2SO4.10H2O     
        IF(I.EQ.28 .AND. INT(var).EQ.1) iSolid(nSLD)=19   
C       ..Na3H(SO4)2    
        IF(I.EQ.29 .AND. INT(var).EQ.1) iSolid(nSLD)=20   
C       ..NaHSO4.H2O    
        IF(I.EQ.30 .AND. INT(var).EQ.1) iSolid(nSLD)=21   
C       ..NaHSO4          
        IF(I.EQ.31 .AND. INT(var).EQ.1) iSolid(nSLD)=22   
C       ..NaH3(SO4)2.H2O  
        IF(I.EQ.32 .AND. INT(var).EQ.1) iSolid(nSLD)=23   
C       ..Na2SO4.(NH4)2SO4.4H2O
        IF(I.EQ.33 .AND. INT(var).EQ.1) iSolid(nSLD)=24   
C       ..NaNO3             
        IF(I.EQ.34 .AND. INT(var).EQ.1) iSolid(nSLD)=25   
C       ..NaNO3.Na2SO4.H2O  
        IF(I.EQ.35 .AND. INT(var).EQ.1) iSolid(nSLD)=26   
C       ..NaCl            
        IF(I.EQ.36 .AND. INT(var).EQ.1) iSolid(nSLD)=27   
C       ..HCl.3H2O
        IF(I.EQ.37 .AND. INT(var).EQ.1) iSolid(nSLD)=28   
C
1     CONTINUE
C
C
C      ---------------------------------
C     | Reduce font size, if requested. |
C      ---------------------------------
      IF(type_size.EQ.1) THEN
        WRITE(*,'(A)') '<FONT SIZE="-1">' 
      ENDIF
C
C     ******
      RETURN
C     ******
C
999   WRITE(*,100)
      STOP
C
100   FORMAT(/1X,'Error reading data posted from web page, or ',/1X,
     >           'input line has zero length. Stop Read_Web_Input.')
      END
C
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
      SUBROUTINE Results(T,MOLAL,ACT,AW,pH2O,pHNO3,pHCl,pHBr,pNH3,
     >                   pH2SO4,nSLD,iSolid,SRATIO)
      IMPLICIT REAL*8 (A-H,O-Z)
C
      PARAMETER(NCmax=3, NAmax=5, NNmax=1, NSmax=29, NGmax=6,
     >          NSpmax=NSmax)
C
C  -- arguments --
      INTEGER Options_Species(5,NSpmax),iSolid(25)
      REAL*8 MOLAL(-NAmax:NCmax),ACT(-NAmax:NCmax)
      REAL*8 pGas(NGmax),Moles_Species(5,NSpmax),
     >       Actcoeffs_Species(5,NSpmax),SRATIO(25)
C
C  -- local --
      REAL*8  SatRatio_Solids(NSmax)
C
C
C     -----------------
C     Initialise arrays
C     -----------------
      DO 10 I=1,5
        DO 11 J=1,NSpmax
          Moles_Species(I,J)=0.D0
          Options_Species(I,J)=0
          ActCoeffs_Species(I,J)=0.D0
          IF(I.EQ.1) SatRatio_Solids(J)=0.D0
11      CONTINUE
10    CONTINUE
C
C
      RH = AW
C
C     -------------------------
C     Convert moles of species
C     and activity coefficients 
C     -------------------------
      IF(NNmax.GT.1) THEN
C     ..conversion code below assumes only
C       neutral species is water
        WRITE(*,100) 
        STOP
      ENDIF
C
C   ..moles of water
      Moles_Species(3,1)=55.508681D0
      Actcoeffs_Species(3,1)=ACT(0)
      DO 1 I=1,2
        IF(I.EQ.1) THEN
          jmax=NCmax
          idum=1
        ELSEIF(I.EQ.2) THEN
          jmax=NAmax        
          idum=-1
        ENDIF
        DO 2 J=1,jmax
          Moles_Species(I,J)=MOLAL(idum*J)
          Actcoeffs_Species(I,J)=ACT(idum*J)
2       CONTINUE
C
1     CONTINUE
C
C     --------------------------
C      Convert partial pressures
C     --------------------------
      pGas(1)=pH2O
      pGas(2)=pHNO3
      pGas(3)=pHCl
      pGas(4)=pNH3
      pGas(5)=pH2SO4
      pGas(6)=pHBr
C
C
C     -------------------------------
C     Assign Option=2 to those solids
C     whose saturation ratios are
C     being checked.
C     -------------------------------
      DO 3 I=1,nSLD      
C     ..set option to be a 'variable'
        Options_Species(4,iSolid(I))=2
C     ..assign saturation ratio to new variable
        SatRatio_Solids(iSolid(I))=SRATIO(I)
3     CONTINUE
C
C
      CALL Print_to_Web(T,RH,Options_Species,pGas,
     >                   Moles_Species,Actcoeffs_Species,
     >                   SatRatio_Solids)
C
C
C     ******
      RETURN
C     ******
C
100   FORMAT(/1X,'Error - only a single neutral species (water) is ',
     >       /1X,'expected. Stop in routine CONVERT.')
C
      END
C
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
      SUBROUTINE Print_to_Web(T,RH,Options_Species,pGas,
     >                        Moles_Species,Actcoeffs_Species,
     >                        SatRatio_Solids)
      IMPLICIT REAL*8(A-H,O-Z)
C
C
      PARAMETER(NCmax=3,  NAmax=5, NNmax=1)
      PARAMETER(NSmax=29, NGmax=6, NSpmax=NSmax)
C
C
C     -- local --
      CHARACTER*22 Names_Species_R(5,NSpmax)
      LOGICAL Sat_Check,ldum
      REAL*8  Mole_Fraction,Molality,Molar_Masses_R(5,NSpmax)

C     -- arguments --
      INTEGER Options_Species(5,NSpmax)
      REAL*8  Moles_Species(5,NSpmax),Actcoeffs_Species(5,NSpmax),
     >        SatRatio_Solids(NSmax),pGas(NGmax)
C
C   ---------------------------------------------------------------
C  | Initialize a reference array. Molar_Masses_R holds the gram   |
C  | molecular mass of each cation, anion, neutral, solid and gas. |              |
C   ---------------------------------------------------------------
      DATA (Molar_Masses_R(1,J),J=1,NCmax)/
     >     1.0079D0,18.0383D0,22.990D0/
      DATA (Molar_Masses_R(2,J),J=1,NAmax)/
     >     97.0655D0,96.058D0,62.005D0,35.453D0,79.904D0/
      DATA (Molar_Masses_R(3,J),J=1,NNmax)/
     >     18.0152D0/
      DATA (Molar_Masses_R(4,J),J=1,NSmax)/                 
     >     18.0152D0,98.073D0,116.089D0,134.104D0,152.119D0,170.134D0,
     >     215.172D0,81.028D0,117.058D0,132.142D0,247.238D0,115.104D0,
     >     80.043D0,292.228D0,372.271D0,195.147D0,53.491D0,142.037D0,
     >     322.189D0,262.092D0,138.070D0,120.055D0,236.144D0,346.240D0,
     >     84.995D0,245.047D0,58.443D0,90.5065D0,80.9119D0/
      DATA (Molar_Masses_R(5,J),J=1,NGmax)/
     >     18.0152D0,63.013D0,36.461D0,17.030D0,98.073D0,80.9119D0/
C
C   ------------------------------------------------------------------
C  | Initialize a reference array. Names_Species_R holds the names of |
C  | each cation, anion, neutral, solid and gas species.              |
C   ------------------------------------------------------------------
      DATA (Names_Species_R(1,I),I=1,NCmax)/'H(aq)','NH4(aq)','Na(aq)'/
      DATA (Names_Species_R(2,I),I=1,NAmax)/'HSO4(aq)','SO4(aq)',
     >'NO3(aq)','Cl(aq)','Br(aq)'/
      DATA (Names_Species_R(3,I),I=1,NNmax)/'H2O(l)'/
      DATA (Names_Species_R(4,I),I=1,NSmax)/'H2O(s)','H2SO4',
     >'H2SO4.H2O','H2SO4.2H2O','H2SO4.3H2O','H2SO4.4H2O','H2SO4.6.5H2O',
     >'HNO3.H2O','HNO3.3H2O','(NH4)2SO4','(NH4)3H(SO4)2','NH4HSO4',
     >'NH4NO3','2NH4NO3.(NH4)2SO4','3NH4NO3.(NH4)2SO4',
     >'NH4NO3.NH4HSO4','NH4Cl','Na2SO4','Na2SO4.10H2O',
     >'Na3H(SO4)2','NaHSO4.H2O','NaHSO4','NaH3(SO4)2.H2O',
     >'Na2SO4.(NH4)2SO4.4H2O','NaNO3','NaNO3.Na2SO4.H2O','NaCl',
     >'HCl.3H2O','HBr'/
      DATA (Names_Species_R(5,I),I=1,NGmax)/'H2O(g)','HNO3(g)','HCl(g)',
     >'NH3(g)','H2SO4(g)','HBr(g)'/
C
C
      Sat_Check = .FALSE.
C      ---------------------------------------------------
C     | If there is a liquid phase, then look through the |
C     | options for solids whose saturation ratios should |
C     | be checked.                                       |
C      ---------------------------------------------------
      DO 10 I=1,NSmax
        IF(Options_Species(4,I).EQ.2) Sat_Check = .TRUE.
10    CONTINUE
C
C
C
C
C     --------------------------------------
C    |              .RS1                    |
C    | Write problem number and iFail value |
C     --------------------------------------
      iProblem=1
      WRITE(*,195) iProblem
C
C
C     ----------------------------------------
C    |              .RS1                      |
C    | into RH and write out with T,P,RH,pH2O |
C     ----------------------------------------
      WRITE(*,200) T,RH,pGas(1)       
C
C
C
C      ---------------------------------------------
C     |                   .RS1                      |
C     | If we have a liquid phase, then write out   |
C     | the numbers of moles of cations anions and  |
C     | neutral species present, together with mole |
C     | fractions and activity coefficients         |
C      ---------------------------------------------
      WRITE(*,210) 
      Sum=Sum_Liquid_Moles(Moles_Species)
      ivar=0
      DO 20 I=1,NCmax+NAmax+NNmax
        ivar=ivar+1
        IF(ivar .LE. NCmax) THEN
C       ..print cation values
          iType = 1
          iSpecies = ivar
        ELSEIF(ivar.GT.NCmax .AND. ivar.LE.(NCmax+NAmax)) THEN
C       ..print anion values
          iType = 2
        iSpecies = ivar-NCmax
        ELSEIF(ivar.GT.(NCmax+NAmax) .AND. 
     >         ivar.LE.(NCmax+NAmax+NNmax)) THEN
C       ..print neutral species values
          iType = 3
          iSpecies = ivar-(NCmax+NAmax)
        ENDIF
C
        IF(Sum.GT.0.D0) THEN
          Mole_Fraction=Moles_Species(iType,iSpecies)/sum
        ELSE
          Mole_Fraction=0.D0
        ENDIF
C
        Molality=Moles_Species(iType,iSpecies)*55.508681D0/
     >             Moles_Species(3,1)
        IF(Actcoeffs_Species(iType,iSpecies).GT.0.D0) THEN
           WRITE(*,220) Names_Species_R(iType,iSpecies),
     >     Moles_Species(iType,iSpecies),
     >     Moles_Species(iType,iSpecies)*
     >           Molar_Masses_R(iType,iSpecies),Molality,
     >           Mole_Fraction,Actcoeffs_Species(iType,iSpecies)
C
        ELSE
           WRITE(*,220) Names_Species_R(iType,iSpecies),
     >     Moles_Species(iType,iSpecies),
     >     Moles_Species(iType,iSpecies)*
     >           Molar_Masses_R(iType,iSpecies),Molality,
     >           Mole_Fraction
        ENDIF
20    CONTINUE
C
C
C
C      -----------------------------------------
C     |                 .RS1                    |
C     | Write out equilibrium partial pressures | 
C     | calculated from liquid phase quantities.|
C      -----------------------------------------
      ldum=.TRUE.
      DO 45 I=2,NGmax
        IF(pGas(I).GT.0.D0) THEN
C       ..put header above first ratio to be printed
          IF(ldum) WRITE(*,310)
          WRITE(*,320) pGas(I),Names_Species_R(5,I)
          ldum=.FALSE.
        ENDIF
45    CONTINUE
C
C
C
C      -----------------------------------------
C     |                 .RS1                    |
C     | Write out equilibrium partial pressure  | 
C     | products calculated from liquid phase   |
C     | quantities, for systems where there is  |
C     | no H+.                                  |
C      -----------------------------------------
      IF(Moles_Species(1,1).EQ.0.D0) THEN
        CALL Print_Pressures_Liquid(T,Moles_Species,ActCoeffs_Species)
      ENDIF
C
C
C
C      ---------------------------------
C     |            .RS1                 |
C     | Write out the saturation ratio  |
C     | checks (if any required)        |
C      --------------------------------- 
      IF(Sat_Check) THEN
        ldum=.TRUE.
        DO 35 I=1,NSmax
          IF(Options_Species(4,I).EQ.2) THEN
            IF(I.EQ.1 .AND. T.GT.273.15D0) THEN
C           ..skip ice if the temperature is above 0oC
              CONTINUE
            ELSE
C           ..only write ratios greater than 0.01
C           ..put header above first ratio to be printed
              IF(ldum) WRITE(*,280)
              WRITE(*,290) SatRatio_Solids(I),Names_Species_R(4,I)
              ldum=.FALSE.
            ENDIF
          ENDIF
35      CONTINUE
      ENDIF
C
C
C      -------------------------------
C     | Write line of '*' to indicate |
C     | end of this calculation.      |
C      -------------------------------      
      WRITE(*,300)
C
C
C     ******
      RETURN
C     ******
C
195   FORMAT(/1X,'Problem no. ',I4)
200   FORMAT(/1X,'T  = ',F6.2,' K',/1X,
     >          'RH = ',F6.4,' (pH2O = ',E11.4,' atm)')
210   FORMAT(//1X,'*** LIQUID PHASE SPECIES ***',/1X,'Species',5X,
     >            'Moles',9X,'Grams',8X,'Molality',5X,'Mole Frac.',
     >            3X,'Act. Coeff.')
220   FORMAT(1X,A9,2X,E12.5,2X,E11.4,2X,E11.4,2X,E11.4,2X,E11.4)
280   FORMAT(/1X,
     > '- Calculated saturation ratios of solids -',/2X,
     > 'Sat. Ratio',5X,'Species')
290   FORMAT(1X,E11.4,5X,A)
300   FORMAT(//1X,75('*'))
310   FORMAT(/1X,'- Calculated partial pressures (over liquid) -',
     >       /1X,'Pressure (atm)',3X,'Species')
320   FORMAT(1X,E11.4,6X,A)
      END
C
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
      SUBROUTINE Print_Pressures_Liquid(T,Moles_Species,
     >                                  ActCoeffs_Species)
      IMPLICIT REAL*8 (A-H,O-Z)
C
C     *********************************************************************
C    * This routine prints equilibrium partial pressure products, over a   *
C    * liquid phase, for the following reactions:                          *
C    *                                                                     *
C    *          2NH4(aq) + SO4(aq) = 2NH3(g) + H2SO4(g)                    *
C    *          NH4(aq) + NO3(aq) = NH3(g) + HNO3(g)                       *
C    *          NH4(aq) + Cl(aq) = NH3(g) + HCl(g)                         *
C    *                                                                     *
C    * The routine should only be called when (a) there is a liquid phase, *
C    * and (b) there is *no* liquid phase H+ (because if there is, then    *
C    * the partial pressures of each gas can be calculated directly if the *
C    * appropriate ions are present).                                      *
C     *********************************************************************
C
      PARAMETER(NSpmax=28)
C
C
C  -- arguments --
      REAL*8 Moles_Species(5,NSpmax),ActCoeffs_Species(5,NSpmax)
C
C  -- local --
      LOGICAL print_header, NH3_H2SO4, NH3_HNO3, NH3_HCl
C
C
C     --------------------------------
C    | Initialise the variables that  |
C    | determine what will be printed.|
C     --------------------------------
      print_header = .TRUE.
      NH3_H2SO4    = .TRUE.
      NH3_HNO3     = .TRUE.
      NH3_HCl      = .TRUE.
C
C
C     ---------------------------------------
C    | If solids are present (from which the |
C    | pressure products can be calculated), |
C    | then switch off printing here.        |
C     ---------------------------------------
      IF(Moles_Species(4,10).GT.0.D0) NH3_H2SO4 = .FALSE.
      IF(Moles_Species(4,13).GT.0.D0) NH3_HNO3  = .FALSE.      
      IF(Moles_Species(4,17).GT.0.D0) NH3_HCl   = .FALSE.      
C   ..we don't need to print anything, so return.
      IF((.NOT.NH3_H2SO4) .AND. (.NOT.NH3_HNO3) .AND.
     >   (.NOT.NH3_HCl)) THEN
C       ******
        RETURN
C       ******
      ENDIF
C
C
C
C     -------------------------------------------
C    | Get sum of all liquid moles, in order to  |
C    | calculate aqueous mole fractions later.   |
C     -------------------------------------------
      Sum=Sum_Liquid_Moles(Moles_Species)
C
C
C
      IF(Moles_Species(1,2).GT.0.D0 .AND. 
     >   Moles_Species(2,2).GT.0.D0 .AND.
     >   NH3_H2SO4) THEN
C       ---------------------------------------------- 
C      | ammonium and sulphate are present. Calculate |
C      | and print the pressure product.              |
C      | 2NH4(aq) + SO4(aq) = 2NH3(g) + H2SO4(g)      |
C       ---------------------------------------------- 
C
        IF(print_header) THEN
C       ..header, if this is first
C         product to be printed.
          WRITE(*,100)
          print_header=.FALSE.
        ENDIF
C
        eql_constant=SATEQL(26,T)/(xKH_NH3(T)**2*xKH_H2SO4(T))
        xNH4=Moles_Species(1,2)/Sum
        xSO4=Moles_Species(2,2)/Sum
        act_product=(xNH4*ActCoeffs_Species(1,2))**2 *
     >               xSO4*ActCoeffs_Species(2,2)
        p_product=eql_constant * act_product
C
        WRITE(*,110) p_product
C
      ENDIF
C
C
      IF(Moles_Species(1,2).GT.0.D0 .AND. 
     >   Moles_Species(2,3).GT.0.D0 .AND.
     >   NH3_HNO3) THEN
C        ---------------------------------------------
C       | ammonium and nitrate are present. Calculate | 
C       | and print the pressure product.             |
C       | NH4(aq) + NO3(aq) = NH3(g) + HNO3(g)        |
C        ---------------------------------------------
C                                                      
        IF(print_header) THEN
C       ..header, if this is first
C         product to be printed.
          WRITE(*,100)
          print_header=.FALSE.
        ENDIF
C
C
        eql_constant=SATEQL(29,T)/(xKH_NH3(T)*xKH_HNO3(T))
        xNH4=Moles_Species(1,2)/Sum
        xNO3=Moles_Species(2,3)/Sum
        act_product=xNH4*ActCoeffs_Species(1,2) *
     >              xNO3*ActCoeffs_Species(2,3)
        p_product=eql_constant * act_product
C
        WRITE(*,120) p_product
C
      ENDIF
C
C
      IF(Moles_Species(1,2).GT.0.D0 .AND. 
     >   Moles_Species(2,4).GT.0.D0 .AND.
     >   NH3_HCl) THEN
C        ----------------------------------------------
C       | ammonium and chloride are present. Calculate |
C       | and print the pressure product.              |
C       | NH4(aq) + Cl(aq) = NH3(g) + HCl(g)           |
C        ---------------------------------------------- 
C
C       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
C       !  H-NH4-Cl solutions not yet enabled for f(T) model !
C       !  Reaction ref. for NH4Cl sat. ratio in SATEQL call !
C       !  is 33 (SATEQL returns -99).                       !
C       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
C
        IF(print_header) THEN
C       ..header, if this is first
C         product to be printed.
          WRITE(*,100)
          print_header=.FALSE.
        ENDIF
C
        eql_constant=SATEQL(33,T)/(xKH_NH3(T)*xKH_HCl(T))
C                            ^-reaction ref for NH4Cl sat. ratio
        xNH4=Moles_Species(1,2)/Sum
        xCl=Moles_Species(2,4)/Sum
        act_product=xNH4*ActCoeffs_Species(1,2) *
     >              xCl*ActCoeffs_Species(2,4)
        p_product=eql_constant * act_product
C
        WRITE(*,130) p_product
C
      ENDIF

C     ******
      RETURN
C     ******
C
100   FORMAT(
     >/1X,'- Calculated partial pressure products (over liquid) -',
     >/1X,'Product (atm)',22X,'Reaction')
110   FORMAT(1X,'pNH3**2 * pH2SO4 =',E12.4E3,',',T37,
     >'2*NH4(aq) + SO4(aq) = 2*NH3(g) + H2SO4(g)')
120   FORMAT(1X,'pNH3 * pHNO3 =',E12.4E3,',',T37,
     >'NH4(aq) + NO3(aq) = NH3(g) + HNO3(g)')
130   FORMAT(1X,'pNH3 * pHCl =',E12.4E3,',',T37,
     >'NH4(aq) + Cl(aq) = NH3(g) + HCl(g)')
C
      END
C
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
      SUBROUTINE Charge_Balance(Molal,Test_Charge_Bal)
      IMPLICIT REAL*8 (A-H,O-Z)
C
C  *=============================================================*
C  * This routine tests for valid numbers of moles of ions, and  *
C  * for charge balance. For inputs that satisfy the test, but   *
C  * whose charge balance is not *exact*, the mole numbers are   *
C  * revised to make this so if the internal parameter 'revise'  *
C  * is set to '.true.'.                                         *
C  *=============================================================*
C
      PARAMETER(NCmax=3, NAmax=5)
C
      LOGICAL Revise
      PARAMETER(Revise = .TRUE., Tol_Charge=1.D-4)
C
C     ---- arguments ----
      LOGICAL Test_Charge_Bal
      REAL*8  Molal(-NAmax:NCmax)
C
      COMMON/CHARGE/CATchrg(NCmax),ANchrg(NAmax),
     >              Vca(NCmax,NAmax),Vac(NAmax,NCmax)
C
      Test_Charge_Bal=.FALSE.     
      sum_Cation=0.D0
      sum_Anion=0.D0
C
C      ------------------------------------- 
C     | Add up the cation and anion charges |
C      ------------------------------------- 
      DO 4 I=1,MAX(NCmax,NAmax)
        IF(I.LE.NCmax) 
     >    sum_Cation=sum_Cation+CATchrg(I)*Molal(I)
        IF(I.LE.NAmax)
     >    sum_Anion=sum_Anion+ANchrg(I)*Molal(-I)
4     CONTINUE
C
C
C      ----------------------------------------------------------
C     | Check for negative sums, allow a system with no ions to  |
C     | pass, and calculate fractional difference between cation |
C     | and anion sums.                                          |
C      ----------------------------------------------------------       
      IF(sum_Cation.LT.0.D0 .OR. sum_Anion.LT.0.D0) THEN
C     ..fail
        var=999.D0  
      ELSEIF(sum_Cation.EQ.0.D0 .AND. sum_Anion.EQ.0.D0) THEN
C     ..allow test to be passed if no cations or anions are present
        var=0.D0
      ELSE
C     ..obtain difference as fraction of total
        var=(sum_Cation-sum_Anion)/(sum_Cation+sum_Anion)
      ENDIF
C
C
C      ------------------------------------------ 
C     | Test fractional difference (var) against |
C     | requested tolerance.                     |
C     | If test is passed, but var NE 0, then    |
C     | rebalance the numbers of moles to        |
C     | attain exact agreement. This is needed   |
C     | so that the Gibbs minimiser can satisfy  |
C     | the charge balance constraint accurately.|
C      ------------------------------------------
      IF(ABS(var).LT.Tol_Charge) THEN
        Test_Charge_Bal=.TRUE.
C
        IF(Revise .AND. (var.NE.0.D0)) THEN
          IF(var.GT.0.D0) THEN
C        .. +ve charge is > -ve charge, so increase moles of -ve ions
C           proportionately to get agreement.
            ratio=sum_Cation/sum_Anion
            imax=NAmax
C         ..rebalance
            DO 5 I=1,imax
              Molal(-I)=Molal(-I) * ratio
5           CONTINUE
          ELSEIF(var.LT.0.D0) THEN
C        .. -ve charge is > +ve charge, so increase moles of +ve ions
C           proportionately, to get agreement.
            ratio=sum_Anion/sum_Cation
            imax=NCmax
C         ..rebalance
            DO 6 I=1,imax
              Molal(I)=Molal(I) * ratio
6           CONTINUE
          ENDIF
        ENDIF
C
      ENDIF
C
C     ******
      RETURN
C     ******
C       
      END
C
C +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
      FUNCTION Sum_Liquid_Moles(Moles_Species)
      IMPLICIT REAL*8 (A-H,O-Z)
C
      PARAMETER(NCmax=3, NAmax=5, NNmax=1)
      PARAMETER(NSpmax=28)
C
C     --- arguments ---
      REAL*8 Moles_Species(5,NSpmax)
C
      Sum_Liquid_Moles=0.D0
C
      DO 1 I=1,3
        IF(I.EQ.1) J_max=NCmax
        IF(I.EQ.2) J_max=NAmax
        IF(I.EQ.3) J_max=NNmax
        DO 2 J=1,J_max
          Sum_Liquid_Moles=Sum_Liquid_Moles + Moles_Species(I,J)
2       CONTINUE
1     CONTINUE
C
      END
C
C >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
C              end of aerotr.for code for web input/output
C <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
C >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
C <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
C
      SUBROUTINE fL1J1(TEMP,MOLAL,EPS,I1,L1,J1)
      IMPLICIT REAL*8 (A-H,O-Z)
C
      PARAMETER(NCmax=3, NAmax=5, NNmax=1, HREL=1.D-2)
C
      REAL*8 MOLAL(-NAmax:NCmax),MFRAC(-NAmax:NCmax),ACT(-NAmax:NCmax),
     >       mlocal(-NAmax:NCmax),TX(5),AWLG(5),L1,J1
C
C-------------------------------------------------------------------------
C     This routine calculates the partial molar enthlpy and heat capacity
C     of water.
C-------------------------------------------------------------------------
C
C
      IF(I1.NE.1) THEN
C     ..no calculations required
        L1=-99.D0
        J1=-99.D0
        RETURN
      ENDIF
C
      DO 2 I=-NAmax,NCmax
        mlocal(I)=MOLAL(I)
2     CONTINUE
C
      Tincr=HREL*TEMP
      TX(1)=TEMP + 2.D0*Tincr
      TX(2)=TEMP + Tincr
      TX(3)=TEMP 
      TX(4)=TEMP - Tincr
      TX(5)=TEMP - 2.D0*Tincr
    
      DO 3 J=1,5
        CALL SPEC(TX(J),mlocal,EPS,MFRAC,ACT,AW)
        AWLG(J)=LOG(AW)         
3     CONTINUE      
C
C       
      DLNAW=(-AWLG(1) + 8.D0*AWLG(2) - 8.D0*AWLG(4) 
     >       + AWLG(5))/(12.D0*Tincr) 
C       
      DDLNAW=(-AWLG(1) + 16.D0*AWLG(2) - 30.D0*AWLG(3)
     >        + 16.D0*AWLG(4) - AWLG(5))/(12.D0*Tincr**2)
C       
      L1 = -8.3144D0*TEMP**2 * DLNAW
      J1 = -8.3144D0*TEMP*(2.D0*DLNAW + TEMP*DDLNAW)
C
      RETURN
      END
C
C +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
      SUBROUTINE fLsJs(TEMP,MOLAL,EPS,Is,Ls,Js)
      IMPLICIT REAL*8 (A-H,O-Z)
C
      PARAMETER(NCmax=3, NAmax=5, NNmax=1, HREL=1.D-2)
C
      REAL*8 MOLAL(-NAmax:NCmax),MFRAC(-NAmax:NCmax),ACT(-NAmax:NCmax),
     >       mlocal(-NAmax:NCmax),TX(5),actSlg(5),Ls,Js
C
C-------------------------------------------------------------------------
C     This routine calculates the partial molar enthlpy and heat capacity
C     of a solute ion specified using the integer variable 'Is'.
C-------------------------------------------------------------------------
C
C
      IF(Is.LT.-NAmax .OR. Is.GT.NCmax .OR. Is.EQ.0) THEN
C     ..no calculations required
        Ls=-99.D0
        Js=-99.D0
        RETURN
      ENDIF
C
      DO 2 I=-NAmax,NCmax
        mlocal(I)=MOLAL(I)
2     CONTINUE
C
      Tincr=HREL*TEMP
      TX(1)=TEMP + 2.D0*Tincr
      TX(2)=TEMP + Tincr
      TX(3)=TEMP 
      TX(4)=TEMP - Tincr
      TX(5)=TEMP - 2.D0*Tincr
    
      DO 3 J=1,5
        CALL SPEC(TX(J),mlocal,EPS,MFRAC,ACT,AW)
        actSlg(J)=LOG(ACT(Is))         
3     CONTINUE      
C
C       
      DLNs=(-actSlg(1) + 8.D0*actSlg(2) - 8.D0*actSlg(4) 
     >       + actSlg(5))/(12.D0*Tincr) 
C       
      DDLNs=(-actSlg(1) + 16.D0*actSlg(2) - 30.D0*actSlg(3)
     >        + 16.D0*actSlg(4) - actSlg(5))/(12.D0*Tincr**2)
C       
      Ls = -8.3144D0*TEMP**2 * DLNs
      Js = -8.3144D0*TEMP*(2.D0*DLNs + TEMP*DDLNs)
C
      RETURN
      END
C
C +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
C >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
C              all the following code should be identical with the
C              routines present in the aerocalc program for the 
C              tropospheric/stratospheric model. (NB: a few redundant
C              declarations removed 18/1/98.)
C <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
C >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
C <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
C
C ===================== SLDS AND GASES (Gibbs program) ===================
C
C  Ref: Solid             Ref: Solid               Ref: Solid   
C    1  H2O                10  (NH4)2SO4            18  Na2SO4           
C    2  H2SO4              11  (NH4)3H(SO4)2        19  Na2SO4.10H2O     
C    3  H2SO4.H2O          12  NH4HSO4              20  Na3H(SO4)2     
C    4  H2SO4.2H2O         13  NH4NO3               21  NaHSO4.H2O     
C    5  H2SO4.3H2O         14  2NH4NO3.(NH4)2SO4    22  NaHSO4           
C    6  H2SO4.4H2O         15  3NH4NO3.(NH4)2SO4    23  NaH3(SO4)2.H2O   
C    7  H2SO4.6.5H2O       16  NH4NO3.NH4HSO4       24  Na2SO4.(NH4)2SO4.4H2O
C    8  HNO3.H2O           17  NH4Cl                25  NaNO3   
C    9  HNO3.3H2O                                   26  NaNO3.Na2SO4.H2O
C                                                   27  NaCl
C                                                   28  HCl.3H2O
C  Ref: Gas                                         29  (HBr) [a] 
C    1  H2O                                       
C    2  HNO3                                        [a] dummy solid to use as
C    3  HCl                                             electrolyte reference
C    4  NH3                                             number on input.
C    5  H2SO4
C    6  HBr
C
C
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
      BLOCK DATA
      IMPLICIT REAL*8 (A-H,O-Z)
C
C ********************************
C                                *
C Last updated: 9/5/98           *
C                                *
C ********************************
C
      PARAMETER(NCOEF=80)
      COMMON/PARAM/COEFF(NCOEF)
C
      PARAMETER(NCmax=3, NAmax=5, NNmax=1,
     >          multNN=NNmax**2,multNNN=NNmax**3,multCA=NCmax*NAmax, 
     >          multCCA=multCA*NCmax,multAAC=multCA*NAmax,
     >          multNCA=multCA*NNmax,multNCCA=multCCA*NNmax)
      PARAMETER(multNAAC=multAAC*NNmax,multCCAA=multCA**2,
     >          multNNca=multNN*multCA,multCCCA=multCCA*NCmax,
     >          multAAAC=multAAC*NAmax)
C      
      COMMON/SPEC1ES/NC,NA,NN,nCAT(NCmax),nAN(NAmax),nNEUT(NNmax)
      COMMON/COEFF1/RHO,A1ca(NCmax,NAmax),A2ca(NCmax,NAmax),
     >       B1ca(NCmax,NAmax),B2ca(NCmax,NAmax),
     >       Wnca(NNmax,NCmax,NAmax),Unca(NNmax,NCmax,NAmax),
     >       Vnca(NNmax,NCmax,NAmax),Wcca(NCmax,NCmax,NAmax),
     >       Waac(NAmax,NAmax,NCmax),Qncca(NNmax,NCmax,NCmax,NAmax),
     >       Qnaac(NNmax,NAmax,NAmax,NCmax),Ucca(NCmax,NCmax,NAmax),
     >       Uaac(NAmax,NAmax,NCmax),Wnn(NNmax,NNmax),Unn(NNmax,NNmax)
      COMMON/COEFF2/Ynnca(NNmax,NNmax,NCmax,NAmax),
     >              Xaaac(NAmax,NAmax,NAmax,NCmax),
     >              Xccca(NCmax,NCmax,NCmax,NAmax),
     >              Zccaa(NCmax,NCmax,NAmax,NAmax),
     >              Cnnn(NNmax,NNmax,NNmax)
      COMMON/CHARGE/CATchrg(NCmax),ANchrg(NAmax),
     >              Vca(NCmax,NAmax),Vac(NAmax,NCmax)
C
      INTEGER SOLID_R(1:30,4),SOLID_S(1:30,4),SOLID_N(1:30)
      COMMON/SLDS/SOLID_R,SOLID_S,SOLID_N
C
      DATA A1ca,A2ca,B1ca,B2ca,Wnca,Unca,Vnca/multCA*0.D0,multCA*0.D0,
     >     multCA*0.D0,multCA*0.D0,multNCA*0.D0,multNCA*0.D0,
     >     multNCA*0.D0/
      DATA Wcca,Waac,Ucca,Uaac,Qncca,Qnaac/multCCA*0.D0,multAAC*0.D0,
     >     multCCA*0.D0,multAAC*0.D0,multNCCA*0.D0,multNAAC*0.D0/
      DATA Wnn,Unn/multNN*0.D0,multNN*0.D0/
      DATA Ynnca,Xaaac,Xccca,Zccaa,Cnnn/multNNCA*0.D0,multAAAC*0.D0,
     >     multCCCA*0.D0,multCCAA*0.D0,multNNN*0.D0/
     >     
C
C   ..Initialize SPECIES, CHARGE and COEFF variables.
C
      DATA NC,NA,NN/3,5,1/
      DATA nCAT,nAN,nNEUT/1,2,3,1,2,3,4,5,1/
      DATA CATchrg,ANchrg/3*1.D0, 1.D0,2.D0,1.D0,1.D0,1.D0/
      DATA ((Vca(I,J),I=1,3),J=1,5)/3*0.5D0,3*0.666667D0,9*0.5D0/
      DATA ((Vac(J,I),J=1,5),I=1,3)/0.5D0,0.333333D0,3*0.5D0,
     >                              0.5D0,0.333333D0,3*0.5D0,
     >                              0.5D0,0.333333D0,3*0.5D0/
C
      DATA RHO/13.D0/
C
C   ..H2SO4:
      DATA (COEFF(I),I=1,48)/0.178334467D2,-0.625268629D1,0.591429323D0,
     > 0.223751841D0,2*0.D0,-0.998416390D1,0.348821776D0,
     >-0.119526172D-1,0.909425662D-2,0.149166944D-3,0.D0,
     >-0.143238371D1,-0.201636224D0,-0.443804232D-1,0.641847819D-2, 
     > 0.296327801D-3,0.D0,-0.207474566D1,0.594737744D0,
     > 0.674052221D-1,0.D0,-0.394845016D-3,0.D0,-0.982408701D2,
     >-0.205401806D2,-0.207137292D1,-0.376521937D-1,-0.139689758D-1,
     > 0.D0,-0.107752155D2,-0.879298257D0,-0.440528485D0,
     >-0.544913927D-1,-0.173541364D-3,0.D0,-0.133603464D2,
     >-0.459479578D1,-0.146220346D1,-0.157872023D0,-0.162230945D-3,
     > 0.D0,0.310121997D1,0.446189009D1,0.975254718D0,0.588748231D-2,
     >-0.901983372D-3,0.D0/
C
C   ..HNO3:
      DATA (COEFF(I),I=49,80)/0.2470883216D2,0.8537056072D1,
     > 0.D0,-0.5686459898D1,-0.2282434425D1,-0.2606309141D0,
     >-0.9318482814D-2,0.D0,-0.3301004522D1,0.1812617746D0, 
     > 0.1138792155D-1,0.3112011367D-2,2*0.D0,0.5981401888D-5,
     > 0.D0,-0.1009606757D1,-0.9191466033D-1,0.3513876226D-1, 
     > 2*0.D0,0.4419544700D-3,0.4469965190D-4,0.D0,0.1807488767D1,
     > 0.D0,-0.2615013798D0,-0.5744350823D-1,0.D0,0.5151608617D-3,
     > 2*0.D0/ 
C
C
C  Liquid Phase Species:
C
C  **Cations**           **Anions**             **Neutrals**
C  1: H                  -1: HSO4               0: H2O
C  2: NH4                -2: SO4
C                        -3: NO3
C                        -4: Cl
C                        -5: Br
C  Solid Phases:
C
C  Ref: Solid             Ref: Solid               Ref: Solid   
C    1  H2O                10  (NH4)2SO4            18  Na2SO4           
C    2  H2SO4              11  (NH4)3H(SO4)2        19  Na2SO4.10H2O     
C    3  H2SO4.H2O          12  NH4HSO4              20  Na3H(SO4)2     
C    4  H2SO4.2H2O         13  NH4NO3               21  NaHSO4.H2O     
C    5  H2SO4.3H2O         14  2NH4NO3.(NH4)2SO4    22  NaHSO4           
C    6  H2SO4.4H2O         15  3NH4NO3.(NH4)2SO4    23  NaH3(SO4)2.H2O   
C    7  H2SO4.6.5H2O       16  NH4NO3.NH4HSO4       24  Na2SO4.(NH4)2SO4.4H2O
C    8  HNO3.H2O           17  NH4Cl                25  NaNO3   
C    9  HNO3.3H2O                                   26  NaNO3.Na2SO4.H2O
C                                                   27  NaCl
C                                                   28  HCl.3H2O
C                                                   29  HBr
C
C
C  Initialise SOLID_R (reference numbers), SOLID_S (stoichiometries)
C  and SOLID_N (number of species associated with each solid phase)
C
C                                      1.ice   2.H2SO4     3.H2SO4.H2O
      DATA ((SOLID_R(I,J),J=1,4),I=1,30)/0,3*99, 1,-2,2*99,  1,-2,0,99,
C         4.H2SO4.2H2O         5.H2SO4.3H2O          6.H2SO4.4H2O 
     >      1,-2,0,99,           1,-2,0,99,            1,-2,0,99,
C         7.H2SO4.6.5H2O       8.HNO3.H2O            9.HNO3.3H2O
     >      1,-2,0,99,           1,-3,0,99,            1,-3,0,99,
C                                                   
C        10.NH42SO4           11.letovicite         12.NH4HSO4
     >      2,-2,2*99,           2,1,-2,99,            2,1,-2,99,
C        13.NH4NO3            14.2NH4NO3.(NH4)2SO4  15.3NH4NO3.(NH4)2SO4
     >      2,-3,2*99,           2,-3,-2,99,           2,-3,-2,99,
C        16.NH4NO3.NH4HSO4    17.NH4Cl              
     >      2,1,-3,-2,           2,-4,2*99,         
C                                                   
C        18.Na2SO4            19.Na2SO4.10H2O       20.Na3H(SO4)2
     >      3,-2,2*99,           3,-2,0,99,            3,1,-2,99,
C        21.NaHSO4.H2O        22.NaHSO4             23.NaH3(SO4)2.H2O
     >      3,1,-2,0,            3,1,-2,99,            3,1,-2,0,
C    24.Na2SO4.(NH4)2SO4.4H2O 25.NaNO3              26.NaNO3.Na2SO4.H2O
     >      3,2,-2,0,            3,-3,2*99,            3,-2,-3,0,
C        27.NaCl              28.HCl.3H2O           29.HBr
     >      3,-4,2*99,           1,-4,0,99,            1,-5,2*99,
     >      4*99/
C
C     NB: in the stoichiometries below H2SO4.6.5H2O is treated as 
C         2(H2SO4).13H2O. The equilibrium constant has been altered
C         to reflect this.
C
C
C                                        1.ice   2.H2SO4     3.H2SO4.H2O
      DATA ((SOLID_S(I,J),J=1,4),I=1,30)/1,3*0,    2,1,0,0,    2,1,1,0,
C         4.H2SO4.2H2O         5.H2SO4.3H2O          6.H2SO4.4H2O 
     >      2,1,2,0,             2,1,3,0,              2,1,4,0,
C         7.H2SO4.6.5H2O       8.HNO3.H2O            9.HNO3.3H2O
     >      4,2,13,0,            1,1,1,0,              1,1,3,0,
C                                                   
C        10.NH42SO4           11.letovicite         12.NH4HSO4
     >      2,1,2*0,             3,1,2,0,              1,1,1,0,
C        13.NH4NO3            14.2NH4NO3.(NH4)2SO4  15.3NH4NO3.(NH4)2SO4
     >      1,1,2*0,             4,2,1,0,              5,3,1,0,
C        16.NH4NO3.NH4HSO4    17.NH4Cl              
     >      2,1,1,1,             1,1,2*0,         
C                                                   
C        18.Na2SO4            19.Na2SO4.10H2O       20.Na3H(SO4)2
     >      2,1,0,0,             2,1,10,0,             3,1,2,0,
C        21.NaHSO4.H2O        22.NaHSO4             23.NaH3(SO4)2.H2O
     >      1,1,1,1,             1,1,1,0,              1,3,2,1,
C    24.Na2SO4.(NH4)2SO4.4H2O 25.NaNO3              26.NaNO3.Na2SO4.H2O
     >      2,2,2,4,             1,1,0,0,              3,1,1,1,
C        27.NaCl              28.HCl.3H2O           29.HBr 
     >      1,1,0,0,             1,1,3,0,              1,1,0,0,
     >      4*99/
C
C
C
C                                 1.ice   2.H2SO4     3.H2SO4.H2O
      DATA (SOLID_N(I),I=1,30)/     1,        2,           3,    
C         4.H2SO4.2H2O         5.H2SO4.3H2O          6.H2SO4.4H2O 
     >         3,                   3,                    3,         
C         7.H2SO4.6.5H2O       8.HNO3.H2O            9.HNO3.3H2O
     >         3,                   3,                    3,       
C                                                   
C        10.NH42SO4           11.letovicite         12.NH4HSO4
     >        2,                    3,                    3,    
C        13.NH4NO3            14.2NH4NO3.(NH4)2SO4  15.3NH4NO3.(NH4)2SO4
     >        2,                    3,                    3,       
C        16.NH4NO3.NH4HSO4    17.NH4Cl              
     >        4,                    2,             
C                                                   
C        18.Na2SO4            19.Na2SO4.10H2O       20.Na3H(SO4)2
     >        2,                    3,                    3,       
C        21.NaHSO4.H2O        22.NaHSO4             23.NaH3(SO4)2.H2O
     >        4,                    3,                    4,       
C    24.Na2SO4.(NH4)2SO4.4H2O 25.NaNO3              26.NaNO3.Na2SO4.H2O
     >        4,                    2,                    4,       
C        27.NaCl              28.HCl.3H2O           29.HBr 
     >        2,                    3,                    2,
     >        0/
C
      END
C
C +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
      FUNCTION SATEQL(nRef,T)
      IMPLICIT REAL*8 (A-H,O-Z)
C
      REAL*8 lnKTr
C
C Calculation of eql. const. for 
C saturation wrt solid phases
C
C ********************************
C                                *
C Last updated: 3/11/97          *
C                                *
C ********************************
C
C
       IF(nREF.EQ.1) THEN
C      ..ice:
         D=273.15D0-T
         DUM=-4.2091D-3*D-0.2152D-5*D**2+0.3233D-7*D**3
     >       +0.3446D-9*D**4+0.1758D-11*D**5+0.765D-14*D**6
         SATEQL=10.D0**DUM
       ELSEIF(nREF.EQ.2) THEN
C      ..H2SO4 (not implemented):
         SATEQL= -99.D0
       ELSEIF(nREF.EQ.3) THEN
C      ..H2SO4.H2O:
         SATEQL=EXP(66.92217-6434.0498/T-0.150405*T)
       ELSEIF(nREF.EQ.4) THEN
C      ..H2SO4.2H2O:
         SATEQL=EXP(0.821887D0-316.594/T)   
       ELSEIF(nREF.EQ.5) THEN
C      ..H2SO4.3H2O:
         SATEQL=EXP(47.45356-6758.7767/T-0.097946*T)
       ELSEIF(nREF.EQ.6) THEN
C      ..H2SO4.4H2O:
         SATEQL=EXP(68.85084-1.014909D4/T-0.140118*T)
       ELSEIF(nREF.EQ.7) THEN
C      ..H2SO4.6.5H2O [as 2(H2SO4).13H2O)]:
         SATEQL=EXP(2*(13.10788-5167.4839/T))
       ELSEIF(nRef.EQ.8) THEN
C      ..HNO3.H2O:
         SATEQL=EXP(29.20239D0-4061.5108D0/T-0.045788D0*T)
       ELSEIF(nREF.EQ.9) THEN
C      ..HNO3.3H2O:
         SATEQL=EXP(10.275284D0-3408.14402D0/T)
       ELSEIF(nREF.EQ.28) THEN
C      ..HCl.3H2O:
         SATEQL=EXP(19.4767D0-2478.44D0/T-0.033934D0*T)
       ELSEIF(nREF.EQ.10) THEN
C      ..(NH4)2SO4 (1995 paper with Chak):
         Tr=298.15D0
         lnKTr=-11.9603D0
         deltaH=6.08393D3
         deltaA=1640.84D0
         deltaB=-6.39321D0
         deltaC=0.
         SATEQL=EQLCON(lnKTr,deltaH,deltaA,deltaB,deltaC,Tr,T)
       ELSEIF(nREF.EQ.11) THEN
C      ..(NH4)3H(SO4)2 (May/June '97 revisions):
         Tr=298.15D0
         lnKTr=-26.0074487D0
         deltaH=-5.3159354331D3
         deltaA=-630.401191D0
         deltaB=0.
         deltaC=0.
         SATEQL=EQLCON(lnKTr,deltaH,deltaA,deltaB,deltaC,Tr,T)
       ELSEIF(nREF.EQ.12) THEN
C      ..NH4HSO4 (May/June '97 revisions):
         Tr=298.15D0
         lnKTr=-11.408434D0
         deltaH=-15.1675681455D3
         deltaA=-242.281346D0
         deltaB=0.
         deltaC=0.
         SATEQL=EQLCON(lnKTr,deltaH,deltaA,deltaB,deltaC,Tr,T)
       ELSEIF(nREF.EQ.13) THEN
C      ..NH4NO3 (May/June '97 revisions):
         IF(T.GE.256.2D0 .AND. T.LE.305.38D0) THEN
           Tr=298.15D0
           lnKTr=-5.52949D0
           deltaH=25.69D3
           deltaA=-101.23D0
           deltaB=0.
           deltaC=0.
         ELSEIF(T.LT.256.2D0) THEN
           Tr=256.2D0
           lnKTr=-7.37371D0
           deltaH=30.41D3
           deltaA=-101.23D0
           deltaB=0.
           deltaC=0.
         ELSEIF(T.GT.305.38D0) THEN
           Tr=305.38D0
           lnKTr=-5.28760D0
           deltaH=23.26D3
           deltaA=-101.23D0
           deltaB=0.
           deltaC=0.
         ENDIF
         SATEQL=EQLCON(lnKTr,deltaH,deltaA,deltaB,deltaC,Tr,T)
       ELSEIF(nREF.EQ.14) THEN
C      ..2NH4NO3.(NH4)2SO4 (May/June '97 revisions):
         Tr=298.15D0
         lnKTr=-23.681D0
         deltaH=58.845D3
         deltaA=0.
         deltaB=0.
         deltaC=0.
         SATEQL=EQLCON(lnKTr,deltaH,deltaA,deltaB,deltaC,Tr,T)
       ELSEIF(nREF.EQ.15) THEN
C      ..3NH4NO3.(NH4)2SO4 (May/June '97 revisions):
         Tr=298.15D0
         lnKTr=-29.422D0
         deltaH=84.860D3
         deltaA=0.
         deltaB=0.
         deltaC=0.
         SATEQL=EQLCON(lnKTr,deltaH,deltaA,deltaB,deltaC,Tr,T)
       ELSEIF(nREF.EQ.16) THEN
C      ..NH4NO3.NH4HSO4 (May/June '97 revisions):
         Tr=298.15D0
         lnKTr=-18.0814D0
         deltaH=8.730D3
         deltaA=0.
         deltaB=0.
         deltaC=0.
         SATEQL=EQLCON(lnKTr,deltaH,deltaA,deltaB,deltaC,Tr,T)
       ELSEIF(nREF.EQ.17) THEN
C      ..NH4Cl [not implemented]:
         SATEQL= -99.D0
       ELSEIF(nREF.EQ.18) THEN
C      ..Na2SO4 [not implemented]:
         SATEQL= -99.D0
       ELSEIF(nREF.EQ.19) THEN
C      ..Na2SO4.10H2O [not implemented]:
         SATEQL= -99.D0
       ELSEIF(nREF.EQ.25) THEN
C      ..NaNO3 [not implemented]:
         SATEQL= -99.D0
       ELSEIF(nREF.EQ.26) THEN
C      ..Na2SO4.NaNO3.H2O [not implemented]:
         SATEQL= -99.D0
       ELSEIF(nREF.EQ.27) THEN
C      ..NaCl [not implemented]:
         SATEQL= -99.D0
       ELSEIF(nREF.EQ.20) THEN
C      ..Na3H(SO4)2 [not implemented]:
         SATEQL= -99.D0
       ELSEIF(nREF.EQ.21) THEN
C      ..NaHSO4.H2O  [not implemented]:
         SATEQL= -99.D0
       ELSEIF(nREF.EQ.22) THEN
C      ..NaHSO4 [not implemented]:
         SATEQL= -99.D0
       ELSEIF(nREF.EQ.23) THEN
C      ..NaH3(SO4)2.H2O [not implemented]:
         SATEQL= -99.D0
       ELSEIF(nREF.EQ.24) THEN
C      ..Na2SO4.(NH4)2SO4.4H2O [not implemented]:
         SATEQL= -99.D0
       ELSEIF(nREF.EQ.29) THEN
C      ..HBr [not implemented, as this 'solid' is just a 
C        dummy for input as an electrolyte]:
         SATEQL= -99.D0
       ENDIF
C
       RETURN
       END
C
C
C +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
         FUNCTION EQLCON(lnKTr,deltaH,deltaA,deltaB,deltaC,Tr,T)
         IMPLICIT REAL*8 (A-H,O-Z)
C
         REAL*8 lnKTr
C
         R=8.3144D0
C
         dum= lnKTr + deltaH/R * (1/Tr-1/T)
     >       +deltaA/R*(Tr/T - 1 + LOG(T/Tr))
     >       +deltaB/(2*R)*(Tr*(Tr/T-1) + T - Tr)
     >       +deltaC/(6*R)*(2*Tr**2*(Tr/T-1) + T**2 - Tr**2)
C
         EQLCON=EXP(dum)
C
         END
C
C +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
      SUBROUTINE PARASET(T)
      IMPLICIT REAL*8 (A-H,O-Z)
C
      PARAMETER(NCOEF=80, TR1=328.15D0)
      PARAMETER(NCmax=3,NAmax=5,NNmax=1)
C
      REAL*8 COEFF(NCOEF)
      COMMON/PARAM/COEFF
C
C
C**********************************************************************
C                                                                     *
C     Last changed: 9/5/98                                            *
C                                                                     *
C**********************************************************************
C Elements required in main program unit follow:                      *
C**********************************************************************
C
      COMMON/COEFF1/RHO,A1ca(NCmax,NAmax),A2ca(NCmax,NAmax),
     >       B1ca(NCmax,NAmax),B2ca(NCmax,NAmax),
     >       Wnca(NNmax,NCmax,NAmax),Unca(NNmax,NCmax,NAmax),
     >       Vnca(NNmax,NCmax,NAmax),Wcca(NCmax,NCmax,NAmax),
     >       Waac(NAmax,NAmax,NCmax),Qncca(NNmax,NCmax,NCmax,NAmax),
     >       Qnaac(NNmax,NAmax,NAmax,NCmax),Ucca(NCmax,NCmax,NAmax),
     >       Uaac(NAmax,NAmax,NCmax),Wnn(NNmax,NNmax),Unn(NNmax,NNmax)
      COMMON/COEFF2/Ynnca(NNmax,NNmax,NCmax,NAmax),
     >              Xaaac(NAmax,NAmax,NAmax,NCmax),
     >              Xccca(NCmax,NCmax,NCmax,NAmax),
     >              Zccaa(NCmax,NCmax,NAmax,NAmax),
     >              Cnnn(NNmax,NNmax,NNmax)
C
C
C ====================================================================
C ====================================================================
C    SECTION (1) FOR H2SO4, HNO3. PARAMETERS FROM THE CARSLAW 
C    ET AL 1995 PAPER:
C
      D1=T-TR1
C
C.....H-HSO4...
C
      A1ca(1,1)=17.D0
      B1ca(1,1)=COEFF(1) + D1*(COEFF(2)*1.D-1 + 
     >                    D1*(0.5D0*COEFF(3)*1.D-2 +
     >                        D1*(COEFF(4)*1.D-3/6.D0 + 
     >                            D1*(COEFF(5)*1.D-3/12.D0 +
     >                                D1*COEFF(6)*1.D-3/20.D0))))
      Wnca(1,1,1)=COEFF(7) + D1*(COEFF(8)*1.D-1 + 
     >                    D1*(0.5D0*COEFF(9)*1.D-2 +
     >                        D1*(COEFF(10)*1.D-3/6.D0 + 
     >                            D1*(COEFF(11)*1.D-3/12.D0 +
     >                                D1*COEFF(12)*1.D-3/20.D0))))
      Unca(1,1,1)=COEFF(13) + D1*(COEFF(14)*1.D-1 + 
     >                    D1*(0.5D0*COEFF(15)*1.D-2 +
     >                        D1*(COEFF(16)*1.D-3/6.D0 + 
     >                            D1*(COEFF(17)*1.D-3/12.D0 +
     >                                D1*COEFF(18)*1.D-3/20.D0))))
      Vnca(1,1,1)=COEFF(19) + D1*(COEFF(20)*1.D-1 + 
     >                    D1*(0.5D0*COEFF(21)*1.D-2 +
     >                        D1*(COEFF(22)*1.D-3/6.D0 + 
     >                            D1*(COEFF(23)*1.D-3/12.D0 +
     >                                D1*COEFF(24)*1.D-3/20.D0))))
C
C.....H-SO4...
C
      A1ca(1,2)=9.5D0
      B1ca(1,2)=COEFF(25) + D1*(COEFF(26)*1.D-1 + 
     >                    D1*(0.5D0*COEFF(27)*1.D-2 +
     >                        D1*(COEFF(28)*1.D-3/6.D0 + 
     >                            D1*(COEFF(29)*1.D-3/12.D0 +
     >                                D1*COEFF(30)*1.D-3/20.D0))))
      Wnca(1,1,2)=COEFF(31) + D1*(COEFF(32)*1.D-1 + 
     >                    D1*(0.5D0*COEFF(33)*1.D-2 +
     >                        D1*(COEFF(34)*1.D-3/6.D0 + 
     >                            D1*(COEFF(35)*1.D-3/12.D0 +
     >                                D1*COEFF(36)*1.D-3/20.D0))))
      Unca(1,1,2)=COEFF(37) + D1*(COEFF(38)*1.D-1 + 
     >                    D1*(0.5D0*COEFF(39)*1.D-2 +
     >                        D1*(COEFF(40)*1.D-3/6.D0 + 
     >                            D1*(COEFF(41)*1.D-3/12.D0 +
     >                                D1*COEFF(42)*1.D-3/20.D0))))
      Vnca(1,1,2)=COEFF(43) + D1*(COEFF(44)*1.D-1 + 
     >                    D1*(0.5D0*COEFF(45)*1.D-2 +
     >                        D1*(COEFF(46)*1.D-3/6.D0 + 
     >                            D1*(COEFF(47)*1.D-3/12.D0 +
     >                                D1*COEFF(48)*1.D-3/20.D0))))
C
C.....HNO3...
C
      D4=T-330.D0
      fn1=1.D-1
      fn2=1.D-2/2.D0
      fn3=1.D-3/6.D0
      fn4=1.D-4/12.D0
      fn5=1.D-5/20.D0
      fn6=1.D-6/30.D0
      fn7=1.D-7/42.D0
C      
      A1ca(1,3)=29.03169775D0
      B1ca(1,3)=coeff(49) + d4*(coeff(50)*fn1 +
     >                d4*(coeff(51)*fn2 +
     >                 d4*(coeff(52)*fn3 +
     >                  d4*(coeff(53)*fn4 +
     >                   d4*(coeff(54)*fn5 +
     >                    d4*(coeff(55)*fn6 +
     >                     d4*coeff(56)*fn7))))))
C
      Wnca(1,1,3)=coeff(57) + d4*(coeff(58)*fn1 +
     >                d4*(coeff(59)*fn2 +
     >                 d4*(coeff(60)*fn3 +
     >                  d4*(coeff(61)*fn4 +
     >                   d4*(coeff(61)*fn5 +
     >                    d4*(coeff(63)*fn6 +
     >                     d4*coeff(64)*fn7))))))
C
      Unca(1,1,3)=coeff(65)+d4*(coeff(66)*fn1 +
     >                d4*(coeff(67)*fn2 +
     >                 d4*(coeff(68)*fn3 +
     >                  d4*(coeff(69)*fn4 +
     >                   d4*(coeff(70)*fn5 +
     >                    d4*(coeff(71)*fn6 +
     >                     d4*coeff(72)*fn7))))))
C
      Vnca(1,1,3)=coeff(73)+d4*(coeff(74)*fn1 +
     >                d4*(coeff(75)*fn2 +
     >                 d4*(coeff(76)*fn3 +
     >                  d4*(coeff(77)*fn4 +
     >                   d4*(coeff(78)*fn5 +
     >                    d4*(coeff(79)*fn6 +
     >                     d4*coeff(80)*fn7))))))
C
C
C
C  ** H-Cl **
C
      D4=T-330.D0
      fn1=1.D-1
      fn2=1.D-2/2.D0
      fn3=1.D-3/6.D0
      fn4=1.D-4/12.D0
      fn5=1.D-5/20.D0
      fn6=1.D-6/30.D0
      fn7=1.D-7/42.D0
C
      A1ca(1,4)=2.43288D0
      B1ca(1,4)=0.637130453D2 + d4*(-0.127961571D1*fn1 +
     >                d4*(-0.855184587D0*fn2 +
     >                 d4*(-0.483462259D0*fn3 +
     >                  d4*(-0.769606918D-1*fn4 +
     >                   d4*(0.227428690D-3*fn5 +
     >                    d4*(0.D0*fn6 +
     >                     d4*0.D0*fn7))))))
      Wnca(1,1,4)=0.D0 + d4*(0.712414208D0*fn1 +
     >                d4*(0.D0*fn2 +
     >                 d4*(-0.253763077D-1*fn3 +
     >                  d4*(-0.725301184D-2*fn4 +
     >                   d4*(0.D0*fn5 +
     >                    d4*(0.D0*fn6 +
     >                     d4*0.D0*fn7))))))
      Unca(1,1,4)=0.757587583D1 + d4*(0.121792108D1*fn1 +
     >                d4*(0.230535499D0*fn2 +
     >                 d4*(0.576304748D-1*fn3 +
     >                  d4*(0.D0*fn4 +
     >                   d4*(0.D0*fn5 +
     >                    d4*(0.D0*fn6 +
     >                     d4*0.D0*fn7))))))
      Vnca(1,1,4)=-0.153591245D2 + d4*(-0.389164072D0*fn1 +
     >                d4*(-0.891045013D-1*fn2 +
     >                 d4*(0.D0*fn3 +
     >                  d4*(0.654055055D-2*fn4 +
     >                   d4*(0.D0*fn5 +
     >                    d4*(0.D0*fn6 +
     >                     d4*0.D0*fn7))))))
C
C
C   ** H - Br ** (added 9/5/98)
      A1ca(1,5) = 6.19176D0
      B1ca(1,5) =     0.327847880D2  + 
     >               d4*(0.919647588D0    *fn1 +
     >               d4*(0.742200662D0    *fn2 +
     >               d4*(0.263972123D0    *fn3 +
     >               d4*(0.D0             *fn4 +
     >               d4*(0.616206313D-3   *fn5 +
     >               d4*(0.D0             *fn6 +
     >               d4* 0.D0             *fn7))))))
      Wnca(1,1,5) =  -0.131504151D2  + 
     >               d4*(0.116060928D1    *fn1 +
     >               d4*(0.312135679D0    *fn2 +
     >               d4*(0.146033030D0    *fn3 +
     >               d4*(0.708445887D-2   *fn4 +
     >               d4*(0.D0             *fn5 +
     >               d4*(0.D0             *fn6 +
     >               d4* 0.D0             *fn7))))))
      Unca(1,1,5) =  -0.699721143D1  +
     >               d4*(0.215281904D1    *fn1 +
     >               d4*(0.902856022D0    *fn2 +
     >               d4*(0.410269215D0    *fn3 +
     >               d4*(0.214487024D-1   *fn4 +
     >               d4*(0.D0             *fn5 +
     >               d4*(0.D0             *fn6 +
     >               d4* 0.D0             *fn7))))))
      Vnca(1,1,5) =  -0.414565596D1  + 
     >               d4*(-0.124431743D1   *fn1 +
     >               d4*(-0.748579644D0   *fn2 +
     >               d4*(-0.319716995D0   *fn3 +
     >               d4*(-0.161436912D-1  *fn4 +
     >               d4*(0.D0             *fn5 +
     >               d4*(0.D0             *fn6 +
     >               d4* 0.D0             *fn7))))))
C
C
C.....HSO4-NO3-H interactions:
      Waac(1,3,1)=-4.280D0
      Uaac(1,3,1)=0.201362D0 + 0.084830D0*(T-273.15D0)
      Waac(3,1,1)=Waac(1,3,1)
      Uaac(3,1,1)=-Uaac(1,3,1)
C
C.....SO4-NO3-H interactions:
      Waac(2,3,1)=-0.033291*(T-273.15D0)
      Waac(3,2,1)=Waac(2,3,1)
C
C.....HSO4-Br-H interactions:
      Waac(1,5,1)=-0.798159D0
      Waac(5,1,1)=Waac(1,5,1)
C
C.....SO4-Br-H interactions:
      Waac(2,5,1)=13.6034584D0 
     >           + 1.D3*39.72218D0*(1.D0/T - 1.D0/298.15D0)
      Waac(5,2,1)=Waac(2,5,1)
C
C
C ====================================================================
C ====================================================================
C    SECTION (2), other model elements.
C
C..NH4-HSO4.. (June '97 revisions)
C
      D1=T-298.15D0
C
      A1ca(2,1)=13.D0
      B1ca(2,1)=0.495871301D+01 + D1*(0.4D0*1.D-1) 
      Wnca(1,2,1)=-0.163376858D+01 + D1*(-0.539327161D-1*1.D-1) 
      Unca(1,2,1)=0.362701710D+00
      Vnca(1,2,1)=0.367166597D+00 + D1*(-0.74241954D-1*1.D-1) 
C
C..NH4-SO4..(Values from 1995 paper with Chak)
C
      TR=298.15D0
      A1ca(2,2)=13.D0
      B1ca(2,2)=  0.1399385D+02 + (T-TR)*(0.D0 - TR*(-0.02103172D0))
     >              + 0.5D0*(T**2 - TR**2)*(-0.02103172D0)
C
      A2ca(2,2)=1.5D0
      B2ca(2,2)= -0.1713243D2 + (T-TR)*(0.8461758D0-TR*(-0.9880561D-3))
     >              + 0.5D0*(T**2 - TR**2)*(-0.9880561D-3)
C
      Wnca(1,2,2)=-0.1904921D1 + (T-TR)*(0.07911813D0-TR*(0.1173141D-3))
     >              + 0.5D0*(T**2 - TR**2)*(0.1173141D-3)
C
      Unca(1,2,2)= 0.2125957D1 + (T-TR)*(-7.160173D-3-TR*(0.3064704D-3))
     >              + 0.5D0*(T**2 - TR**2)*(0.3064704D-3)
C
      Vnca(1,2,2)= -0.2291087D1+(T-TR)*(-0.03772486D0-TR*(-1.865885D-4))
     >              + 0.5D0*(T**2 - TR**2)*(-1.865885D-4)
C
C
C
C   ..Now NH4NO3 (298.15K & thermal parameters revised May '97)
C
      TR=298.15D0
      BL1=0.283192D0
      BL2=-0.0716627D0
      WL= -0.00696723D0
      UL=  0.171736D-2 
      VL=  0.221706D-2 
C
      BJ1= -0.352077D-02
      BJ2= 0.D0
      WJ=  0.489933D-5 
      UJ=  0.D0
      VJ= -0.359406D-4 
C
      d2B1=BJ1 - (2/TR)*BL1
      d2B2=BJ2 - (2/TR)*BL2
      d2W=WJ - (2/TR)*WL
      d2U=UJ - (2/TR)*UL
      d2V=VJ - (2/TR)*VL
C
      A1ca(2,3)=7.D0
      A2ca(2,3)=13.D0
      B1ca(2,3)=0.130466D2+(T-TR)*(BL1-TR*d2B1)+0.5D0*(T**2-TR**2)*d2B1 
      B2ca(2,3)=-0.162254D2+(T-TR)*(BL2-TR*d2B2)+0.5D0*(T**2-TR**2)*d2B2 
      Wnca(1,2,3)=0.616136D0+(T-TR)*(WL-TR*d2W)+0.5D0*(T**2-TR**2)*d2W
      Unca(1,2,3)=-0.403564D-1+(T-TR)*(UL-TR*d2U)+0.5D0*(T**2-TR**2)*d2U
      Vnca(1,2,3)=-0.680507D0+(T-TR)*(VL-TR*d2V)+0.5D0*(T**2-TR**2)*d2V 
C
C
C   ..Now the mixture parameters:
C
C     ..HSO4-SO4-NH4 interactions (revised June '97):
        D1=T-298.15D0
        Waac(1,2,2)=-8.17099584D0
        Uaac(1,2,2)=-11.3871516D0 + D1*(0.633515474D0*1.D-1) 
        Waac(2,1,2)=Waac(1,2,2)
        Uaac(2,1,2)=-Uaac(1,2,2)
C
C
C     ..H-NH4-HSO4 interactions (revised June '97):
        D1=T-298.15D0
        Wcca(1,2,1)=-16.4807039D0 + D1*(0.494791553D0*1.D-1) 
        Qncca(1,1,2,1)=7.084464D0 
        Wcca(2,1,1)=Wcca(1,2,1)
        Qncca(1,2,1,1)=Qncca(1,1,2,1)
C
C
C     ..H-NH4-SO4 interactions (revised June '97):
        Wcca(1,2,2)=-8.33078088D0
        Qncca(1,1,2,2)=2.60624775D0
        Wcca(2,1,2)=Wcca(1,2,2)
        Qncca(1,2,1,2)=Qncca(1,1,2,2)
C
C     ..H-NH4-NO3 interactions (revised May and June '97):
        Wcca(1,2,3)=-0.40728D1 + 0.42031D-01*(T-298.15D0)
     >              -0.32245D-03*(T-298.15D0)**2
        Qncca(1,1,2,3)=0.75481D0
C     ..note NH4 first, as this corresponds to order in fitting program.
        Ucca(2,1,3)=0.10055D1
        Wcca(2,1,3)    = Wcca(1,2,3)
        Qncca(1,2,1,3) = Qncca(1,1,2,3)
        Ucca(1,2,3)    = -Ucca(2,1,3)
C
C     ..SO4-NO3-NH4 interactions (revised 2 June '97):
        Waac(2,3,2)=0.32359D0+0.16567D0*(T-298.15D0)
        Qnaac(1,2,3,2)=2.6800D0-0.13322D0*(T-298.15D0)
C     ..sign reversal below as NO3 came before SO4 in the NH4-NO3-SO4 fit.
        Uaac(2,3,2)=-0.31582D0
        Waac(3,2,2)    = Waac(2,3,2)
        Qnaac(1,3,2,2) = Qnaac(1,2,3,2)
        Uaac(3,2,2)    =-Uaac(2,3,2)
C
C     ..HSO4-NO3-NH4 interactions (June '97):
        Waac(1,3,2)=-3.07081D0
        Waac(3,1,2)=Waac(1,3,2)
C
C
      RETURN
      END
C
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
      SUBROUTINE SPEC(T,MOLAL,EPS,MFRAC,ACT,AW)
      IMPLICIT REAL*8 (A-H,O-Z)
C
      PARAMETER(NCmax=3,NAmax=5,NNmax=1)
C
      REAL*8 MOLAL(-NAmax:NCmax),MFRAC(-NAmax:NCmax),ACT(-NAmax:NCmax)
      INTEGER iACTC(NCmax),iACTA(NAmax),iACTN(NNmax)
      REAL*8 mHtotal,mSO4total,LOW
      REAL*8 xCAT(NCmax),xAN(NAmax),xNEUT(NNmax),
     >       actC(NCmax),actA(NAmax),actN(NNmax)
C
C   ..set to 1 to ensure all activity coeffs are calculated
      DATA iACTA/NAmax*1/iACTC/NCmax*1/
C
C
C   ..set the model parameters for the present temperature:
C
      CALL PARASET(T)
C
C
C   ..if there is no H+ and no SO42- in the system, calculate species mole 
C     fractions and obtain activity coefficients by direct calls to the 
C     routines:
C
C   ..Checking for presence of both H+ and SO42-:
      IF(MOLAL(-1).EQ.0.D0 .AND. 
     >   (MOLAL(1).EQ.0.D0 .OR. MOLAL(-2).EQ.0.D0)) THEN
C     ..no bisulphate, and either H+ or SO42- also absent:
        iSpec=0
      ELSE
C     ..we must calculate the equilibrium:
        iSpec=1
C       ensure all HSO4- is converted to H+ and SO42-:
        MOLAL(-2)=MOLAL(-2)+MOLAL(-1)
        MOLAL(1)=MOLAL(1)+MOLAL(-1)
        MOLAL(-1)=0.D0
      ENDIF
C
      IF(Ispec.EQ.0) THEN
C     ..No speciation to calculate      
        SUM=0.D0
        DO 1 I=-NAmax,NCmax
          IF(I.NE.0) SUM=SUM + MOLAL(I)
1       CONTINUE
        var=55.508681D0+SUM
C
        DO 2 I=-NAmax,NCmax
          IF(I.LT.0) xAN(-I)=MOLAL(I)/var
          IF(I.EQ.0) xNEUT(1)=55.508681D0/var
          IF(I.GT.0) xCAT(I)=MOLAL(I)/var
2       CONTINUE
C
      ELSE
C
C     ..set range of initial estimates for free H+ (HSO4 = SO4 + H)
C
        HIGH=MOLAL(1)
        mHtotal=MOLAL(1)
        mSO4total=MOLAL(-2)
        var=mHtotal - mSO4total
        IF(var .GE. 0.D0) THEN
          LOW=VAR
        ELSE
          LOW=0.D0
        ENDIF
C
C     ..determine sulphate/bisulphate speciation iteratively:
        CALL ZBRENT(LOW,HIGH,mHtotal,mSO4total,MOLAL,T,EPS,xCAT,xAN,
     >              xNEUT,ACT)
C
      ENDIF
C
C     Calculate the activity coefficients, and water activity:
      DO 3 I=1,NCmax
        IF(Ispec.EQ.1 .AND. I.EQ.1) THEN
C       ..we already have the activity coefficient of H+
          CONTINUE
        ELSE
          CALL nCATION(xCAT,xAN,xNEUT,T,ACT(I),I,iACTC,actC,.FALSE.)
        ENDIF
3     CONTINUE
C
      DO 4 J=-1,-NAmax,-1
        IF(Ispec.EQ.1 .AND. J.GE.-2) THEN
C       ..we already have the activity coefficients of HSO4 & SO4
          CONTINUE
        ELSE
          CALL nANION(xCAT,xAN,xNEUT,T,ACT(J),-J,iACTA,actA,.FALSE.)
        ENDIF
4     CONTINUE
C
      CALL nNEUTRAL(xCAT,xAN,xNEUT,T,ACT(0),1,iACTN,actN,AW,G,
     >              .FALSE.)
C
C
C   ..assign values of the mole fractions to be returned to the main
C     program:
      DO 5 I=-NAmax,NCmax
        IF(I.LT.0) MFRAC(I)=xAN(-I)
        IF(I.EQ.0) MFRAC(I)=xNEUT(1)
        IF(I.GT.0) MFRAC(I)=xCAT(I)
5     CONTINUE
C
      RETURN
      END
C
C +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
      SUBROUTINE SALTS(NS,mSalt,refSalt,MOLAL)
      IMPLICIT REAL*8 (A-H,O-Z)      
C
      PARAMETER(NCmax=3, NAmax=5)
C
      INTEGER refSalt(10)
      REAL*8 mSalt(10),MOLAL(-NAmax:NCmax)
C
      INTEGER SOLID_R(1:30,4),SOLID_S(1:30,4),SOLID_N(1:30)
      COMMON/SLDS/SOLID_R,SOLID_S,SOLID_N
C
C   ..initialize molalities to zero
      DO 1 I=-NAmax,NCmax
        MOLAL(I)=0.D0
1     CONTINUE
C
C
      DO 2 I=1,NS
        nref=refSalt(I)
        DO 3 J=1,SOLID_N(nref)
          IF(SOLID_R(nref,J).NE.0 .AND. SOLID_R(nref,J).NE.99) THEN
            MOLAL(SOLID_R(nref,J))=MOLAL(SOLID_R(nref,J))
     >                            +mSalt(I)*SOLID_S(nref,J)
          ELSEIF(SOLID_R(nref,J) .EQ. 99) THEN
            WRITE(2,100) nref,J
            STOP 'INITIALIZE'
          ENDIF
3       CONTINUE
2     CONTINUE
C 
      RETURN
100   FORMAT(1X,'Solid = ',I2,2X,'Species = ',I2,'. Stop in SALTS')
      END
C
C +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
      SUBROUTINE CALCSAT(nSLD,iSolid,T,MFRAC,ACT,SRATIO)
      IMPLICIT REAL*8 (A-H,O-Z)
C
      PARAMETER(NCmax=3, NAmax=5)
C
      REAL*8 MULT,MFRAC(-NAmax:NCmax),ACT(-NAmax:NCmax),SRATIO(25)
C
      INTEGER iSolid(25)
      INTEGER SOLID_R(1:30,4),SOLID_S(1:30,4),SOLID_N(1:30)
      COMMON/SLDS/SOLID_R,SOLID_S,SOLID_N
C
      DO 1 I=1,nSLD
        nRef=iSOLID(I)
        IF(SOLID_N(nRef).EQ.0) THEN
C       ..no solid phase is assigned for this reference number:
          SRATIO(I)=-99.D0
        ELSE
          MULT=1.D0
          DO 2 J=1,SOLID_N(nRef)
            iRef=SOLID_R(nRef,J)
            xSpecies=MFRAC(iRef)
            IF(xSpecies.GT.0.D0) THEN
              MULT=MULT*(xSpecies*ACT(iRef))**SOLID_S(nRef,J)
            ELSE
C           ..species iRef is not present for this solution:
              SRATIO(I)=0.D0
              GOTO 1000
            ENDIF
2         CONTINUE
          SRATIO(I)=MULT/SATEQL(nRef,T)
1000      CONTINUE
        ENDIF
1     CONTINUE
C
      RETURN
      END
C
C +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C

      FUNCTION ASSOC(T)
      IMPLICIT REAL*8 (A-H,O-Z)
C
C     Dickson's equation follows below......
C
      DUM=562.694864456D0 - 102.5154D0*LOG(T) - 1.117033D-4*T**2
     >    + 0.2477538D0*T - 13273.75D0/T
      ASSOC=10.D0**(-DUM)
C
      END
C
C +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
      FUNCTION FC05(mHtotal,mSO4total,mH,MOLAL,T,CONST,xCAT,xAN,xNEUT,
     >              ACT)
      IMPLICIT REAL*8 (A-H,O-Z)
C
      PARAMETER(NCmax=3,NAmax=5,NNmax=1)
C
      INTEGER iACTC(NCmax),iACTA(NAmax)
      REAL*8 xCAT(NCmax),xAN(NAmax),xNEUT(NNmax),actC(NCmax),actA(NAmax)
      REAL*8 MOLAL(-NAmax:NCmax),ACT(-NAmax:NCmax)
      REAL*8 mHtotal,mSO4total,mH,KSTAR
C
C
C   ..adjust molalities of HSO4, SO4 and H for new free H+ estimate:
      MOLAL(-1)=mHtotal-mH
      MOLAL(-2)=mSO4total-MOLAL(-1)
      MOLAL(1)=mH

      SUM=0.D0
      DO 1 I=-NAmax,NCmax
        IF(I.NE.0) SUM=SUM + MOLAL(I)
1     CONTINUE
      var=55.508681D0+SUM
C
      DO 2 I=-NAmax,NCmax
        IF(I.LT.0) xAN(-I)=MOLAL(I)/var
        IF(I.EQ.0) xNEUT(1)=55.508681D0/var
        IF(I.GT.0) xCAT(I)=MOLAL(I)/var
2     CONTINUE
C
      CALL nCATION(xCAT,xAN,xNEUT,T,fH,1,iACTC,actC,.FALSE.)
      CALL nANION(xCAT,xAN,xNEUT,T,fHSO4,1,iACTA,actA,.FALSE.)
      CALL nANION(xCAT,xAN,xNEUT,T,fSO4,2,iACTA,actA,.FALSE.)
C
      KSTAR=CONST*fH*fSO4/fHSO4
      ACT(1)=fH
      ACT(-1)=fHSO4
      ACT(-2)=fSO4
C
C     now function F
C
      FC05=xCAT(1)*(xAN(1)+xAN(2))*KSTAR/(1.D0+xCAT(1)*KSTAR) - xAN(1)
C
      END
C
C +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
      SUBROUTINE ZBRENT(LOW,HIGH,mHtotal,mSO4total,MOLAL,T,EPS,xCAT,xAN,
     >                  xNEUT,ACT)
      IMPLICIT REAL*8(A-H,O-Z)
C
      PARAMETER(NCmax=3,NAmax=5,NNmax=1)
      PARAMETER(ITMAX=50,TOL=1.D-15)
C
      REAL*8 xCAT(NCmax),xAN(NAmax),xNEUT(NNmax)
      REAL*8 MOLAL(-NAmax:NCmax),ACT(-NAmax:NCmax)
      REAL*8 LOW,mHtotal,mSO4total
C
      A=LOW
      B=HIGH
      CONST=ASSOC(T) * 55.508681D0
C
      FA=FC05(mHtotal,mSO4total,A,MOLAL,T,CONST,xCAT,xAN,xNEUT,ACT)
      FB=FC05(mHtotal,mSO4total,B,MOLAL,T,CONST,xCAT,xAN,xNEUT,ACT)
C
      IF(FA*FB.GT.0.D0) THEN
        WRITE(2,'(1X,''MUST BRACKET ROOT! STOP IN ZBRENT'')')
        WRITE(2,'(1X,'' A,FA, B,FB:'',4(2X,E11.4))') A,FA,B,FB
        STOP 'ZBRENT'
      ENDIF
C
      C=B
      FC=FB
      DO 2 ITER=1,ITMAX
        IF(FB*FC .GT. 0.D0) THEN
          C=A
          FC=FA
          D=B-A
          E=D
        ENDIF
        IF(ABS(FC) .LT. ABS(FB)) THEN
          A=B
          B=C
          C=A
          FA=FB
          FB=FC
          FC=FA
        ENDIF
C
        TOL1=2.D0*EPS*ABS(B)+5.0D-1*TOL
        XM=5.0D-1*(C-B)
        IF((ABS(XM).LE.TOL1) .OR. (FB.EQ.0.D0)) THEN
C       ..current value of 'B' is best estimate. Return, accepting current
C         speciation.
          MOLAL(1)=B
          RETURN
        ENDIF
C
        IF((ABS(E).GE.TOL1) .AND. (ABS(FA).GT.ABS(FB))) THEN
          S=FB/FA
          IF(A.EQ.C) THEN
            P=2.D0*XM*S
            Q=1.D0-S
          ELSE
            Q=FA/FC
            R=FB/FC
            P=S*(2.D0*XM*Q*(Q-R)-(B-A)*(R-1.D0))
            Q=(Q-1.D0)*(R-1.D0)*(S-1.D0)
          ENDIF
          IF(P .GT. 0.D0) Q=-Q
          P=ABS(P)
          IF(2.D0*P .LT. MIN(3.D0*XM*Q-ABS(TOL1*Q),ABS(E*Q))) THEN
            E=D
            D=P/Q
          ELSE
            D=XM
            E=D
          ENDIF
        ELSE
          D=XM
          E=D
        ENDIF
        A=B
        FA=FB
        IF(ABS(D).GT.TOL1) THEN
          B=B+D
        ELSE
          B=B+SIGN(TOL1,XM)
        ENDIF
C
      FB=FC05(mHtotal,mSO4total,B,MOLAL,T,CONST,xCAT,xAN,xNEUT,ACT)
C
2     CONTINUE
C
      WRITE(*,*) ITMAX
      STOP 'ZBRENT MAX'
C
C
200   FORMAT(1X/' !ERROR: ITERATION LIMIT REACHED.'/
     >          ' ITMAX2= ',I2/
     >          ' *** TERMINATED IN ZBRENT ***')
C
C This routine, ZBRENT, is based upon the routine ZBRENT from the book 
C Numerical Recipes in FORTRAN (Cambridge University Press), Copyright
C (C) 1986, 1992 by Numerical Recipes Software. Used by permission. Use of 
C this routine other than as an integral part of the aerotr.for program
C requires an additional license from Numerical Recipes Software. Further
C distribution in any form is prohibited.
C
      END
C
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
      FUNCTION HNO3g(T,xH,xNO3,fH,fNO3)
      IMPLICIT REAL*8 (A-H,O-Z)
C
C     -------------------------------------
C     Calculate the partial pressure (atm).
C     -------------------------------------
C
      HNO3g=xH*fH*xNO3*fNO3/xKH_HNO3(T)    
C
      END
C
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
      FUNCTION HClg(T,xH,xCl,fH,fCl)
      IMPLICIT REAL*8 (A-H,O-Z)
C
C     -------------------------------------
C     Calculate the partial pressure (atm).
C     -------------------------------------
C
      HClg=xH*fH*xCl*fCl/xKH_HCl(T)                    
C
      END
C
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
      FUNCTION HBrg(T,xH,xBr,fH,fBr)
      IMPLICIT REAL*8 (A-H,O-Z)
C
C     -------------------------------------
C     Calculate the partial pressure (atm).
C     -------------------------------------
C
      HBrg=xH*fH*xBr*fBr/xKH_HBr(T)                    
C
      END
C
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
      REAL*8 FUNCTION NH3g(T,xH,xNH4,fH,fNH4)
      IMPLICIT REAL*8 (A-H,O-Z)
C
C     -------------------------------------
C     Calculate the partial pressure (atm).
C     -------------------------------------
C
      IF(xH .LE. 0.D0) THEN
        NH3g = 0.D0
      ELSE
        NH3g = xNH4*fNH4/(xH*fH)/xKH_NH3(T)
      ENDIF
C
      END
C
C +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
      FUNCTION H2SO4g(T,xH,xSO4,fH,fSO4)
      IMPLICIT REAL*8 (A-H,O-Z)
C
C     -------------------------------------
C     Calculate the partial pressure (atm).
C     -------------------------------------  
C
      H2SO4g =  (xH*fH)**2*xSO4*fSO4 /xKH_H2SO4(T)
C
      END
C
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
      FUNCTION xKH_HNO3(T)
      IMPLICIT REAL*8 (A-H,O-Z) 
C
      DUM = 385.972199D0 - 3020.3522D0/T - 71.001998D0*LOG(T)
     >     + 0.131442311D0*T - 0.420928363D-4*T**2
      xKH_HNO3=EXP(DUM)
C
      END
C
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
      FUNCTION xKH_HCl(T)
      IMPLICIT REAL*8 (A-H,O-Z)
C
      DUM = 6.4954D0 - 9.0027D3*(1.D0/298.15D0 - 1.D0/T) 
     >     - 65.346D0*(298.15D0/T - 1.D0 + LOG(T/298.15D0))
     >     + 0.078178D0*((298.15D0/T - 1.D0)*298.15D0 + T - 298.15D0)
      xKH_HCl=EXP(DUM)
C
      END
C
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
      FUNCTION xKH_HBr(T) 
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8 lnKH298
C
C     Mole fraction KH of HBr, from final fit to H2SO4-HBr-H2O
C     system, April 1998.
C     NB: deltaCp = delA + delB*T for the reaction.
C
      DATA lnKH298, delHr/12.5062D0, -85.15D3/ 
      DATA R/8.3144D0/
C
      delB = -1.74315D0
      delA = -170.94D0 - 298.15D0*delB

      DUM = lnKH298 + delHr/R*(1.D0/298.15D0 - 1.D0/T) 
     >       + (delA/R)*(298.15D0/T - 1.D0 + LOG(T/298.15D0))
     >       + (delB/2/R)*((298.15D0/T - 1.D0)*298.15D0 + T - 298.15D0)
      xKH_HBr = EXP(DUM)
C
      END
C
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
      FUNCTION xKH_NH3(T) 
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8 lnKa298
C
C***  In this updated routine I use an expression for the equilibrium
C     constant NH3(g) = NH4+(aq) - H+(aq) directly, rather than the
C     division KH/Ka which I used before. Future improvements may include
C     the variation of deltaCp[NH4+(aq)] with temperature (5/7/96).
C    
C   ..Equilibrium constant (units: atm-1) at 298.15 K is obtained from
C     the Clegg and Brimblecombe KH, and Clegg and Whitfield Ka:
         lnKa298 = LOG(60.72D0/5.6937D-10)
C   ..Delta H from Wagman [NH4+(aq)], and JANAF [NH3(g)]:
         deltaH = -132.15D3 + 45.898D3
C   ..Delta Cp from Roux et al. [NH4+(aq)], and JANAF [NH3(g)]:
         deltaCp = 70.D0 - 35.652D0
C
      Tr = 298.15D0
      R = 8.3144D0
C
C   ..mole fraction value at temperature T:
      xKH_NH3 = EXP(lnKa298 + (deltaH/R)*(1.D0/Tr-1/T)
     >              + (deltaCp/R)*(Tr/T-(1.D0+LOG(Tr/T))))
C
      END
C
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
      FUNCTION xKH_H2SO4(T) 
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8 lnKH298
C
      DATA lnKH298, delHr, delA, delB/
     >     27.4749403461D0, -180.383169D3, -1877.569397D0, 4.86846D0/ 
      DATA R/8.3144D0/
C
      T2=T
      T1=298.15D0
      dum=lnKH298 + delA/R*LOG(T2/T1) + delB/(2*R)*(T2-T1) 
     >     + (1.D0/R)*(-delHr + delA*T1 + delB/2*T1**2)
     >               *(1.D0/T2 - 1.D0/T1)
      xKH_H2SO4 = EXP(dum)
C
      END
C
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      FUNCTION H2Og(T,aw)
      IMPLICIT REAL*8 (A-H,O-Z)
C
C     Below 2oC an expression for pH2O0 based on Goff-Gratch is used
C     but modified so that the heat capacities of (supercooled) pure
C     water are those implied by an extrapolation of Hill's equation
C     of state. Above 2oC, an empirical fit to pH2O0 given in the CRC 
C     Handbook is used, valid to 100oC. At the 'join', the fit gives 
C     6.96289E-3 atm, and the lower temperature expression yields 
C     6.96157E-3 atm, a difference of 0.019%.
C
C
      REAL*8 LNK1,LNK2
      DATA LNK1, delHr/1.809429005D0, 45.077441D3/
      DATA R/8.3144D0/
      DATA Agas,Bgas,Cgas,Dgas,Egas/33.269811D0,0.00113261D0,
     >                              -1.09982D-5,
     >                              3.573575D-8,0.D0/
C     Hill-based liquid phase heat capacities:
      DATA Aliq,Bliq,Cliq,Dliq,Eliq/295.1612D0,-1.540498D0,
     >                              2.7023D-3,
     >                              0.0D0,0.0D0/
C
C
C
      IF(T .LE. 275.15D0) THEN
C     ..below 2 degrees C, used the Goff-Gratch based expression:
C
        delA=Agas - Aliq
        delB=Bgas - Bliq
        delC=Cgas - Cliq
        delD=Dgas - Dliq
        delE=Egas - Eliq
C
        T2=T
        T1=273.15D0
        DUM=delA/R*LOG(T2/T1) + delB/(2*R)*(T2-T1) 
     >     + delC/(6*R)*(T2**2-T1**2)
     >     + delD/(12*R)*(T2**3-T1**3)
     >     + delE/(2*R)*(1.D0/T2**2-1.D0/T1**2)
     >     + (1.D0/R)*(-delHr + delA*T1 + delB/2*T1**2 + delC/3*T1**3
     >               + delD/4*T1**4 - delE/T1)*(1.D0/T2 - 1.D0/T1)
C
        LNK2 = DUM + LNK1
        pH2O0 = EXP(LNK2) * 0.000986923D0
C                           ^- converts to atm.
C
      ELSE
C     ..else use empirical fit to CRC tabulated vapour pressures:
C
        DUM=23.54872D0 - 6459.987931D0/T - 0.022791752D0*T 
     >     + 1.6290826D-5*T**2
        pH2O0 = EXP(DUM)
      ENDIF
C
      H2Og = aw * pH2O0
C
      END
C
C +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
      SUBROUTINE nCATION(xCAT,xAN,xNEUT,TEMP,CATact,CATref,iACTC,actC,
     >                   PrintVar)
      IMPLICIT REAL*8(A-H,O-Z)
      LOGICAL PrintVar
C
      PARAMETER(NCmax=3,NAmax=5,NNmax=1, nUNIT=3)
C
      INTEGER PICK,CATref,iACTC(NCmax)
      REAL*8 Ix,MEc,MF,xCAT(NCmax),xAN(NAmax),xNEUT(NNmax),actC(NCmax)
C
      COMMON/SPEC1ES/NC,NA,NN,nCAT(NCmax),nAN(NAmax),nNEUT(NNmax)
      COMMON/COEFF1/RHO,A1ca(NCmax,NAmax),A2ca(NCmax,NAmax),
     >       B1ca(NCmax,NAmax),B2ca(NCmax,NAmax),
     >       Wnca(NNmax,NCmax,NAmax),Unca(NNmax,NCmax,NAmax),
     >       Vnca(NNmax,NCmax,NAmax),Wcca(NCmax,NCmax,NAmax),
     >       Waac(NAmax,NAmax,NCmax),Qncca(NNmax,NCmax,NCmax,NAmax),
     >       Qnaac(NNmax,NAmax,NAmax,NCmax),Ucca(NCmax,NCmax,NAmax),
     >       Uaac(NAmax,NAmax,NCmax),Wnn(NNmax,NNmax),Unn(NNmax,NNmax)
      COMMON/COEFF2/Ynnca(NNmax,NNmax,NCmax,NAmax),
     >              Xaaac(NAmax,NAmax,NAmax,NCmax),
     >              Xccca(NCmax,NCmax,NCmax,NAmax),
     >              Zccaa(NCmax,NCmax,NAmax,NAmax),
     >              Cnnn(NNmax,NNmax,NNmax)
      COMMON/CHARGE/CATchrg(NCmax),ANchrg(NAmax),
     >              Vca(NCmax,NAmax),Vac(NAmax,NCmax)
C
C
C:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-
C                                                             -
C Subroutine nCATION, includes multiple neutral species       -
C                                                             -
C Calculates cation activity coefficient, mole fraction basis -
C                                                             -
C:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-
C
C
C (1) Calculate Debye-Huckel constant Ax, 
C     and ionic strength Ix.
C
      Ax=AFT(TEMP)
      Ix=FUNCI(NC,NA,xCAT,xAN,CATchrg,ANchrg)
C
C
C (2) If CATref=0, then calculate the act coeffs of all cations
C     given in the indicator variable iACTC. If CATref NE 0, calculate 
C     the act coeff of this selected cation [nCAT(CATref)] only. 
C
      DO 1000 ICOUNT=1,NC
        IF(CATref.NE.0) THEN
          nSELECT=PICK(CATref,1)
        ELSEIF(iACTC(ICOUNT).EQ.1) THEN
          nSELECT=ICOUNT
          CATact=0.D0
        ELSEIF(iACTC(ICOUNT).EQ.0) THEN
C       ..activity coefficient not required, go to end of do loop.
          GOTO 1000
        ENDIF
C
C
C   ..Useful substitutions:
C
      SQRTIX=SQRT(Ix)
      zM=CATchrg(nSELECT)
C   ..and mole fraction functions 'F' and MF:
      FFdum=FF(xCAT,xAN,CATchrg,ANchrg)
      MF=FFdum*(1.D0 - (zM/2.D0)*FFdum)
C
C   #(01)# Calculate Debye-Huckel (long-range) terms.
C          ** Act coeff contribution is +SUMDH ** 
C
      SUMB1Ma=0.D0
      SUMB2Ma=0.D0
      SUMca1=0.D0
      SUMca2=0.D0
      DO 30 J=1,NA
        SUMB1Ma=SUMB1Ma+xAN(J)*B1ca(nSELECT,J)
     >                  *GFNC(A1ca(nSELECT,J)*SQRTIX)
        SUMB2Ma=SUMB2Ma+xAN(J)*B2ca(nSELECT,J)
     >                  *GFNC(A2ca(nSELECT,J)*SQRTIX)
C
        DO 31 I=1,NC
          SUMca1=SUMca1+xCAT(I)*xAN(J)*B1ca(I,J)*(zM**2
     >          *GFNC(A1ca(I,J)*SQRTIX)/(2.D0*Ix) + (1.D0-zM**2
     >          /(2.D0*Ix))*EXP(-A1ca(I,J)*SQRTIX))
          SUMca2=SUMca2+xCAT(I)*xAN(J)*B2ca(I,J)*(zM**2
     >          *GFNC(A2ca(I,J)*SQRTIX)/(2.D0*Ix) + (1.D0-zM**2
     >          /(2.D0*Ix))*EXP(-A2ca(I,J)*SQRTIX))
31      CONTINUE
30    CONTINUE
      SUMDH=-zM**2*Ax*((2.D0/RHO)*LOG(1.D0+RHO*SQRTIX)
     >      +SQRTIX*(1.D0-2.D0*Ix/zM**2)/(1.D0+RHO*SQRTIX))
     >      +SUMB1Ma+SUMB2Ma-SUMca1-SUMca2
C
C     
C   #(02)# Calculate unsymmetrical mixing terms. 
C          ** Act coeff contributions are +UNSYMc, -UNSYMcc and -UNSYMaa ** 
C
      UNSYMc=0.D0
      UNSYMcc=0.D0
      izM=INT(ZM)
      DO 32 I1=1,NC
        IF(I1.NE.nSELECT) THEN
          izI1=INT(CATchrg(I1))
          thetaMc=EFUNC(izM,izI1,Ax,Ix)
          thetadMc=EDFUNC(izM,izI1,Ax,Ix,thetaMc)
          UNSYMc=UNSYMc+2*xCAT(I1)*(thetaMc-xCAT(nSELECT)*
     >                       (thetaMc + thetadMc*(Ix-zM**2/2)))
          DO 33 I2=I1+1,NC
            IF(I2.NE.nSELECT .AND. I1.NE.NC) THEN         
              izI2=INT(CATchrg(I2))
              thetacc=EFUNC(izI1,izI2,Ax,Ix)
              thetadcc=EDFUNC(izI1,izI2,Ax,Ix,thetacc)
              UNSYMcc=UNSYMcc
     >             +2*xCAT(I1)*xCAT(I2)*(thetacc+thetadcc*(Ix-zM**2/2))
            ENDIF
33        CONTINUE
        ENDIF
32    CONTINUE  
C
      UNSYMaa=0.D0
      DO 34 J1=1,NA-1
       izJ1=INT(ANchrg(J1))
        DO 35 J2=J1+1,NA
          izJ2=INT(ANchrg(J2))
          thetaaa=EFUNC(izJ1,izJ2,Ax,Ix)
          thetadaa=EDFUNC(izJ1,izJ2,Ax,Ix,thetaaa)
            UNSYMaa=UNSYMaa
     >            +2*xAN(J1)*xAN(J2)*(thetaaa+thetadaa*(Ix-zM**2/2))
35      CONTINUE
34    CONTINUE  
C
C
C
C   #(3)# Calculate summation for Wcca parameters: 
C         ** Act coeff contribution is +SUMWCCA **
C
      SUMWcca=0.D0
      DO 1 J=1,NA
        SUMC=0.D0
        DO 2 I=1,NC
          IF(I.NE.nSELECT) THEN
            SUMC=SUMC+xCAT(I)*Wcca(nSELECT,I,J)
          ENDIF
2       CONTINUE
C
        SUMCC=0.D0
        DO 3 I1=1,NC-1
          DO 4 I2=I1+1,NC
            SUMCC=SUMCC+xCAT(I1)*xCAT(I2)*Wcca(I1,I2,J)
4         CONTINUE
3       CONTINUE
C
        SUMWcca=SUMWcca+2*Ea(J,xAN,ANchrg)*(SUMC-SUMCC)
1     CONTINUE
C
C
C   #(4)# Calculate summation for Waac parameters: 
C         ** Act coeff contribution is -SUMWAAC **
C
      SUMWaac=0.D0
      DO 5 I=1,NC
        SUMAA=0.D0
        DO 6 J1=1,NA-1
          DO 7 J2=J1+1,NA
            SUMAA=SUMAA+xAN(J1)*xAN(J2)*Waac(J1,J2,I)
7         CONTINUE
6       CONTINUE
        SUMWaac=SUMWaac+2*(Ec(I,xCAT,CATchrg)-
     >                   MEc(nSELECT,I,xCAT,CATchrg))*SUMAA
5     CONTINUE
C
C
C   #(5)# Calculate summation for Ucca parameters: 
C         ** Act coeff contribution is +SUMUCCA **
C
      SUMUcca=0.D0
      DO 8 J=1,NA
        SUMC=0.D0
        DO 9 I=1,NC
          IF(I.NE.nSELECT) THEN
            SUMC=SUMC+xCAT(I)*(2.D0*xCAT(nSELECT)
     >           /Vca(nSELECT,J)-xCAT(I)/Vca(I,J))*Ucca(nSELECT,I,J)
          ENDIF
9       CONTINUE
C
        SUMCC=0.D0
        DO 10 I1=1,NC-1
          DO 11 I2=I1+1,NC
            SUMCC=SUMCC+xCAT(I1)*xCAT(I2)*(xCAT(I1)/Vca(I1,J)
     >                  -xCAT(I2)/Vca(I2,J))*Ucca(I1,I2,J)
11        CONTINUE
10      CONTINUE
        SUMUcca=SUMUcca+2*Ea(J,xAN,ANchrg)*(SUMC-2.D0*SUMCC)
8     CONTINUE
C
C
C   #(6)# Calculate summation for Uaac parameters: 
C         ** Act coeff contribution is -SUMUAAC **
C
      SUMUaac=0.D0
      DO 12 I=1,NC
        SUMaa=0.D0
        DO 13 J1=1,NA-1
          DO 14 J2=J1+1,NA
            SUMaa=SUMaa+xAN(J1)*xAN(J2)
     >            *(xAN(J1)/Vac(J1,I)-xAN(J2)/Vac(J2,I))*Uaac(J1,J2,I)
14        CONTINUE
13      CONTINUE
        SUMUaac=SUMUaac+2*(2.D0*Ec(I,xCAT,CATchrg)
     >                     -MEc(nSELECT,I,xCAT,CATchrg))*SUMaa
12    CONTINUE
C
C
C   #(7)# Calculate summation for Wnca parameters: 
C         ** Act coeff contribution is +SUMWnca **
C
      SUMWnca=0.D0
      DO 36 K=1,NN       
        SUMca=0.D0
        DO 15 J=1,NA
          SUMC=0.D0
          DO 16 I=1,NC
            SUMC=SUMC+Ec(I,xCAT,CATchrg)*(CATchrg(I)+ANchrg(J))
     >                /(CATchrg(I)*ANchrg(J))*Wnca(K,I,J)
16        CONTINUE
          SUMca=SUMca+Ea(J,xAN,ANchrg)*((zM+ANchrg(J))/ANchrg(J)
     >                   *Wnca(K,nSELECT,J)-SUMC*(zM/2.D0+1.D0/FFdum))
15      CONTINUE
        SUMWnca=SUMWnca + xNEUT(K)*SUMca
36    CONTINUE
C
C
C   #(8)# Calculate summation for Unca parameters: 
C         ** Act coeff contribution is +SUMUnca **
C      
      SUMUnca=0.D0
      DO 37 K=1,NN
        SUMca=0.D0
        DO 17 J=1,NA
          SUMC=0.D0
          DO 18 I=1,NC
            SUMC=SUMC+2.D0*xCAT(I)*(CATchrg(I)+ANchrg(J))**2
     >           /(CATchrg(I)*ANchrg(J))*Unca(K,I,J)
18        CONTINUE
          DUM=(zM+ANchrg(J))**2/(zM*ANchrg(J))*Unca(K,nSELECT,J)
          SUMca=SUMca+xAN(J)*(DUM - SUMC)
17      CONTINUE
        SUMUnca=SUMUnca + xNEUT(K)*SUMca
37    CONTINUE
C
C
C   #(9)# Calculate summation for Vnca parameters: 
C         ** Act coeff contribution is +SUMVnca **
C
      SUMVnca=0.D0
      DO 38 K=1,NN
        SUMca=0.D0
        DO 19 J=1,NA
          SUMC=0.D0
          DO 20 I=1,NC
            SUMC=SUMC+3.D0*xCAT(I)*Vnca(K,I,J)
20        CONTINUE
          SUMca=SUMca+xAN(J)*(Vnca(K,nSELECT,J)-SUMC)
19      CONTINUE
        SUMVnca=SUMVnca + 4.D0*xNEUT(K)**2*SUMca
38    CONTINUE
C
C
C   #(10)# Calculate summation for Qncca parameters: 
C         ** Act coeff contribution is +SUMQncca **
C
      SUMQncca=0.D0
      DO 39 K=1,NN
        SUMcca=0.D0
        DO 21 J=1,NA
          SUMC=0.D0
          DO 22 I=1,NC
            IF(I.NE.nSELECT) THEN
              SUMC=SUMC+xCAT(I)*Qncca(K,nSELECT,I,J)
            ENDIF
22        CONTINUE
C
          SUMCC=0.D0
          DO 23 I1=1,NC-1
            DO 24 I2=I1+1,NC
              SUMCC=SUMCC+xCAT(I1)*xCAT(I2)*Qncca(K,I1,I2,J)
24          CONTINUE
23        CONTINUE
C
          SUMcca=SUMcca+Ea(J,xAN,ANchrg)*(SUMC-2*SUMCC)
21      CONTINUE
        SUMQncca=SUMQncca + 4.D0*xNEUT(K)*SUMcca      
39    CONTINUE
C
C
C   #(11)# Calculate summation for Qnaac parameters: 
C          ** Act coeff contribution is -SUMQnaac **
C
      SUMQnaac=0.D0
      DO 40 K=1,NN
        SUMaac=0.D0
        DO 25 I=1,NC
          SUMAA=0.D0
          DO 26 J1=1,NA-1
            DO 27 J2=J1+1,NA
              SUMAA=SUMAA+xAN(J1)*xAN(J2)*Qnaac(K,J1,J2,I)
27          CONTINUE
26        CONTINUE
          SUMaac=SUMaac+(2.D0*Ec(I,xCAT,CATchrg)-
     >                       MEc(nSELECT,I,xCAT,CATchrg))*SUMAA
25      CONTINUE
        SUMQnaac=SUMQnaac + 4.D0*xNEUT(K)*SUMaac      
40    CONTINUE
C
C
C   #(12)# Calculate summation for Ynnca parameters: 
C         ** Act coeff contribution is +SUMYnnca **
C
      SUMYnnca=0.D0
      DO 43 K1=1,NN-1
        DO 44 K2=K1+1,NN
          SUMa=0.D0
          DO 45 J=1,NA
            SUMc=0.D0
            DO 46 I=1,NC
              SUMc=SUMc+Ec(I,xCAT,CATchrg)*(CATchrg(I)+ANchrg(J))
     >                /(CATchrg(I)*ANchrg(J))*Ynnca(K1,K2,I,J)
46          CONTINUE
            SUMa=SUMa+Ea(J,xAN,ANchrg)*((zM+ANchrg(J))/ANchrg(J)
     >                *Ynnca(K1,K2,nSELECT,J) - (zM/2+2/FFdum)*SUMc)
45        CONTINUE
          SUMYnnca=SUMYnnca+xNEUT(K1)*xNEUT(K2)*SUMa
44      CONTINUE
43    CONTINUE
C
C
C
C   #(13)# Calculate summation for Xccca parameters: 
C         ** Act coeff contribution is +SUMXccca **
C
      SUMXccca=0.D0
      DO 47 J=1,NA
        SUMcc=0.D0
        DO 48 I1=1,NC-1
          DO 49 I2=I1+1,NC
            IF(I1.NE.nSELECT .AND. I2.NE.nSELECT) THEN
              SUMcc=SUMcc+xCAT(I1)*xCAT(I2)*Xccca(nSELECT,I1,I2,J)
            ENDIF
49        CONTINUE
48      CONTINUE
        SUMccc=0.D0
        DO 50 I1=1,NC-2
          DO 51 I2=I1+1,NC-1
            DO 52 I3=I2+1,NC
              SUMccc=SUMccc+xCAT(I1)*xCAT(I2)*xCAT(I3)*Xccca(I1,I2,I3,J)
52          CONTINUE
51        CONTINUE
50      CONTINUE
        SUMXccca=SUMXccca+Ea(J,xAN,ANchrg)*(SUMcc-2*SUMccc)
47    CONTINUE
C
C
C
C   #(14)# Calculate summation for Xaaac parameters: 
C         ** Act coeff contribution is -SUMXaaac **
C
      SUMXaaac=0.D0
      DO 53 I=1,NC
        SUMaaa=0.D0
        DO 54 J1=1,NA-2
          DO 55 J2=J1+1,NA-1
            DO 56 J3=J2+1,NA
              SUMaaa=SUMaaa+xAN(J1)*xAN(J2)*xAN(J3)*Xaaac(J1,J2,J3,I)
56          CONTINUE
55        CONTINUE
54      CONTINUE
        SUMXaaac=SUMXaaac+(2*Ec(I,xCAT,CATchrg)
     >                     -MEc(nSELECT,I,xCAT,CATchrg))*SUMaaa
53    CONTINUE
C
C
C
C   #(15)# Calculate summation for Zccaa parameters: 
C         ** Act coeff contribution is +SUMZccaa **
C
      SUMZccaa=0.D0
      DO 57 J1=1,NA-1
        DO 58 J2=J1+1,NA
          SUMc=0.D0
          DO 59 I=1,NC
            IF(I.NE.nSELECT) THEN
              SUMc=SUMc+xCAT(I)*Zccaa(nSELECT,I,J1,J2)
            ENDIF
59        CONTINUE
          SUMcc=0.D0
          DO 60 I1=1,NC-1
            DO 61 I2=I1+1,NC
              SUMcc=SUMcc+xCAT(I1)*xCAT(I2)*Zccaa(I1,I2,J1,J2)
61          CONTINUE
60        CONTINUE
          SUMZccaa=SUMZccaa+xAN(J1)*xAN(J2)*(FFdum*SUMc-(3*FFdum-MF)*
     >                                       SUMcc)
58      CONTINUE
57    CONTINUE
C
C
C   #(16)# Calculate summation for Wnn and Unn parameters: 
C          ** Act coeff contribution is -SUMnn **
      SUMnn=0.D0
      DO 41 K1=1,NN-1
        DO 42 K2=K1+1,NN
          SUMnn=SUMnn+xNEUT(K1)*xNEUT(K2)*(Wnn(K1,K2)+2*(xNEUT(K1)-
     >                                       xNEUT(K2))*Unn(K1,K2))
42      CONTINUE
41    CONTINUE
C
C
C
C   #(17)# Calculate summation for Cnnn parameters: 
C         ** Act coeff contribution is -SUMCnnn **
C
      SUMCnnn=0.D0
      DO 62 K1=1,NN-2
        DO 63 K2=K1+1,NN-1
          DO 64 K3=K2+1,NN
            SUMCnnn=SUMCnnn+2*xNEUT(K1)*xNEUT(K2)*xNEUT(K3)
     >                       *Cnnn(K1,K2,K3)
64        CONTINUE
63      CONTINUE
62    CONTINUE
C
C
C
C   #(18)# Calculate correction to infinite dilution reference state:
C          ** Act coeff contribution is -CORRact **
C          ** Expression below is formulated assuming ref. state of 
C             infinite dilution with respect to a pure (single) solvent
C             which *must* be neutral '1' **
C
      CORRact=0.D0
      DO 28 J=1,NA
        SUMC=0.D0
        DO 29 I=1,NC
          IF(I.NE.nSELECT) THEN
            SUMC=SUMC+Ec(I,xCAT,CATchrg)*(CATchrg(I)+ANchrg(J))/
     >                (CATchrg(I)*ANchrg(J))*Wnca(1,I,J)
          ENDIF
29      CONTINUE
        SUMC=SUMC*zM/2.D0
        CORRact=CORRact+
     >           Ea(J,xAN,ANchrg)*((1.D0-Ec(nSELECT,xCAT,CATchrg)/2.D0)
     >          *(zM+ANchrg(J))/ANchrg(J)*Wnca(1,nSELECT,J) - SUMC)
28    CONTINUE
C
C
C!!!! reference state correction set to zero for the tests !!!!!
C!!!! CORRact=0.D0
C!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
C
C
      IF(PrintVar) THEN
        WRITE(nUNIT,'(1X/,'' (CATION): Ax = '',E10.3)') Ax
        WRITE(nUNIT,'(1X,''Ix = '',E10.3)') Ix
        WRITE(nUNIT,100) SUMDH,UNSYMc,-UNSYMcc,-UNSYMaa,SUMWcca,
     >       -SUMWaac,SUMUcca,-SUMUaac,SUMWnca,SUMUnca,SUMVnca,
     >        SUMQncca,-SUMQnaac,SUMYnnca,SUMXccca,-SUMXaaac,
     >        SUMZccaa,-SUMnn,-SUMCnnn,-CORRact
      ENDIF
C
C
      DUMact=EXP(SUMDH + UNSYMc - UNSYMcc - UNSYMaa + SUMWCCA - SUMWAAC
     >                      + SUMUCCA - SUMUAAC + SUMWnca + SUMUnca
     >                      + SUMVnca + SUMQncca - SUMQnaac + SUMYnnca 
     >                      + SUMXccca - SUMXaaac + SUMZccaa - SUMnn
     >                      - SUMCnnn - CORRact)
C
C     #(14)# act coeffs of successive cations on run thru DO loop
      IF(CATref.EQ.0) THEN
        actC(nSELECT)=DUMact
      ELSE
C     #(15)# ACTC = activity coeff of selected (NREF) cation
        CATact=DUMact
C       ******
        RETURN
C       ******
      ENDIF
C     
1000  CONTINUE
C
      RETURN
100   FORMAT(1X,'SUMDH = ',E13.6,2X,'UNSYMc = ',E13.6,2X,'UNSYMcc = ',
     >       E13.6,2X,'UNSYMaa = ',E13.6,2X,'SUMWcca = ',E13.6,/2X,
     >      'SUMWaac = ',E13.6,2X,'SUMUcca = ',E13.6,2X,'SUMUaac = ',
     >       E13.6,2X,'SUMWnca = ',E13.6,2X,'SUMUnca = ',E13.6,/2X,
     >      'SUMVnca = ',E13.6,2X,'SUMQncca = ',E13.6,2X,'SUMQnaac = ',
     >       E13.6,2X,'SUMYnnca = ',E13.6,/2X,'SUMXccca = ',E13.6,2X,
     >      'SUMXaaac = ',E13.6,2X,'SUMZccaa = ',E13.6,2X,'SUMnn = ',
     >       E13.6,2X,'SUMCnnn = ',E13.6,/2X,'CORRact = ',E13.6)          
      END
C
C AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
C
      SUBROUTINE nANION(xCAT,xAN,xNEUT,TEMP,ANact,ANref,iACTA,actA,
     >                  PrintVar)
      IMPLICIT REAL*8(A-H,O-Z)
C
      LOGICAL PrintVar
      PARAMETER(NCmax=3,NAmax=5,NNmax=1, nUNIT=3)
C
      INTEGER PICK,ANref,iACTA(NAmax)
      REAL*8 Ix,xCAT(NCmax),xAN(NAmax),xNEUT(NNmax),actA(NAmax)
C
      COMMON/SPEC1ES/NC,NA,NN,nCAT(NCmax),nAN(NAmax),nNEUT(NNmax)
      COMMON/COEFF1/RHO,A1ca(NCmax,NAmax),A2ca(NCmax,NAmax),
     >       B1ca(NCmax,NAmax),B2ca(NCmax,NAmax),
     >       Wnca(NNmax,NCmax,NAmax),Unca(NNmax,NCmax,NAmax),
     >       Vnca(NNmax,NCmax,NAmax),Wcca(NCmax,NCmax,NAmax),
     >       Waac(NAmax,NAmax,NCmax),Qncca(NNmax,NCmax,NCmax,NAmax),
     >       Qnaac(NNmax,NAmax,NAmax,NCmax),Ucca(NCmax,NCmax,NAmax),
     >       Uaac(NAmax,NAmax,NCmax),Wnn(NNmax,NNmax),Unn(NNmax,NNmax)
      COMMON/COEFF2/Ynnca(NNmax,NNmax,NCmax,NAmax),
     >              Xaaac(NAmax,NAmax,NAmax,NCmax),
     >              Xccca(NCmax,NCmax,NCmax,NAmax),
     >              Zccaa(NCmax,NCmax,NAmax,NAmax),
     >              Cnnn(NNmax,NNmax,NNmax)
      COMMON/CHARGE/CATchrg(NCmax),ANchrg(NAmax),
     >              Vca(NCmax,NAmax),Vac(NAmax,NCmax)
C
C
C:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-
C                                                             -
C Subroutine nANION, includes multiple neutral species       -
C                                                             -
C Calculates cation activity coefficient, mole fraction basis -
C                                                             -
C:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-
C
C
C (1) Calculate Debye-Huckel constant Ax, 
C     and ionic strength Ix.
C
      Ax=AFT(TEMP)
      Ix=FUNCI(NC,NA,xCAT,xAN,CATchrg,ANchrg)

C
C
C (2) If ANref=0, then calculate the act coeffs of all anions
C     given in the indicator variable iACTA. If ANref NE 0, calculate 
C     the act coeff of this selected anion [nAN(ANref)] only. 
C
      DO 1000 ICOUNT=1,NA
        IF(ANref.NE.0) THEN
          nSELECT=PICK(ANref,-1)
        ELSEIF(iACTA(ICOUNT).EQ.1) THEN
          nSELECT=ICOUNT
          ANact=0.D0
        ELSEIF(iACTA(ICOUNT).EQ.0) THEN
C       ..activity coefficient not required, go to end of do loop.
          GOTO 1000
        ENDIF
C
C
C   ..Useful substitutions:
C
      SQRTIX=SQRT(Ix)
      zX=ANchrg(nSELECT)
C   ..and mole fraction functions 'F' and XF:
      FFdum=FF(xCAT,xAN,CATchrg,ANchrg)
      XF=FFdum*(1.D0 - (zX/2.D0)*FFdum)

C
C   #(01)# Calculate Debye-Huckel (long-range) terms.
C          ** Act coeff contribution is +SUMDH ** 
C
      SUMB1cX=0.D0
      SUMB2cX=0.D0
      SUMca1=0.D0
      SUMca2=0.D0
      DO 30 I=1,NC
        SUMB1cX=SUMB1cX+xCAT(I)*B1ca(I,nSELECT)
     >                  *GFNC(A1ca(I,nSELECT)*SQRTIX)
        SUMB2cX=SUMB2cX+xCAT(I)*B2ca(I,nSELECT)
     >                  *GFNC(A2ca(I,nSELECT)*SQRTIX)
C
        DO 31 J=1,NA
          SUMca1=SUMca1+xCAT(I)*xAN(J)*B1ca(I,J)*(zX**2
     >          *GFNC(A1ca(I,J)*SQRTIX)/(2.D0*Ix) + (1.D0-zX**2
     >          /(2.D0*Ix))*EXP(-A1ca(I,J)*SQRTIX))
          SUMca2=SUMca2+xCAT(I)*xAN(J)*B2ca(I,J)*(zX**2
     >          *GFNC(A2ca(I,J)*SQRTIX)/(2.D0*Ix) + (1.D0-zX**2
     >          /(2.D0*Ix))*EXP(-A2ca(I,J)*SQRTIX))
31      CONTINUE
30    CONTINUE
      SUMDH=-zX**2*Ax*((2.D0/13.D0)*LOG(1.D0+RHO*SQRTIX)
     >      +SQRTIX*(1.D0-2.D0*Ix/zX**2)/(1.D0+RHO*SQRTIX))
     >      +SUMB1cX+SUMB2cX-SUMca1-SUMca2
C
C     
C   #(02)# Calculate unsymmetrical mixing terms. 
C          ** Act coeff contributions are +SUMC, -SUMCC and -SUMAA ** 
C
      UNSYMa=0.D0
      UNSYMaa=0.D0
      izX=INT(zX)
      DO 32 J1=1,NA
        IF(J1.NE.nSELECT) THEN
          izJ1=INT(ANchrg(J1))
          thetacX=EFUNC(izX,izJ1,Ax,Ix)
          thetadcX=EDFUNC(izX,izJ1,Ax,Ix,thetacX)
          UNSYMa=UNSYMa+2*xAN(J1)*(thetacX-xAN(nSELECT)*
     >                       (thetacX + thetadcX*(Ix-zX**2/2)))
          DO 33 J2=J1+1,NA
            IF(J2.NE.nSELECT .AND. J1.NE.NA) THEN         
              izJ2=INT(ANchrg(J2))
              thetaaa=EFUNC(izJ1,izJ2,Ax,Ix)
              thetadaa=EDFUNC(izJ1,izJ2,Ax,Ix,thetaaa)
              UNSYMaa=UNSYMaa
     >             +2*xAN(J1)*xAN(J2)*(thetaaa+thetadaa*(Ix-zX**2/2))
            ENDIF
33        CONTINUE
        ENDIF
32    CONTINUE  
C
      UNSYMcc=0.D0
      DO 34 I1=1,NC-1
        izI1=INT(CATchrg(I1))
        DO 35 I2=I1+1,NC
          izI2=INT(CATchrg(I2))
          thetacc=EFUNC(izI1,izI2,Ax,Ix)
          thetadcc=EDFUNC(izI1,izI2,Ax,Ix,thetacc)
            UNSYMcc=UNSYMcc
     >            +2*xCAT(I1)*xCAT(I2)*(thetacc+thetadcc*(Ix-zX**2/2))
35      CONTINUE
34    CONTINUE  
C
C
C
C   #(3)# Calculate summation for Waac parameters: 
C         ** Act coeff contribution is +SUMWAAC **
      SUMWAAC=0.D0
      DO 1 I=1,NC
        SUMA=0.D0
        DO 2 J=1,NA
          IF(J.NE.nSELECT) THEN
            SUMA=SUMA+xAN(J)*Waac(nSELECT,J,I)
          ENDIF
2       CONTINUE
C
        SUMAA=0.D0
        DO 3 J1=1,NA-1
          DO 4 J2=J1+1,NA
            SUMAA=SUMAA+xAN(J1)*xAN(J2)*Waac(J1,J2,I)
4         CONTINUE
3       CONTINUE
C
        SUMWAAC=SUMWAAC+2*Ec(I,xCAT,CATchrg)*(SUMA-SUMAA)
1     CONTINUE
C
C
C   #(4)# Calculate summation for Wcca parameters: 
C         ** Act coeff contribution is -SUMWcca **
C
      SUMWCCA=0.D0
      DO 5 J=1,NA
        SUMCC=0.D0
        DO 6 I1=1,NC-1
          DO 7 I2=I1+1,NC
            SUMCC=SUMCC+xCAT(I1)*xCAT(I2)*Wcca(I1,I2,J)
7         CONTINUE
6       CONTINUE
        SUMWcca=SUMWcca+2*(Ea(J,xAN,ANchrg)-
     >                   XEa(nSELECT,J,xAN,ANchrg))*SUMCC
5     CONTINUE
C
C
C   #(5)# Calculate summation for Uaac parameters: 
C         ** Act coeff contribution is +SUMUaac **
C
      SUMUaac=0.D0
      DO 8 I=1,NC
        SUMA=0.D0
        DO 9 J=1,NA
          IF(J.NE.nSELECT) THEN
            SUMA=SUMA+xAN(J)*(2.D0*xAN(nSELECT)
     >           /Vac(nSELECT,I)-xAN(J)/Vac(J,I))*Uaac(nSELECT,J,I)
          ENDIF
9       CONTINUE
C
        SUMAA=0.D0
        DO 10 J1=1,NA-1
          DO 11 J2=J1+1,NA
            SUMAA=SUMAA+xAN(J1)*xAN(J2)*(xAN(J1)/Vac(J1,I)
     >                  -xAN(J2)/Vac(J2,I))*Uaac(J1,J2,I)
11        CONTINUE
10      CONTINUE
        SUMUaac=SUMUaac+2*Ec(I,xCAT,CATchrg)*(SUMA-2.D0*SUMAA)
8     CONTINUE
C
C
C   #(6)# Calculate summation for Ucca parameters: 
C         ** Act coeff contribution is -SUMUcca **
C
      SUMUcca=0.D0
      DO 12 J=1,NA
        SUMCC=0.D0
        DO 13 I1=1,NC-1
          DO 14 I2=I1+1,NC
            SUMCC=SUMCC+xCAT(I1)*xCAT(I2)
     >            *(xCAT(I1)/Vca(I1,J)-xCAT(I2)/Vca(I2,J))*Ucca(I1,I2,J)
14        CONTINUE
13      CONTINUE
        SUMUcca=SUMUcca+2*(2.D0*Ea(J,xAN,ANchrg)
     >                     -XEa(nSELECT,J,xAN,ANchrg))*SUMCC
12    CONTINUE
C
C
C   #(7)# Calculate summation for Wnca parameters: 
C         ** Act coeff contribution is +SUMWnca **
C
      SUMWnca=0.D0 
      DO 36 K=1,NN
        SUMca=0.D0      
        DO 15 I=1,NC
          SUMA=0.D0
          DO 16 J=1,NA
            SUMA=SUMA+Ea(J,xAN,ANchrg)*(CATchrg(I)+ANchrg(J))
     >                /(CATchrg(I)*ANchrg(J))*Wnca(K,I,J)
16        CONTINUE
          SUMca=SUMca+Ec(I,xCAT,CATchrg)*((CATchrg(I)+zX)/CATchrg(I)
     >               *Wnca(K,I,nSELECT) - SUMA*(zX/2.D0+1.D0/FFdum))
15      CONTINUE
        SUMWnca=SUMWnca + xNEUT(K)*SUMca
36    CONTINUE
C
C
C   #(8)# Calculate summation for Unca parameters: 
C         ** Act coeff contribution is +SUMUnca **
C      
      SUMUnca=0.D0
      DO 37 K=1,NN
        SUMca=0.D0
        DO 17 I=1,NC
          SUMA=0.D0
          DO 18 J=1,NA
            SUMA=SUMA+2.D0*xAN(J)*(CATchrg(I)+ANchrg(J))**2
     >                /(CATchrg(I)*ANchrg(J))*Unca(K,I,J)
18        CONTINUE
          DUM=(CATchrg(I)+zX)**2
     >        /(zX*CATchrg(I))*Unca(K,I,nSELECT)
          SUMca=SUMca+xCAT(I)*(DUM - SUMA)
17      CONTINUE
        SUMUnca=SUMUnca + xNEUT(K)*SUMca
37    CONTINUE
C
C
C   #(9)# Calculate summation for Vnca parameters: 
C         ** Act coeff contribution is +SUMVnca **
C
      SUMVnca=0.D0
      DO 38 K=1,NN
        SUMca=0.D0  
        DO 19 I=1,NC
          SUMA=0.D0
          DO 20 J=1,NA
            SUMA=SUMA+3.D0*xAN(J)*Vnca(K,I,J)
20        CONTINUE
          SUMca=SUMca+xCAT(I)*(Vnca(K,I,nSELECT)-SUMA)
19      CONTINUE
        SUMVnca=SUMVnca + 4.D0*xNEUT(K)**2*SUMca
38    CONTINUE
C
C
C   #(10)# Calculate summation for Qnaac parameters: 
C         ** Act coeff contribution is +SUMQnaac **
C
      SUMQnaac=0.D0
      DO 39 K=1,NN
        SUMaac=0.D0
        DO 21 I=1,NC
          SUMA=0.D0
          DO 22 J=1,NA
            IF(J.NE.nSELECT) THEN
              SUMA=SUMA+xAN(J)*Qnaac(K,nSELECT,J,I)
            ENDIF
22        CONTINUE
C  
          SUMAA=0.D0
          DO 23 J1=1,NA-1
            DO 24 J2=J1+1,NA
              SUMAA=SUMAA+xAN(J1)*xAN(J2)*Qnaac(K,J1,J2,I)
24          CONTINUE
23        CONTINUE
C
          SUMaac=SUMaac+Ec(I,xCAT,CATchrg)*(SUMA-2*SUMAA)
21      CONTINUE
        SUMQnaac=SUMQnaac + 4.D0*xNEUT(K)*SUMaac      
39    CONTINUE
C
C
C   #(11)# Calculate summation for Qncca parameters: 
C          ** Act coeff contribution is -SUMQncca **
C
      SUMQncca=0.D0
      DO 40 K=1,NN
        SUMcca=0.D0
        DO 25 J=1,NA
          SUMCC=0.D0
          DO 26 I1=1,NC-1
            DO 27 I2=I1+1,NC
              SUMCC=SUMCC+xCAT(I1)*xCAT(I2)*Qncca(K,I1,I2,J)
27          CONTINUE
26        CONTINUE
          SUMcca=SUMcca+(2.D0*Ea(J,xAN,ANchrg)-
     >                       XEa(nSELECT,J,xAN,ANchrg))*SUMCC
25      CONTINUE
        SUMQncca=SUMQncca + 4.D0*xNEUT(K)*SUMcca      
40    CONTINUE
C
C
C   #(12)# Calculate summation for Ynnca parameters: 
C         ** Act coeff contribution is +SUMYnnca **
C
      SUMYnnca=0.D0
      DO 43 K1=1,NN-1
        DO 44 K2=K1+1,NN
          SUMc=0.D0
          DO 45 I=1,NC
            SUMa=0.D0
            DO 46 J=1,NA
              SUMa=SUMa+Ea(J,xAN,ANchrg)*(CATchrg(I)+ANchrg(J))
     >                /(CATchrg(I)*ANchrg(J))*Ynnca(K1,K2,I,J)
46          CONTINUE
            SUMc=SUMc+Ec(I,xCAT,CATchrg)*((zX+CATchrg(I))/CATchrg(I)
     >                *Ynnca(K1,K2,I,nSELECT) - (zX/2+2/FFdum)*SUMa)
45        CONTINUE
          SUMYnnca=SUMYnnca+xNEUT(K1)*xNEUT(K2)*SUMc
44      CONTINUE
43    CONTINUE
C
C
C
C   #(13)# Calculate summation for Xaaac parameters: 
C         ** Act coeff contribution is +SUMXaaac **
C
      SUMXaaac=0.D0
      DO 47 I=1,NC
        SUMaa=0.D0
        DO 48 J1=1,NA-1
          DO 49 J2=J1+1,NA
            IF(J1.NE.nSELECT .AND. J2.NE.nSELECT) THEN
              SUMaa=SUMaa+xAN(J1)*xAN(J2)*Xaaac(nSELECT,J1,J2,I)
            ENDIF
49        CONTINUE
48      CONTINUE
        SUMaaa=0.D0
        DO 50 J1=1,NA-2
          DO 51 J2=J1+1,NA-1
            DO 52 J3=J2+1,NA
              SUMaaa=SUMaaa+xAN(J1)*xAN(J2)*xAN(J3)*Xaaac(J1,J2,J3,I)
52          CONTINUE
51        CONTINUE
50      CONTINUE
        SUMXaaac=SUMXaaac+Ec(I,xCAT,CATchrg)*(SUMaa-2*SUMaaa)
47    CONTINUE
C
C
C
C   #(14)# Calculate summation for Xccca parameters: 
C         ** Act coeff contribution is -SUMXccca **
C
      SUMXccca=0.D0
      DO 53 J=1,NA
        SUMccc=0.D0
        DO 54 I1=1,NC-2
          DO 55 I2=I1+1,NC-1
            DO 56 I3=I2+1,NC
              SUMccc=SUMccc+xCAT(I1)*xCAT(I2)*xCAT(I3)*Xccca(I1,I2,I3,J)
56          CONTINUE
55        CONTINUE
54      CONTINUE
        SUMXccca=SUMXccca+(2*Ea(J,xAN,ANchrg)
     >                     -XEa(nSELECT,J,xAN,ANchrg))*SUMccc
53    CONTINUE
C
C
C
C   #(15)# Calculate summation for Zccaa parameters: 
C         ** Act coeff contribution is +SUMZccaa **
C
      SUMZccaa=0.D0
      DO 57 I1=1,NC-1
        DO 58 I2=I1+1,NC
          SUMa=0.D0
          DO 59 J=1,NA
            IF(J.NE.nSELECT) THEN
              SUMa=SUMa+xAN(J)*Zccaa(I1,I2,nSELECT,J)
            ENDIF
59        CONTINUE
          SUMaa=0.D0
          DO 60 J1=1,NA-1
            DO 61 J2=J1+1,NA
              SUMaa=SUMaa+xAN(J1)*xAN(J2)*Zccaa(I1,I2,J1,J2)
61          CONTINUE
60        CONTINUE
          SUMZccaa=SUMZccaa+xCAT(I1)*xCAT(I2)*(FFdum*SUMa-(3*FFdum-XF)*
     >                                       SUMaa)
58      CONTINUE
57    CONTINUE
C
C
C
C   #(16)# Calculate summation for Wnn and Unn parameters: 
C          ** Act coeff contribution is -SUMnn **
      SUMnn=0.D0
      DO 41 K1=1,NN-1
        DO 42 K2=K1+1,NN
          SUMnn=SUMnn+xNEUT(K1)*xNEUT(K2)*(Wnn(K1,K2)+2*(xNEUT(K1)-
     >                                       xNEUT(K2))*Unn(K1,K2))
42      CONTINUE
41    CONTINUE
C
C
C
C
C   #(17)# Calculate summation for Cnnn parameters: 
C         ** Act coeff contribution is -SUMCnnn **
C
      SUMCnnn=0.D0
      DO 62 K1=1,NN-2
        DO 63 K2=K1+1,NN-1
          DO 64 K3=K2+1,NN
            SUMCnnn=SUMCnnn+2*xNEUT(K1)*xNEUT(K2)*xNEUT(K3)
     >                       *Cnnn(K1,K2,K3)
64        CONTINUE
63      CONTINUE
62    CONTINUE
C
C
C
C   #(12)# Calculate correction to infinite dilution reference state:
C          ** Act coeff contribution is -CORRact **
C          ** Expression below is formulated assuming ref. state of 
C             infinite dilution with respect to a pure (single) solvent
C             which *must* be neutral '1' **
C
      CORRact=0.D0
      DO 28 I=1,NC
        SUMA=0.D0
        DO 29 J=1,NA
          IF(J.NE.nSELECT) THEN
            SUMA=SUMA+Ea(J,xAN,ANchrg)*(CATchrg(I)+ANchrg(J))/
     >                (CATchrg(I)*ANchrg(J))*Wnca(1,I,J)
          ENDIF
29      CONTINUE
        SUMA=SUMA*zX/2.D0
        CORRact=CORRact+
     >           Ec(I,xCAT,CATchrg)*((1.D0-Ea(nSELECT,xAN,ANchrg)/2.D0)
     >          *(CATchrg(I)+zX)/CATchrg(I)*Wnca(1,I,nSELECT) - SUMA)
28    CONTINUE
C
C
C!!!! reference state correction set to zero for the tests !!!!!
C!!!! CORRact=0.D0
C!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
C
C
      IF(PrintVar) THEN
        WRITE(nUNIT,'(1X/,'' (ANION): Ax = '',E10.3)') Ax
        WRITE(nUNIT,'(1X,''Ix = '',E10.3)') Ix
        WRITE(nUNIT,100) SUMDH,UNSYMa,-UNSYMaa,-UNSYMcc,-SUMWcca,
     >        SUMWaac,-SUMUcca,SUMUaac,SUMWnca,SUMUnca,SUMVnca,
     >       -SUMQncca,SUMQnaac,SUMYnnca,SUMXaaac,-SUMXccca,SUMZccaa,
     >       -SUMnn,-SUMCnnn,-CORRact
      ENDIF
C
C
      DUMact=EXP(SUMDH+UNSYMa-UNSYMaa-UNSYMcc - SUMWCCA + SUMWAAC
     >                      - SUMUCCA + SUMUAAC + SUMWnca + SUMUnca
     >                      + SUMVnca - SUMQncca + SUMQnaac + SUMYnnca
     >                      + SUMXaaac - SUMXccca + SUMZccaa - SUMnn
     >                      - SUMCnnn - CORRact)
C
C
C     #(14)# act coeffs of successive anions on run thru DO loop
      IF(ANref.EQ.0) THEN
        actA(nSELECT)=DUMact
      ELSE
C     #(15)# ACTA = activity coeff of selected (NREF) anion
        ANact=DUMact
C       ******
        RETURN
C       ******
      ENDIF
C     
1000  CONTINUE
C
      RETURN
100   FORMAT(1X,'SUMDH = ',E13.6,2X,'UNSYMa = ',E13.6,2X,'UNSYMaa = ',
     >       E13.6,2X,'UNSYMcc = ',E13.6,2X,'SUMWcca = ',E13.6,/2X,
     >      'SUMWaac = ',E13.6,2X,'SUMUcca = ',E13.6,2X,'SUMUaac = ',
     >       E13.6,2X,'SUMWnca = ',E13.6,2X,'SUMUnca = ',E13.6,/2X,
     >      'SUMVnca = ',E13.6,2X,'SUMQncca = ',E13.6,2X,'SUMQnaac = ',
     >       E13.6,2X,'SUMYnnca = ',E13.6,/2X,'SUMXaaac = ',E13.6,2X,
     >      'SUMXccca = ',E13.6,2X,'SUMZccaa = ',E13.6,2X,'SUMnn = ',
     >       E13.6,2X,'SUMCnnn = ',E13.6,/2X,'CORRact = ',E13.6)          
      END
C
C +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
      SUBROUTINE nNEUTRAL(xCAT,xAN,xNEUT,TEMP,NEUTact,NEUTref,iACTN,
     >                    actN,AW,G,PrintVar)
      IMPLICIT REAL*8(A-H,O-Z)
C
      LOGICAL PrintVar
      PARAMETER(NCmax=3,NAmax=5,NNmax=1, nUNIT=3)
C
      INTEGER PICK,NEUTref,iACTN(NNmax)
      REAL*8 Ix,xCAT(NCmax),xAN(NAmax),xNEUT(NNmax),actN(NNmax),
     >       NEUTact
C
      COMMON/SPEC1ES/NC,NA,NN,nCAT(NCmax),nAN(NAmax),nNEUT(NNmax)
      COMMON/COEFF1/RHO,A1ca(NCmax,NAmax),A2ca(NCmax,NAmax),
     >       B1ca(NCmax,NAmax),B2ca(NCmax,NAmax),
     >       Wnca(NNmax,NCmax,NAmax),Unca(NNmax,NCmax,NAmax),
     >       Vnca(NNmax,NCmax,NAmax),Wcca(NCmax,NCmax,NAmax),
     >       Waac(NAmax,NAmax,NCmax),Qncca(NNmax,NCmax,NCmax,NAmax),
     >       Qnaac(NNmax,NAmax,NAmax,NCmax),Ucca(NCmax,NCmax,NAmax),
     >       Uaac(NAmax,NAmax,NCmax),Wnn(NNmax,NNmax),Unn(NNmax,NNmax)
      COMMON/COEFF2/Ynnca(NNmax,NNmax,NCmax,NAmax),
     >              Xaaac(NAmax,NAmax,NAmax,NCmax),
     >              Xccca(NCmax,NCmax,NCmax,NAmax),
     >              Zccaa(NCmax,NCmax,NAmax,NAmax),
     >              Cnnn(NNmax,NNmax,NNmax)
      COMMON/CHARGE/CATchrg(NCmax),ANchrg(NAmax),
     >              Vca(NCmax,NAmax),Vac(NAmax,NCmax)
C
C
C:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-
C                                                             -
C Subroutine NEUTRAL                                          -
C                                                             -
C Calculates activities of neutral species including water    -
C                                                             -
C NB: for all calculations the first neutral species must be  -
C     the solvent.                                            -
C                                                             -
C:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-:-
C
      DO 1000 ICOUNT=1,NN
        IF(NEUTref.NE.0) THEN
          nSELECT=PICK(NEUTref,0)
        ELSEIF(iACTN(ICOUNT).EQ.1) THEN
          nSELECT=ICOUNT
          NEUTact=0.D0
        ELSEIF(iACTN(ICOUNT).EQ.0) THEN
C       ..activity coefficient not required, go to end of do loop.
          GOTO 1000
        ENDIF
C
C
C NB: it is possible that for neutral solutes the extended DH function
C     below will not be needed. If so, it can be put into an IF block.
C
C (1) Calculate Debye-Huckel constant Ax, 
C     and ionic strength Ix.
C
      Ax=AFT(TEMP)
      Ix=FUNCI(NC,NA,xCAT,xAN,CATchrg,ANchrg)
C
C   ..Useful substitutions:
C
      SQRTIX=SQRT(Ix)
C
C   #(01)# Calculate Debye-Huckel (long-range) terms.
C          ** Act coeff contribution is +SUMDH ** 
C
      SUMca1=0.D0
      SUMca2=0.D0
      DO 1 I=1,NC
        DO 2 J=1,NA
          SUMca1=SUMca1+xCAT(I)*xAN(J)*B1ca(I,J)*EXP(-A1ca(I,J)*SQRTIX)
          SUMca2=SUMca2+xCAT(I)*xAN(J)*B2ca(I,J)*EXP(-A2ca(I,J)*SQRTIX)
2       CONTINUE
1     CONTINUE
      SUMDH=2*Ax*Ix**1.5D0/(1.D0+RHO*SQRTIX) - SUMca1 - SUMca2
C
C
C   #(02)# Calculate unsymmetrical mixing terms. 
C          ** Act coeff contributions are -UNSYMaa and -UNSYMcc ** 
C
      UNSYMaa=0.D0
      DO 3 J1=1,NA-1
        izJ1=INT(ANchrg(J1))
        DO 4 J2=J1+1,NA
          izJ2=INT(ANchrg(J2))
          thetaaa=EFUNC(izJ1,izJ2,Ax,Ix)
          thetadaa=EDFUNC(izJ1,izJ2,Ax,Ix,thetaaa)
          UNSYMaa=UNSYMaa
     >            +2*xAN(J1)*xAN(J2)*(thetaaa+thetadaa*Ix)
4       CONTINUE
3     CONTINUE  
C
      UNSYMcc=0.D0
      DO 5 I1=1,NC-1
        izI1=INT(CATchrg(I1))
        DO 6 I2=I1+1,NC
          izI2=INT(CATchrg(I2))
          thetacc=EFUNC(izI1,izI2,Ax,Ix)
          thetadcc=EDFUNC(izI1,izI2,Ax,Ix,thetacc)
          UNSYMcc=UNSYMcc
     >            +2*xCAT(I1)*xCAT(I2)*(thetacc+thetadcc*Ix)
6       CONTINUE
5     CONTINUE  
C
C
C
C   ..mole fraction function 'F' for later use:
      FFdum=FF(xCAT,xAN,CATchrg,ANchrg)
C
C
C   #(3)# Calculate summation for Waac parameters: 
C         ** Act coeff contribution is -SUMWaac **
      SUMWaac=0.D0
      DO 7 I=1,NC
        SUMaa=0.D0
        DO 8 J1=1,NA-1
          DO 9 J2=J1+1,NA
            SUMaa=SUMaa+xAN(J1)*xAN(J2)*Waac(J1,J2,I)
9         CONTINUE
8       CONTINUE
C
        SUMWaac=SUMWaac+2*EC(I,xCAT,CATchrg)*SUMaa
7     CONTINUE
C
C
C   #(4)# Calculate summation for Wcca parameters: 
C         ** Act coeff contribution is -SUMWcca **
C
      SUMWcca=0.D0
      DO 10 J=1,NA
        SUMcc=0.D0
        DO 11 I1=1,NC-1
          DO 12 I2=I1+1,NC
            SUMcc=SUMcc+xCAT(I1)*xCAT(I2)*Wcca(I1,I2,J)
12        CONTINUE
11      CONTINUE
        SUMWcca=SUMWcca+2*EA(J,xAN,ANchrg)*SUMcc
10    CONTINUE
C
C
C   #(5)# Calculate summation for Uaac parameters: 
C         ** Act coeff contribution is -SUMUaac **
C
      SUMUaac=0.D0
      DO 13 I=1,NC
        SUMaa=0.D0
        DO 14 J1=1,NA-1
          DO 15 J2=J1+1,NA
            SUMaa=SUMaa+xAN(J1)*xAN(J2)
     >            *(xAN(J1)/Vac(J1,I)-xAN(J2)/Vac(J2,I))*Uaac(J1,J2,I)
15        CONTINUE
14      CONTINUE
C
        SUMUaac=SUMUaac+4*EC(I,xCAT,CATchrg)*SUMaa
13    CONTINUE
C
C
C   #(6)# Calculate summation for Ucca parameters: 
C         ** Act coeff contribution is -SUMUcca **
C
      SUMUcca=0.D0
      DO 16 J=1,NA
        SUMcc=0.D0
        DO 17 I1=1,NC-1
          DO 18 I2=I1+1,NC
            SUMcc=SUMcc+xCAT(I1)*xCAT(I2)
     >            *(xCAT(I1)/Vca(I1,J)-xCAT(I2)/Vca(I2,J))*Ucca(I1,I2,J)
18        CONTINUE
17      CONTINUE
        SUMUcca=SUMUcca+4*EA(J,xAN,ANchrg)*SUMcc
16    CONTINUE
C
C
C   #(7)# Calculate summation for Wnca parameters: 
C         ** Act coeff contribution is +SUMWnca **
C
      SUMWnca=0.D0       
      DO 19 I=1,NC
        DO 20 J=1,NA
          SUMn=0.D0
          DO 31 K=1,NN
            SUMn=SUMn+xNEUT(K)*Wnca(K,I,J)
31        CONTINUE
          DUM=(CATchrg(I)+ANchrg(J))/(CATchrg(I)*ANchrg(J))
     >        *(Wnca(nSELECT,I,J)-SUMn)
          SUMWnca=SUMWnca+(1.D0/FFdum)*EC(I,xCAT,CATchrg)
     >                                *EA(J,xAN,ANchrg)*DUM
20      CONTINUE
19    CONTINUE
C
C
C   #(8)# Calculate summation for Unca parameters: 
C         ** Act coeff contribution is +SUMUnca **
C      
      SUMUnca=0.D0
      DO 21 I=1,NC
        DO 22 J=1,NA
          SUMn=0.D0
          DO 32 K=1,NN
            SUMn=SUMn+xNEUT(K)*Unca(K,I,J)
32        CONTINUE
          DUM=(CATchrg(I)+ANchrg(J))**2/(CATchrg(I)*ANchrg(J))
     >        *(Unca(nSELECT,I,J)-2.D0*SUMn)
          SUMUnca=SUMUnca+xCAT(I)*xAN(J)*DUM
22      CONTINUE
21    CONTINUE
C
C
C   #(9)# Calculate summation for Vnca parameters: 
C         ** Act coeff contribution is +SUMVnca **
C
      SUMVnca=0.D0
      DO 23 I=1,NC
        DO 24 J=1,NA
          SUMn=0.D0
          DO 33 K=1,NN
            SUMn=SUMn+xNEUT(K)**2*Vnca(K,I,J)
33        CONTINUE
          SUMVnca=SUMVnca+4*xCAT(I)*xAN(J)*(2*xNEUT(nSELECT)
     >                     *Vnca(nSELECT,I,J)-3*SUMn)
24      CONTINUE
23    CONTINUE
C
C
C   #(10)# Calculate summation for Qnaac parameters: 
C         ** Act coeff contribution is +SUMQnaac **
C
      SUMQnaac=0.D0
      DO 25 I=1,NC
        SUMaa=0.D0
        DO 26 J1=1,NA-1
          DO 27 J2=J1+1,NA
            SUMn=0.D0
            DO 34 K=1,NN
              SUMn=SUMn+xNEUT(K)*Qnaac(K,J1,J2,I)
34          CONTINUE
            SUMaa=SUMaa+xAN(J1)*xAN(J2)*(Qnaac(nSELECT,J1,J2,I)-2*SUMn)
27        CONTINUE
26      CONTINUE
C
        SUMQnaac=SUMQnaac+4*EC(I,xCAT,CATchrg)*SUMaa
25    CONTINUE
C
C
C   #(11)# Calculate summation for Qncca parameters: 
C          ** Act coeff contribution is +SUMQncca **
C
      SUMQncca=0.D0
      DO 28 J=1,NA
        SUMcc=0.D0
        DO 29 I1=1,NC-1
          DO 30 I2=I1+1,NC
            SUMn=0.D0
            DO 35 K=1,NN
              SUMn=SUMn+xNEUT(K)*Qncca(K,I1,I2,J)
35          CONTINUE
            SUMcc=SUMcc+xCAT(I1)*xCAT(I2)
     >                  *(Qncca(nSELECT,I1,I2,J)-2*SUMn)
30        CONTINUE
29      CONTINUE
C
        SUMQncca=SUMQncca+4*EA(J,xAN,ANchrg)*SUMcc
28    CONTINUE
C
C
C
C   #(12)# Calculate summation for Ynnca parameters: 
C          ** Act coeff contribution is +SUMYnnca **
C
      SUMYnnca=0.D0
      DO 44 I=1,NC
        DO 45 J=1,NA
          SUMn=0.D0
          DO 46 K=1,NN
            IF(K.NE.nSELECT) THEN
              SUMn=SUMn+xNEUT(K)*Ynnca(nSELECT,K,I,J)
            ENDIF
46        CONTINUE
C
          SUMnn=0.D0
          DO 47 K1=1,NN-1
            DO 48 K2=K1+1,NN
              SUMnn=SUMnn+xNEUT(K1)*xNEUT(K2)*Ynnca(K1,K2,I,J)
48          CONTINUE
47        CONTINUE
C
          SUMYnnca=SUMYnnca+(1.D0/FFdum)*Ec(I,xCAT,CATchrg)*
     >             Ea(J,xAN,ANchrg)*(CATchrg(I)+ANchrg(J))/(CATchrg(I)*
     >             ANchrg(J))*(SUMn-2*SUMnn)
45      CONTINUE
44    CONTINUE
C
C
C
C   #(13)# Calculate summation for Xccca parameters: 
C          ** Act coeff contribution is -SUMXccca **
C
      SUMXccca=0.D0
      DO 49 J=1,NA
        SUMccc=0.D0
        DO 50 I1=1,NC-2
          DO 51 I2=I1+1,NC-1
            DO 52 I3=I2+1,NC
              SUMccc=SUMccc+xCAT(I1)*xCAT(I2)*xCAT(I3)*Xccca(I1,I2,I3,J)
52          CONTINUE
51        CONTINUE
50      CONTINUE
C
        SUMXccca=SUMXccca+2*Ea(J,xAN,ANchrg)*SUMccc
49    CONTINUE
C
C
C
C   #(14)# Calculate summation for Xaaac parameters: 
C          ** Act coeff contribution is -SUMXaaac **
C
      SUMXaaac=0.D0
      DO 53 I=1,NC
        SUMaaa=0.D0
        DO 54 J1=1,NA-2
          DO 55 J2=J1+1,NA-1
            DO 56 J3=J2+1,NA
              SUMaaa=SUMaaa+xAN(J1)*xAN(J2)*xAN(J3)*Xaaac(J1,J2,J3,I)
56          CONTINUE
55        CONTINUE
54      CONTINUE
C
        SUMXaaac=SUMXaaac+2*Ec(I,xCAT,CATchrg)*SUMaaa
53    CONTINUE
C
C
C
C   #(15)# Calculate summation for Zccaa parameters: 
C          ** Act coeff contribution is -SUMZccaa **
C
      SUMZccaa=0.D0
      DO 57 I1=1,NC-1
        DO 58 I2=I1+1,NC
          DO 59 J1=1,NA-1
            DO 60 J2=J1+1,NA
              SUMZccaa=SUMZccaa+2*FFdum*xCAT(I1)*xCAT(I2)*xAN(J1)*
     >                          xAN(J2)*Zccaa(I1,I2,J1,J2)
60          CONTINUE
59        CONTINUE
58      CONTINUE
57    CONTINUE
C
C
C   #(16)# Calculate first summation for Wnn and Unn parameters: 
C          ** Act coeff contribution is +SUMnn1 **
C
      SUMnn1=0.D0
      DO 41 K1=1,NN
        IF(K1.NE.nSELECT) THEN
          DUM=1.D0-xNEUT(nSELECT)
          SUMnn1=SUMnn1+xNEUT(K1)*(DUM*Wnn(nSELECT,K1) + 
     >    (2*(xNEUT(nSELECT)-xNEUT(K1))*DUM+xNEUT(K1))*Unn(nSELECT,K1))
        ENDIF
41    CONTINUE
C
C
C   #(17)# Calculate second summation for Wnn, Unn and Cnnn parameters: 
C          ** Act coeff contribution is -SUMnn2 **
C
      SUMnn2=0.D0
      DO 42 K1=1,NN-1
        DO 43 K2=K1+1,NN
          IF(K1.NE.nSELECT .AND. K2.NE.nSELECT) THEN
            SUMnn2=SUMnn2+xNEUT(K1)*xNEUT(K2)*(Wnn(K1,K2) + 
     >      2*(xNEUT(K1)-xNEUT(K2))*Unn(K1,K2)-(1.D0-2*xNEUT(nSELECT))*
     >         Cnnn(nSELECT,K1,K2))
          ENDIF
43      CONTINUE
42    CONTINUE
C
C
C   #(18)# Calculate summation for Cnnn parameters: 
C          ** Act coeff contribution is -SUMCnnn **
C
      SUMCnnn=0.D0
      DO 61 K1=1,NN-2
        DO 62 K2=K1+1,NN-1
          DO 63 K3=K2+1,NN
            IF(K1.NE.nSELECT .AND. K2.NE.nSELECT .AND. K3.NE.nSELECT)
     >        SUMCnnn=SUMCnnn+2*xNEUT(K1)*xNEUT(K2)*xNEUT(K3)*
     >                        Cnnn(K1,K2,K3)
63        CONTINUE
62      CONTINUE
61    CONTINUE
C
C
C
C   #(19)# Reference state correction for neutral solutes:
C       ** Act coeff contribution is -CORRact **
C       ** Expression below is formulated assuming ref. state of 
C          infinite dilution with respect to a pure (single) solvent
C          which *must* be neutral '1' **

      IF(nSELECT.NE.1) THEN     
        CORRact=Wnn(nSELECT,1)-Unn(nSELECT,1)
      ELSE
        CORRact=0.D0
      ENDIF
C
C!!!! reference state correction set to zero for the tests !!!!!
C!!!! CORRact=0.D0
C!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
C
C
      IF(PrintVar) THEN
        WRITE(nUNIT,'(1X/,'' (NEUTRAL): Ax = '',E14.7)') Ax
        WRITE(nUNIT,'(1X,''Ix = '',E10.3)') Ix
        WRITE(nUNIT,100) SUMDH,-UNSYMaa,-UNSYMcc,-SUMWcca,-SUMWaac,
     >       -SUMUcca,-SUMUaac,SUMWnca,SUMUnca,SUMVnca,SUMQncca,
     >        SUMQnaac,SUMYnnca,-SUMXccca,-SUMXaaac,-SUMZccaa,SUMnn1,
     >       -SUMnn2,-SUMCnnn,-CORRact
      ENDIF
C
C
C
      DUMact=EXP(SUMDH-UNSYMaa-UNSYMcc - SUMWcca - SUMWaac
     >                      - SUMUcca - SUMUaac + SUMWnca + SUMUnca
     >                      + SUMVnca + SUMQncca + SUMQnaac + SUMYnnca
     >                      - SUMXccca - SUMXaaac - SUMZccaa  + SUMnn1
     >                      - SUMnn2 - SUMCnnn - CORRact)
C
C
C   #(15)# act coeffs of successive neutrals run thru DO loop
      IF(NEUTref.EQ.0) THEN
        actN(nSELECT)=DUMact
        IF(nSELECT.EQ.1) THEN
          aw=xNEUT(1)*DUMact
          g=LOG(aw)/LOG(xNEUT(1))
        ENDIF
      ELSE
C   #(16)# NEUTact = activity coeff of selected (NREF) neutral
        NEUTact=DUMact
        IF(nSELECT.EQ.1) THEN
          aw=xNEUT(1)*DUMact
          g=LOG(aw)/LOG(xNEUT(1))
        ENDIF
C       ******
        RETURN
C       ******
      ENDIF
C     
1000  CONTINUE
C
      RETURN
100   FORMAT(1X,'SUMDH = ',E13.6,2X,'UNSYMaa = ',E13.6,2X,'UNSYMcc = ',
     >       E13.6,2X,'SUMWcca = ',E13.6,2X,'SUMWaac = ',E13.6,/2X,
     >      'SUMUcca = ',E13.6,2X,'SUMUaac = ',E13.6,2X,'SUMWnca = ',
     >       E13.6,2X,'SUMUnca = ',E13.6,2X,'SUMVnca = ',E13.6,/2X,
     >      'SUMQncca = ',E13.6,2X,'SUMQnaac = ',E13.6,2X,'SUMYnnca = ',
     >      E13.6,2X,'SUMXccca = ',E13.6,2X,'SUMXaaac = ',E13.6,/2X,
     >      'SUMZccaa = ',E13.6,2X,'SUMnn1 = ',E13.6,2X,'SUMnn2 = ',
     >      E13.6,2X,'SUMCnnn = ',E13.6,2X,'CORRact = ',E13.6)
C
      END
C
C +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
      FUNCTION AFT(T)
      IMPLICIT REAL*8 (A-H,O-Z)
C
C     function yields molality based DH coefficient for absolute
C     temperature T (K). Values below 273.15 K are an extrapolation.
C
      PARAMETER(R=8.3144D0, Tr = 273.15D0,
     >                      AphiTr = 0.376421485D0,
     >                      AlRTr = 0.600305325D0,
     >                      AjR = 1.677818299D0,
     >                      dAjRdT = 0.1753875763D0,
     >                      N = 17,
     >                      XMIN = 234.15D0, XMAX = 373.15D0)
C
C--------------------------------------------------------------------
C               Tr must be >=234.15K, 
C               Aphi = DH constant at Tr,
C               AlRTr = AL / (RTr) at Tr,
C               AjR = AJ / R at Tr,
C               dAjRdT = d (AJ/R) / dT at Tr.
C--------------------------------------------------------------------
C
      REAL*8 AA(0:18),TX(0:18)
C
      DATA AA/0.797256081240D+00,0.573389669896D-01,
     > 0.977632177788D-03, 0.489973732417D-02,-0.313151784342D-02,
     > 0.179145971002D-02,-0.920584241844D-03, 0.443862726879D-03,
     >-0.203661129991D-03, 0.900924147948D-04,-0.388189392385D-04,
     > 0.164245088592D-04,-0.686031972567D-05, 0.283455806377D-05,
     >-0.115641433004D-05, 0.461489672579D-06,-0.177069754948D-06,
     > 0.612464488231D-07,-0.175689013085D-07/
C
C   ..for T>Tr, polynomial reproduces Archer's values:
      If(T .GE. Tr) then
        X=(2.D0*T-XMAX-XMIN)/(XMAX-XMIN)
        TX(0)=1.D0
        TX(1)=X
        AFT=0.5*AA(0)*TX(0) + AA(1)*TX(1)
        DO 1 I=1,N
          TX(I+1)=2*X*TX(I) - TX(I-1)
          AFT=AFT + AA(I+1)*TX(I+1)
1       CONTINUE
C
C   ..for T<Tr, we have an empirical extrapolation:
      Else
        AlTr = AlRTr * R * Tr
        AjTr = AjR * R
        dAJdT = dAjRdT * R 
        a=1.45824246467D0
        x=(dAJdT - a)/(2.D0*Tr)
C
        AFT=AphiTr + AlTr/(4*R)*(1.D0/Tr-1.D0/T) 
     >      + AjTr/(4.D0*R)*(LOG(T/Tr)+Tr/T-1.D0)
     >      + a/(8.D0*R)*(T-Tr**2/T-2.D0*Tr*LOG(T/Tr))
     >      + x/(4.D0*R)*(T**2/6.D0+Tr**2/2.D0-Tr**2*LOG(T/Tr)
     >                    -2.D0/3.D0*Tr**3/T)
C
      Endif
C
C   ..convert to mole fraction basis:
      AFT=AFT * 7.4504148D0
C
      END
C
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
      FUNCTION GFNC(X)
      IMPLICIT REAL*8 (A-H,O-Z)
      IF(X.EQ.0.D0) THEN
        GFNC=0.D0
      ELSE
        GFNC=2.D0*(1.D0-(1.D0+X)*EXP(-X))/X**2
      ENDIF
      END
C
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
      FUNCTION EFUNC(ICHARG,JCHARG,A,XION)
       IMPLICIT REAL*8(A-H,O-Z)
       REAL*8 XIJ(3),J0(3)
C
       IF((ICHARG.EQ.JCHARG) .OR. (XION.LE.1.D-30)) THEN
         EFUNC=0.D0
       ELSE
         DUM=6.D0*A*SQRT(XION)
         XIJ(1)=ICHARG*JCHARG*DUM
         XIJ(2)=ICHARG**2*DUM
         XIJ(3)=JCHARG**2*DUM
         DO 1 I=1,3
           J0(I)=XIJ(I)/(4.D0+4.581D0*XIJ(I)**(-7.238D-1)*EXP(-1.2D-2*
     >           XIJ(I)**5.28D-1))
1        CONTINUE
         EFUNC=ICHARG*JCHARG/(4.D0*XION)*(J0(1)-5.D-1*J0(2)-5.D-1*J0(3))
       ENDIF
      END
C
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
      FUNCTION EDFUNC(ICHARG,JCHARG,A,XION,ETHETA)
       IMPLICIT REAL*8(A-H,O-Z)
       REAL*8 XIJ(3),J1(3)
C
       IF((ICHARG.EQ.JCHARG) .OR. (XION.LE.1.D-30)) THEN
          EDFUNC=0.D0
       ELSE
          DUM1=6.D0*A*SQRT(XION)
          XIJ(1)=ICHARG*JCHARG*DUM1
          XIJ(2)=ICHARG**2*DUM1
          XIJ(3)=JCHARG**2*DUM1
          DO 1 I=1,3
            DUM=-1.2D-2*XIJ(I)**5.28D-1
            J1(I)=(4.D0+4.581D0*XIJ(I)**(-7.238D-1)*EXP(DUM)*(1.D0+
     >            7.238D-1+1.2D-2*5.28D-1*XIJ(I)**5.28D-1))/(4.D0+
     >            4.581D0*XIJ(I)**(-7.238D-1)*EXP(DUM))**2
1         CONTINUE
          EDFUNC=ICHARG*JCHARG/(8.D0*XION**2)*(XIJ(1)*J1(1)-5.D-1*
     >           XIJ(2)*J1(2)-5.D-1*XIJ(3)*J1(3))
     >           -ETHETA/XION
       ENDIF
      END
C
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
      FUNCTION FUNCI(NC,NA,xCAT,xAN,CATchrg,ANchrg)
       IMPLICIT REAL*8(A-H,O-Z)
C
       PARAMETER(NCmax=3,NAmax=5)
       REAL*8 CATchrg(NCmax),ANchrg(NAmax),xCAT(NCmax),xAN(NAmax)
C
C -------------------------------------------------------------
C                                                             -
C Calculates the ionic strength of the solution.              -
C                                                             -
C Used by: subr ANION, CATION & OSMC                          -
C                                                             -
C -------------------------------------------------------------
C
       FUNCI=0.D0
       DO 1 I=1,NC
         FUNCI=FUNCI+5.0D-1*xCAT(I)*CATchrg(I)**2
1      CONTINUE
       DO 2 J=1,NA
         FUNCI=FUNCI+5.0D-1*xAN(J)*ANchrg(J)**2
 2     CONTINUE
C
      END
C
C +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
      FUNCTION FF(xCAT,xAN,CATchrg,ANchrg)
      IMPLICIT REAL*8(A-H,O-Z)
C
      PARAMETER(NCmax=3,NAmax=5,NNmax=1)
      REAL*8 CATchrg(NCmax),ANchrg(NAmax),xCAT(NCmax),xAN(NAmax)
C
      COMMON/SPEC1ES/NC,NA,NN,nCAT(NCmax),nAN(NAmax),nNEUT(NNmax)
C
C -------------------------------------------------------------
C                                                             -
C Calculates the mole fraction function F                     -
C                                                             -
C -------------------------------------------------------------
C
       FF=0.D0
       DO 1 I=1,NC
         FF=FF+0.5D0*xCAT(I)*CATchrg(I)
1      CONTINUE
       DO 2 J=1,NA
         FF=FF+0.5D0*xAN(J)*ANchrg(J)
 2     CONTINUE
C
       FF=1.D0/FF
C
      END
C
C +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
      FUNCTION Ec(CATref,xCAT,CATchrg)
      IMPLICIT REAL*8(A-H,O-Z)
C
      PARAMETER(NCmax=3,NAmax=5,NNmax=1)
      REAL*8 CATchrg(NCmax),xCAT(NCmax)
      INTEGER CATref
      COMMON/SPEC1ES/NC,NA,NN,nCAT(NCmax),nAN(NAmax),nNEUT(NNmax)
C
C -------------------------------------------------------------
C                                                             -
C Calculates the function Ec                                  -
C                                                             -
C -------------------------------------------------------------
C
       DUM=0.D0
       DO 1 I=1,NC
         DUM=DUM+xCAT(I)*CATchrg(I)
1      CONTINUE
C
       EC=xCAT(CATref)*CATchrg(CATref)/DUM
C
       END
C
C +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
      FUNCTION Ea(ANref,xAN,ANchrg)
      IMPLICIT REAL*8(A-H,O-Z)
C
      PARAMETER(NCmax=3,NAmax=5,NNmax=1)
      REAL*8 ANchrg(NAmax),xAN(NAmax)
      INTEGER ANref
      COMMON/SPEC1ES/NC,NA,NN,nCAT(NCmax),nAN(NAmax),nNEUT(NNmax)
C
C -------------------------------------------------------------
C                                                             -
C Calculates the function Ea                                  -
C                                                             -
C -------------------------------------------------------------
C
       DUM=0.D0
       DO 1 I=1,NA
         DUM=DUM+xAN(I)*ANchrg(I)
1      CONTINUE
C
       EA=xAN(ANref)*ANchrg(ANref)/DUM
C
       END
C
C +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
      REAL*8 FUNCTION MEc(nSELECT,CATref,xCAT,CATchrg)
      IMPLICIT REAL*8(A-H,O-Z)
C
      PARAMETER(NCmax=3,NAmax=5,NNmax=1)
      REAL*8 CATchrg(NCmax),xCAT(NCmax)
      INTEGER CATref
      COMMON/SPEC1ES/NC,NA,NN,nCAT(NCmax),nAN(NAmax),nNEUT(NNmax)
C
C -------------------------------------------------------------
C                                                             -
C Calculates the function MEc                                 -
C                                                             -
C -------------------------------------------------------------
C
       DUM=0.D0
       DO 1 I=1,NC
         DUM=DUM+xCAT(I)*CATchrg(I)
1      CONTINUE
C
       IF(nSELECT.EQ. CATref) THEN
C    ..first the c = M case:
         MEc=CATchrg(nSELECT)/DUM*(1.D0-Ec(nSELECT,xCAT,CATchrg))
       ELSE
C    ..now the c NE M case:
         MEc=-CATchrg(nSELECT)*Ec(CATref,xCAT,CATchrg)/DUM
       ENDIF
C
       END
C
C +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
      REAL*8 FUNCTION XEa(nSELECT,ANref,xAN,ANchrg)
      IMPLICIT REAL*8(A-H,O-Z)
C
      PARAMETER(NCmax=3,NAmax=5,NNmax=1)
      REAL*8 ANchrg(NAmax),xAN(NAmax)
      INTEGER ANref
      COMMON/SPEC1ES/NC,NA,NN,nCAT(NCmax),nAN(NAmax),nNEUT(NNmax)
C
C -------------------------------------------------------------
C                                                             -
C Calculates the function XEa                                 -
C                                                             -
C -------------------------------------------------------------
C
       DUM=0.D0
       DO 1 J=1,NA
         DUM=DUM+xAN(J)*ANchrg(J)
1      CONTINUE
C
       IF(nSELECT .EQ. ANref) THEN
C    ..first the a = X case:
         XEa=ANchrg(nSELECT)/DUM*(1.D0-Ea(nSELECT,xAN,ANchrg))
       ELSE
C    ..now the a NE X case:
         XEa=-ANchrg(nSELECT)*Ea(ANref,xAN,ANchrg)/DUM
       ENDIF
C
       END
C
C +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
      INTEGER FUNCTION PICK(NREF,INDXSP)
      IMPLICIT REAL*8 (A-H,O-Z)  
C
      PARAMETER(NCmax=3,NAmax=5,NNmax=1)
C
      COMMON/SPEC1ES/NC,NA,NN,nCAT(NCmax),nAN(NAmax),nNEUT(NNmax)
C
      PICK=-99
      IF(INDXSP.EQ.0) THEN
        DO 1 I=1,NN
          IF(NNEUT(I).EQ.NREF) THEN 
            PICK=I
            RETURN
          ENDIF
1       CONTINUE
      ELSE IF(INDXSP.GE.1) THEN
        DO 2 J=1,NC
          IF(NCAT(J).EQ.NREF) THEN
            PICK=J
            RETURN
          ENDIF
2       CONTINUE
      ELSE IF(INDXSP.LE.-1) THEN
        DO 3 K=1,NA
          IF(NAN(K).EQ.NREF) THEN
            PICK=K
            RETURN
          ENDIF
3       CONTINUE
      ENDIF
C
      END
C
C *******************************************************************
C
      LOGICAL FUNCTION PEEK(NREF,INDXSP)
      IMPLICIT REAL*8 (A-H,O-Z)
C
      PARAMETER(NCmax=3,NAmax=5,NNmax=1)
C
      COMMON/SPEC1ES/NC,NA,NN,nCAT(NCmax),nAN(NAmax),nNEUT(NNmax)
C
      PEEK=.FALSE.
      IF(INDXSP.EQ.0) THEN
        DO 1 I=1,NN
          IF(NNEUT(I).EQ.NREF) THEN
            PEEK=.TRUE.
            RETURN
          ENDIF
1       CONTINUE
      ELSE IF(INDXSP.GE.1) THEN
        DO 2 J=1,NC
          IF(NCAT(J).EQ.NREF) THEN
            PEEK=.TRUE.
            RETURN
          ENDIF
2       CONTINUE
      ELSE IF(INDXSP.LE.-1) THEN
        DO 3 K=1,NA
          IF(NAN(K).EQ.NREF) THEN
            PEEK=.TRUE.
            RETURN
          ENDIF
3       CONTINUE
      ENDIF
C
      END
C
C *******************************************************************
C
C                         end of program
C
C                     NOTES ON DATA SOURCES
C
C The model uses the equations of Pitzer, Simonson and Clegg to describe 
C the activity coefficients of the aqueous ions, see:
C 
C S. L. Clegg, K. S. Pitzer and P. Brimblecombe (1992) Thermodynamics 
C of Multicomponent, Miscible, Ionic Solutions. II. Mixtures Including 
C Unsymmetrical Electrolytes. J. Phys. Chem. 96, 9470-9479; 98, 1368;
C 99, 6755.
C 
C Briefly, the activity coefficient equations contain empirical 
C parameters for the interaction of pairs and triplets of species in the 
C solution. These parameters are determined (as functions of temperature) 
C by fitting to thermodynamic data for activities or thermal properties. 
C The present model 2 is considered a 'draft' only, since it is expected 
C that further extensions to cover supersaturated solutions will be made, 
C also more attention will be given to the variation of the properties of 
C H2SO4-(NH4)2SO4-H2O with temperature. Data for this important mixture 
C are sparse, and the model parameterisation away from 298 K is rendered 
C difficult by the fact that HSO4 - SO4 speciation is not well 
C constrained by measurements. The data on which the model parameterisations 
C are based are summarised below:
C
C  -----------
C | H2SO4-H2O | activities, osmotic coefficients, enthalpies, heat 
C  -----------  capacities, saturations with respect to ice and solid 
C phase hydrates from 328 K to <200K. Maximum concentration 40 mol kg-1. 
C See [1].
C
C  ----------
C | HNO3-H2O | activities, osmotic coefficients, vapour pressures, freezing 
C  ----------  points, and thermal data, originally for the entire 
C composition range and temperatures from about 220 K to >373 K [2]. The 
C present representation of activity coefficients has been limited to 
C molalities < 88 mol kg-1, and temperatures below 330 K [3].
C
C  ---------------
C | (NH4)2SO4-H2O | osmotic coefficients, freezing points and boiling 
C  ---------------  points, thermal data, and electrodynamic balance 
C measurements (yielding water activities for supersaturated solutions) 
C from about 255 K to >373 K [4]. Note that the treatment for 
C supersaturated solutions (to about 25 mol kg-1) is not likely to be 
C valid above 323 K. Also, a further model for subsaturated solutions 
C has recently been produced, which incorporates newer data [5].
C
C  ------------
C | NH4NO3-H2O | osmotic coefficients up to the saturation concentration 
C  ------------  of 26 mol kg-1 at 25oC [6], electrodynamic balance 
C water activities to about 100 mol kg-1 at 298 K [7,8], apparent molar 
C enthalpies [9,10] and heat capacities [10-14]. Some high concentration 
C enthalpies were estimated from vapour pressure data of Othmer and 
C Frolich [15]. Where the heat capacity data were for t>25oC (and 
C generally high concentration) the measurements were used to estimate 
C apparent molar heat capacities at 25oC. Only values at this temperature 
C were fitted. The fit should give an adequate representation of aqueous 
C properties for tropospheric temperatures, even for highly 
C supersaturated solutions. Saturation with respect to NH4NO3 is 
C represented by the model from 256.3 K (eutectic) to 328 K. 
C
C  ----------------
C | HNO3-H2SO4-H2O | extensive vapour pressure and solid phase saturation 
C  ----------------  data at and below 273.15 K, see [3].
C
C  ---------------
C | HCl-H2SO4-H2O | vapour pressure data for a wide range of temperatures 
C  ---------------  showed that the mixture parameters in the model can 
C be set to zero. See [3].
C
C  --------------
C | HCl-HNO3-H2O | at the time of preparation of our model (1), there 
C  --------------  were no satisfactory data for this system, and the ternary
C interaction parameters for these ions were set to zero. Measurements by 
C Elrod et al. (Faraday Disc. 100, 269-278, 1995) for HCl-HNO3-H2SO4-H2O 
C suggest that the model will nonetheless give reasonable results for this 
C system.
C
C  ---------------------
C | H2SO4-(NH4)2SO4-H2O | osmotic coefficients (25oC and 50oC),  
C  ---------------------  electrodynamic balance measurements (room 
C temperature) for supersaturated solutions, and solid phase saturations 
C from 273 K to 323 K. Most of these data are referenced in [37]. The fit 
C used here differs from that in [37] as follows: (NH4)2SO4-H2O 
C parameters are from [4]; the fit is extended here to temperatures other 
C than 25oC using measurements not available at the time of the original 
C work [5]. 
C
C  -----------------
C | HNO3-NH4NO3-H2O | solid phase saturations (wrt NH4NO3) from 273 to 
C  -----------------  303 K and 44 mol kg-1 HNO3 [38].
C
C  ----------------------
C | (NH4)2SO4-NH4NO3-H2O | electrodynamic balance data at room 
C  ----------------------  temperature for supersaturated solutions [8], 
C recently revised by Clegg and Chan (unpublished). Also, solid phase 
C saturations (four solid phases) from 273 K to 313 K taken from [38] 
C and original literature sources. 
C
C  --------------------
C | NH4HSO4-NH4NO3-H2O | solid phase saturations (three phases) at 298 K, 
C  --------------------  plus values for the double salt NH4NO3.NH4HSO4 
C from 253 K to 323 K [38]. 
C
C
C                           Bibliography
C
C 1.  Clegg S. L. and Brimblecombe P. (1995) Application of a 
C multicomponent thermodynamic model to activities and thermal properties 
C of 0 - 40 mol kg-1 aqueous sulphuric acid from <200 K to 328 K. J. 
C Chem. Eng. Data 40, 43-64. 
C
C 2.  Clegg S. L. and Brimblecombe P. (1990) Equilibrium partial 
C pressures, and mean activity and osmotic coefficients of 0-100% nitric 
C acid as a function of temperature. J. Phys. Chem. 94, 5369-5380; 96, 
C 6854. 
C
C 3.  Carslaw K. S., Clegg S. L. and Brimblecombe P. (1995) A 
C thermodynamic model of the system HCl-HNO3-H2SO4-H2O, including 
C solubilities of HBr, from <200 K to 328 K. J. Phys. Chem. 99, 11557-
C 11574. 
C
C 4.  Clegg S. L., Ho S. S., Chan C. K. and Brimblecombe P. (1995) 
C Thermodynamic properties of aqueous (NH4)2SO4 to high supersaturation 
C as a function of temperature. J. Chem. Eng. Data 40, 1079-1090. 
C
C 5.  Clegg S. L., Milioto S. and Palmer D. A. Osmotic and activity 
C coefficients of aqueous (NH4)2SO4 as a function of temperature, and 
C (NH4)2SO4-H2SO4 mixtures at 298.15 K and 323.15 K. Submitted to J. 
C Chem. Eng. Data 
C
C 7.  Tang I. N. and Munkelwitz H. R. (1994) Water activities, densities, 
C and refractive indices of aqueous sulphates and sodium nitrate droplets 
C of atmospheric importance. J. Geophys. Res. 99, 18801-18808. 
C
C 8.  Chan C. K., Flagan R. C. and Seinfeld J. H. (1992) Water activities 
C of NH4NO3/(NH4)2SO4 solutions. Atmos. Env. 26A, 1661-1673. 
C
C 9.  Vanderzee C. E., Waugh D. H. and Haas N. C. (1980) Enthalpies of 
C dilution and relative apparent model enathalpies of aqueous ammonium 
C nitrate. The case of a weakly hydrolysed (dissociated) salt. J. Chem.  
C Thermo. 12, 21-25. 
C
C 10. Parker V. B. (1965) National Bureau of Standards Reference Data 
C Series. Thermal Properties of Aqueous Uni-univalent Electrolytes. U.S. 
C Gov. Printing Office, Washington. 
C
C 11. Roux A., Musbally G. M., Perron G., Desnoyers J. E., Singh P. P., 
C Woolley E. M. and Hepler L. G. (1978) Apparent molal heat capacities 
C and volumes of aqueous electrolytes at 25oC: NaClO3, NaClO4, NaNO3, 
C NaBrO3, NaIO3, KClO3, KBrO3, KIO3, NH4NO3, NH4Cl, and NH4ClO4. Can. J. 
C Chem. 56, 24-28. 
C
C 12. Epikhin Yu. A., Bazlova I. V. and Karapet'yants M. Kh. (1977) 
C Changes in the volume and heat capacity in aqueous salt solutions. IV. 
C The ammonium chloride-ammonium nitrate-water system. Russ. J. Phys. 
C Chem. 51, 676-677. 
C
C 13. Gladushko V. I., Bochenko G. A., Prokof'eva G. N., Privalko V. P. 
C and Vinarcik J. (1985) Heat capacity of the ternary ammonium nitrate-
C nitric acid-water system. Inzh. Fiz. Zh. 48, 90-91. 
C
C 14. Sorina G. A., Blinova M. B. and Tsekhanskaya Yu. V. (1983) Boiling 
C points of aqueous ammonium nitrate solutions under pressure. J. Appl. 
C Chem. USSR 56, 1754-1757. 
C
C 15. Othmer D. F. and Frohlich G. J. (1960) Correlating vapour pressures 
C and heats of solution for the ammonium nitrate-water system: the 
C enthalpy concentration diagram. AIChE J. 6, 210-214. 
C
C 37. Clegg S. L. and Brimblecombe P. (1995) A generalised multicomponent 
C thermodynamic model applied to the (NH4)2SO4-H2SO4-H2O system to high 
C supersaturation and low relative humidity at 298.15 K. J. Aerosol Sci. 
C 26, 19-38. 
C
C 38. Silcock H. L. (1979) In Solubilities of Inorganic and Organic 
C Compounds Vol. 3, Pergamon, Oxford. 
C
C ----------------------- END OF FILE ------------------------------------
