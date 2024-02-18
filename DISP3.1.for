      SUBROUTINE  DISP(U,KSTEP,KINC,TIME,NODE,NOEL,JDOF,COORDS)
      INCLUDE 'ABA_PARAM.INC'
c
      DIMENSION U(3),TIME(2),COORDS(3)
      double precision tim
      double precision n_period
      integer::period
      parameter ( T  = 2.0D0 ,
     *           AMP = 0.08     )
      tim=time(1)
      period = tim / T
      n_period = mod(tim,T)
      T_1 = 0.25 * T 
      T_2 = 0.75 * T       
      if (period.LE.500) then
         if(n_period .LE. T_1) then
              U(1) = AMP * (1 * (tim - period * T) + 0)
              U(2) = AMP * 1 * tim
              U(3) = 0
         else if(n_period .gt. T_1 .and. n_period .LE. T_2) then      
              U(1) = AMP * (-1 * (tim - period * T) + 1.0)  
              U(2) = AMP * -1 * tim
              U(3) = 0
         else if(n_period .gt. T_2 .and. n_period .LE. T) then 
              U(1) = AMP * (1 * (tim - period * T) - 2.0)
              U(2) = AMP * 1 * tim
              U(3) = 0
         end if
      else
          call xit
          RETURN
      end if
      RETURN
      END
c______________________________________________________________________
C     USER INPUT FOR ADAPTIVE MESH CONSTRAINT
C
      SUBROUTINE UMESHMOTION(UREF,ULOCAL,NODE,NNDOF,
     $     LNODETYPE,ALOCAL,NDIM,TIME,DTIME,PNEWDT,
     $     KSTEP,KINC,KMESHSWEEP,JMATYP,JGVBLOCK,LSMOOTH)
C
      include 'ABA_PARAM.INC'

      CHARACTER*80 PARTNAME
      DIMENSION ARRAY(3), JPOS(15)
      DIMENSION ULOCAL(NDIM)
      DIMENSION ALOCAL(NDIM,*), TIME(2)
      DIMENSION JMATYP(*), JGVBLOCK(*)
      parameter (NELEMMAX=600000)
      DIMENSION JELEMLIST(NELEMMAX), JELEMTYPE(NELEMMAX)
      parameter (nUM=600000) 
      common /wear/ DEPTH(nUM),DEPTH1(nUM),lvalidinc,lvalidsweep
      data lvalidinc /-1/ 
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc      
      NELEMS = NELEMMAX
      JRCD = 0 
      CALL GETNODETOELEMCONN(NODE,NELEMS,JELEMLIST,JELEMTYPE,
     $     JRCD,JGVBLOCK)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      luseendisp = 1
      IF (KINC.EQ.1.AND.kmeshsweep.eq.0.and.lvalidinc.eq.-1) then
          lvalidinc = 2
          lvalidsweep = 1
          luseendisp = 0    
      else if (kinc .eq. lvalidinc) then
          do k1 = 1,nUM
             DEPTH1(k1) = DEPTH(k1)
             DEPTH(k1) = 0.0d0 
          end do
        lvalidinc = kinc + 1
        lvalidsweep = 1
      end if
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc    
      if (kmeshsweep.eq.lvalidsweep) then
        do k1 = 1,nUM
           DEPTH(k1) = 0.0d0 
        end do
        lvalidsweep = kmeshsweep + 1          
      end if
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc  
C      LOCNUM = 0
      PARTNAME = ' '     
      CALL GETPARTINFO(NODE,0,PARTNAME,LOCNUM,JRCD)
      jtyp = 0
      CALL GETVRMAVGATNODE(NODE,JTYP,'CSTRESS',ARRAY,JRCD,
     $     JELEMLIST,NELEMS,JMATYP,JGVBLOCK)    
      CPRESS = ARRAY(1)
      CSHEAR = SQRT(ARRAY(2)**2+ARRAY(3)**2)
      CALL GETVRMAVGATNODE(NODE,JTYP,'CDISP',ARRAY,JRCD,
     $     JELEMLIST,NELEMS,JMATYP,JGVBLOCK)    
      CSLIP = SQRT(ARRAY(2)**2+ARRAY(3)**2) 
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      UREF = 1.0D-7
      DEPTH(LOCNUM) = UREF * CSLIP * CPRESS
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      if (luseendisp.ne.0) then
         ULOCAL(3) = ULOCAL(3) - DEPTH1(LOCNUM)
      end if
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc    
      RETURN
      END
c-----------------------------------------------------------------------------
      subroutine fric_coef (
C Write only - 
     *   fCoef, fCoefDeriv, 
C Read only - 
     *   nBlock, nProps, nTemp, nFields, 
     *   jFlags, rData, 
     *   surfInt, surfSlv, surfMst, 
     *   props, slipRate, pressure, 
     *   tempAvg, fieldAvg )
C
      include 'aba_param.inc'
C
      dimension fCoef(nBlock), 
     *   fCoefDeriv(nBlock,3),
     *   props(nProps),
     *   slipRate(nBlock),
     *   pressure(nBlock),
     *   tempAvg(nBlock), 
     *   fieldAvg(nBlock,nFields)
C
      parameter( iKStep   = 1,
     *           iKInc    = 2,
     *           nFlags   = 2 )
C
      parameter( iTimStep = 1,
     *           iTimGlb  = 2,
     *           iDTimCur = 3,
     *           nData    = 3 )
C
      dimension jFlags(nFlags), rData(nData)
C
      character*80 surfInt, surfSlv, surfMst 
C
         xMuk = props(1)
         xMus = props(2)
         beta = props(3)
      DO  k = 1, nBlock
         fCoef(k)=xMuk + (xMus - xMuk) * exp(-beta * slipRate(k))
      END DO

      return
      end
