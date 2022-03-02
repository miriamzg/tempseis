
      subroutine calcmisfit( predicted_data, observed_data, ndata,
     1       weight, nwaves, totmisfit, fname, amp_obs_real,
     1       amp_obs_imag, write_best )

      include 'homti_param.inc'

      complex*8 observed_data(maxdata,maxwave),
     1          predicted_data(maxdata,maxwave),
     1          observed_phase,
     1          predicted_phase

      real*4 weight(maxwave),
     1       abs_value_obs, abs_value_synt
      real*4 real_mft, imag_mft
      real*4 totmisfit, misfitval, aval_real, aval_imag, aval,
     1       mft_abs, mft_fft
      real*4 mft_phase, mft_absolute, amplitude
      real*4 totmisfit_P,totmisfit_S,totmisfit_W
      integer n_P, n_S, n_W

      integer ndata(maxwave), nwaves
      character*40 fname(maxwave)
      character*40 instring, string_r, string1, string2, delim
      character wavetype

      real*4 amp_obs_real(maxwave),
     1           amp_obs_imag(maxwave)
      real*4 a(15)
      real*4 m
      real*4 misfit_P_array(maxwave),misfit_S_array(maxwave),
     1           misfit_W_array(maxwave)
      real*4 median_P, median_S, median_W
      logical write_best

      pi = 4.D0*DATAN(1.D0)

!          misfit between observed and predicted
      totmisfit=0.0

      totmisfit_P = 0.0
      totmisfit_S = 0.0
      totmisfit_W = 0.0

      n_P = 0
      n_S = 0
      n_W = 0

      do iw=1,nwaves
        instring = fname(iw)
        call reverse_string(instring, string_r)
        delim = "_"
        call split_string(string_r, string1, string2, delim)
        wavetype = string2(1:1)

        aval=0.0
        aval_real = 0.0
        aval_imag = 0.0
        mft_fft = 0.0
        mft_abs = 0.0
        mft_phase = 0.0
        mft_absolute = 0.0

        do i=1,ndata(iw)
!			----------------------------
!			Original code, do not edit
          real_mft = (real(predicted_data(i,iw)) - 
     1                real(observed_data(i,iw)))**2
          imag_mft = (imag(predicted_data(i,iw)) - 
     1                imag(observed_data(i,iw)))**2
!			----------------------------
          aval = aval + (real_mft + imag_mft)

        end do

        aval = weight(iw) * aval / (amp_obs_real(iw) + amp_obs_imag(iw))
        write(123,*) "STATION: ", fname(iw)(15:26) , " Misfit: ", aval

        if (write_best) then
          write(111,*) "STATION: ", fname(iw)(15:24) , " Misfit: ", aval
        end if

!		*** ORGINAL LINE ***
        if (wavetype .eq. "P") then
          totmisfit_P = totmisfit_P + aval
          n_P = n_P + 1
          misfit_P_array(n_P) = aval


        end if
        if (wavetype .eq. "S") then
          totmisfit_S = totmisfit_S + aval
          n_S = n_S + 1
          misfit_S_array(n_S) = aval
        end if
        if (wavetype .eq. "W") then
          totmisfit_W = totmisfit_W + aval
          n_W = n_W + 1
          misfit_W_array(n_W) = aval
        end if
      end do

      if (n_P .eq. 0) then
        median_P = 0.
      else
        call median(misfit_P_array,n_P,0,m)
        median_P = m
      end if
!---------
      if (n_S .eq. 0) then
        median_S = 0.
      else
        call median(misfit_S_array,n_S,0,m)
        median_S = m
      end if
!---------
      if (n_W .eq. 0) then
        median_W = 0
      else
        call median(misfit_W_array,n_W,0,m)
        median_W = m
      end if
!---------


      totmisfit = median_P + median_S + median_W

      write(*,*) "median P: ", median_P, " median S: ", median_S, 
     1           " median W: ", median_W
      write(*,*) "n P: ", n_P, " n S: ", n_S, " n W: ", n_W

      return
      end

!------------------------------------------------------------------------
      subroutine reverse_string(instring, outstring)
      character*40 instring, outstring
      character temp
      integer i,length


      length = len_trim(instring) ! ignores trailing blanks.
!                 use len(string) to reverse those as well

      outstring = instring
      do i = 1, length/2
        temp = outstring(i:i)
        outstring(i:i) = outstring(length+1-i:length+1-i)
        outstring(length+1-i:length+1-i) = temp
      end do

      end subroutine reverse_string


!------------------------------------------------------------------------
      ! split a string into 2 either side of a delimiter token
      SUBROUTINE split_string(instring, string1, string2, delim)
      CHARACTER(40) :: instring,delim
      CHARACTER(40),INTENT(OUT):: string1,string2
      INTEGER :: index

      instring = adjustl(instring)

      index = SCAN(instring, delim)
      string1 = instring(1:index-1)
      string2 = instring(index+1:)

      END SUBROUTINE split_string

!------------------------------------------------------------------------
      SUBROUTINE MEDIAN(X,N,IWRITE,XMED)
!
!     PURPOSE--THIS SUBROUTINE COMPUTES THE
!              SAMPLE MEDIAN
!              OF THE DATA IN THE INPUT VECTOR X.
!              THE SAMPLE MEDIAN = THAT VALUE SUCH THAT HALF THE
!              DATA SET IS BELOW IT AND HALF ABOVE IT.
!     INPUT  ARGUMENTS--X      = THE SINGLE PRECISION VECTOR OF
!                                (UNSORTED OR SORTED) OBSERVATIONS.
!                     --N      = THE INTEGER NUMBER OF OBSERVATIONS
!                                IN THE VECTOR X.
!                     --IWRITE = AN INTEGER FLAG CODE WHICH
!                                (IF SET TO 0) WILL SUPPRESS
!                                THE PRINTING OF THE
!                                SAMPLE MEDIAN
!                                AS IT IS COMPUTED;
!                                OR (IF SET TO SOME INTEGER
!                                VALUE NOT EQUAL TO 0),
!                                LIKE, SAY, 1) WILL CAUSE
!                                THE PRINTING OF THE
!                                SAMPLE MEDIAN
!                                AT THE TIME IT IS COMPUTED.
!     OUTPUT ARGUMENTS--XMED   = THE SINGLE PRECISION VALUE OF THE
!                                COMPUTED SAMPLE MEDIAN.
!     OUTPUT--THE COMPUTED SINGLE PRECISION VALUE OF THE
!             SAMPLE MEDIAN.
!     PRINTING--NONE, UNLESS IWRITE HAS BEEN SET TO A NON-ZERO
!               INTEGER, OR UNLESS AN INPUT ARGUMENT ERROR
!               CONDITION EXISTS.
!     RESTRICTIONS--THE MAXIMUM ALLOWABLE VALUE OF N
!                   FOR THIS SUBROUTINE IS 15000.
!     OTHER DATAPAC   SUBROUTINES NEEDED--SORT.
!     FORTRAN LIBRARY SUBROUTINES NEEDED--NONE.
!     MODE OF INTERNAL OPERATIONS--SINGLE PRECISION.
!     LANGUAGE--ANSI FORTRAN.
!     REFERENCES--KENDALL AND STUART, THE ADVANCED THEORY OF
!                 STATISTICS, VOLUME 1, EDITION 2, 1963, PAGE 326.
!               --KENDALL AND STUART, THE ADVANCED THEORY OF
!                 STATISTICS, VOLUME 2, EDITION 1, 1961, PAGE 49.
!               --DAVID, ORDER STATISTICS, 1970, PAGE 139.
!               --SNEDECOR AND COCHRAN, STATISTICAL METHODS,
!                 EDITION 6, 1967, PAGE 123.
!               --DIXON AND MASSEY, INTRODUCTION TO STATISTICAL
!                 ANALYSIS, EDITION 2, 1957, PAGE 70.
!     WRITTEN BY--JAMES J. FILLIBEN
!                 STATISTICAL ENGINEERING LABORATORY (205.03)
!                 NATIONAL BUREAU OF STANDARDS
!                 WASHINGTON, D. C. 20234
!                 PHONE:  301-921-2315
!     ORIGINAL VERSION--JUNE      1972.
!     UPDATED         --SEPTEMBER 1975.
!     UPDATED         --NOVEMBER  1975.
!     UPDATED         --FEBRUARY  1976.
!
!---------------------------------------------------------------------
      DIMENSION X(1)
      DIMENSION Y(15000)
      COMMON /BLOCK2/ WS(15000)
      EQUIVALENCE (Y(1),WS(1))
      IPR=6
      IUPPER=15000
!
!     CHECK THE INPUT ARGUMENTS FOR ERRORS
!
      IF(N.LT.1.OR.N.GT.IUPPER)GOTO50
      IF(N.EQ.1)GOTO55
      HOLD=X(1)
      DO60I=2,N
        IF(X(I).NE.HOLD)GOTO90
  60  CONTINUE
      WRITE(IPR, 9)HOLD
      XMED=X(1)
      GOTO101
  50  WRITE(IPR,17)IUPPER
      WRITE(IPR,47)N
      RETURN
  55  WRITE(IPR,18)
      XMED=X(1)
      GOTO101
   90 CONTINUE
    9 FORMAT("***** NON-FATAL DIAGNOSTIC--THE FIRST  INPUT ARGUMENT", 
     1"(A VECTOR) TO THE MEDIAN SUBROUTINE HAS ALL ELEMENTS = ",
     1 E15.8," *****")
   17 FORMAT("***** FATAL ERROR--THE SECOND INPUT ARGUMENT TO THE",
     1"MEDIAN SUBROUTINE IS OUTSIDE THE ALLOWABLE", I6,
     1"INTERVAL *****")
   18 FORMAT("***** NON-FATAL DIAGNOSTIC--THE SECOND INPUT ARGUMENT TO",
     1" THE MEDIAN SUBROUTINE HAS THE VALUE 1 *****")
   47 FORMAT( 35H***** THE VALUE OF THE ARGUMENT IS ,I8   ,6H *****)
!
!-----START POINT-----------------------------------------------------
!
      CALL SORT(X,N,Y)
      IFLAG=N-(N/2)*2
      NMID=N/2
      NMIDP1=NMID+1
      IF(IFLAG.EQ.0)XMED=(Y(NMID)+Y(NMIDP1))/2.0
      IF(IFLAG.EQ.1)XMED=Y(NMIDP1)

  101 IF(IWRITE.EQ.0)RETURN
      WRITE(IPR,999)
      WRITE(IPR,105)N,XMED
  105 FORMAT(1H ,25HTHE SAMPLE MEDIAN OF THE ,
     1       I6,17H OBSERVATIONS IS ,E5.8)
  999 FORMAT(1H )
      RETURN
      END



!------------------------------------------------------------------------
      SUBROUTINE SORT(X,N,Y)
!
!     PURPOSE--THIS SUBROUTINE SORTS (IN ASCENDING ORDER)
!              THE N ELEMENTS OF THE SINGLE PRECISION VECTOR X
!              AND PUTS THE RESULTING N SORTED VALUES INTO THE
!              SINGLE PRECISION VECTOR Y.
!     INPUT  ARGUMENTS--X      = THE SINGLE PRECISION VECTOR OF
!                                OBSERVATIONS TO BE SORTED.
!                     --N      = THE INTEGER NUMBER OF OBSERVATIONS
!                                IN THE VECTOR X.
!     OUTPUT ARGUMENTS--Y      = THE SINGLE PRECISION VECTOR
!                                INTO WHICH THE SORTED DATA VALUES
!                                FROM X WILL BE PLACED.
!     OUTPUT--THE SINGLE PRECISION VECTOR Y
!             CONTAINING THE SORTED
!             (IN ASCENDING ORDER) VALUES
!             OF THE SINGLE PRECISION VECTOR X.
!     PRINTING--NONE UNLESS AN INPUT ARGUMENT ERROR CONDITION EXISTS.
!     RESTRICTIONS--THE DIMENSIONS OF THE VECTORS IL AND IU
!                   (DEFINED AND USED INTERNALLY WITHIN
!                   THIS SUBROUTINE) DICTATE THE MAXIMUM
!                   ALLOWABLE VALUE OF N FOR THIS SUBROUTINE.
!                   IF IL AND IU EACH HAVE DIMENSION K,
!                   THEN N MAY NOT EXCEED 2**(K+1) - 1.
!                   FOR THIS SUBROUTINE AS WRITTEN, THE DIMENSIONS
!                   OF IL AND IU HAVE BEEN SET TO 36,
!                   THUS THE MAXIMUM ALLOWABLE VALUE OF N IS
!                   APPROXIMATELY 137 BILLION.
!                   SINCE THIS EXCEEDS THE MAXIMUM ALLOWABLE
!                   VALUE FOR AN INTEGER VARIABLE IN MANY COMPUTERS,
!                   AND SINCE A SORT OF 137 BILLION ELEMENTS
!                   IS PRESENTLY IMPRACTICAL AND UNLIKELY,
!                   THEN THERE IS NO PRACTICAL RESTRICTION
!                   ON THE MAXIMUM VALUE OF N FOR THIS SUBROUTINE.
!                   (IN LIGHT OF THE ABOVE, NO CHECK OF THE
!                   UPPER LIMIT OF N HAS BEEN INCORPORATED
!                   INTO THIS SUBROUTINE.)
!     OTHER DATAPAC   SUBROUTINES NEEDED--NONE.
!     FORTRAN LIBRARY SUBROUTINES NEEDED--NONE.
!     MODE OF INTERNAL OPERATIONS--SINGLE PRECISION.
!     LANGUAGE--ANSI FORTRAN.
!     COMMENT--THE SMALLEST ELEMENT OF THE VECTOR X
!              WILL BE PLACED IN THE FIRST POSITION
!              OF THE VECTOR Y,
!              THE SECOND SMALLEST ELEMENT IN THE VECTOR X
!              WILL BE PLACED IN THE SECOND POSITION
!              OF THE VECTOR Y, ETC.
!     COMMENT--THE INPUT VECTOR X REMAINS UNALTERED.
!     COMMENT--IF THE ANALYST DESIRES A SORT 'IN PLACE',
!              THIS IS DONE BY HAVING THE SAME
!              OUTPUT VECTOR AS INPUT VECTOR IN THE CALLING SEQUENCE.
!              THUS, FOR EXAMPLE, THE CALLING SEQUENCE
!              CALL SORT(X,N,X)
!              IS ALLOWABLE AND WILL RESULT IN
!              THE DESIRED 'IN-PLACE' SORT.
!     COMMENT--THE SORTING ALGORTHM USED HEREIN
!              IS THE BINARY SORT.
!              THIS ALGORTHIM IS EXTREMELY FAST AS THE
!              FOLLOWING TIME TRIALS INDICATE.
!              THESE TIME TRIALS WERE CARRIED OUT ON THE
!              UNIVAC 1108 EXEC 8 SYSTEM AT NBS
!              IN AUGUST OF 1974.
!              BY WAY OF COMPARISON, THE TIME TRIAL VALUES
!              FOR THE EASY-TO-PROGRAM BUT EXTREMELY
!              INEFFICIENT BUBBLE SORT ALGORITHM HAVE
!              ALSO BEEN INCLUDED--
!              NUMBER OF RANDOM        BINARY SORT       BUBBLE SORT
!               NUMBERS SORTED
!                N = 10                 .002 SEC          .002 SEC
!                N = 100                .011 SEC          .045 SEC
!                N = 1000               .141 SEC         4.332 SEC
!                N = 3000               .476 SEC        37.683 SEC
!                N = 10000             1.887 SEC      NOT COMPUTED
!     REFERENCES--CACM MARCH 1969, PAGE 186 (BINARY SORT ALGORITHM
!                 BY RICHARD C. SINGLETON).
!               --CACM JANUARY 1970, PAGE 54.
!               --CACM OCTOBER 1970, PAGE 624.
!               --JACM JANUARY 1961, PAGE 41.
!     WRITTEN BY--JAMES J. FILLIBEN
!                 STATISTICAL ENGINEERING LABORATORY (205.03)
!                 NATIONAL BUREAU OF STANDARDS
!                 WASHINGTON, D. C. 20234
!                 PHONE--301-921-2315
!     ORIGINAL VERSION--JUNE      1972.
!     UPDATED         --NOVEMBER  1975.
!
!---------------------------------------------------------------------
!
      DIMENSION X(1),Y(1)
      DIMENSION IU(36),IL(36)

      IPR=6
!
!     CHECK THE INPUT ARGUMENTS FOR ERRORS
!
      IF(N.LT.1)GOTO50
      IF(N.EQ.1)GOTO55
      HOLD=X(1)
      DO60I=2,N
        IF(X(I).NE.HOLD)GOTO90
   60 CONTINUE
      WRITE(IPR, 9)HOLD
      DO61I=1,N
        Y(I)=X(I)
   61 CONTINUE
      RETURN
   50 WRITE(IPR,15)
      WRITE(IPR,47)N
      RETURN
   55 WRITE(IPR,18)
      Y(1)=X(1)
      RETURN
   90 CONTINUE
    9 FORMAT("***** NON-FATAL DIAGNOSTIC--THE FIRST  INPUT ARGUMENT ",
     1"(A VECTOR) TO THE SORT   SUBROUTINE HAS ALL ELEMENTS = ",
     1 E15.8," *****")
   15 FORMAT("***** FATAL ERROR--THE SECOND INPUT ARGUMENT TO THE SORT",
     1" SUBROUTINE IS NON-POSITIVE *****")
   18 FORMAT("***** NON-FATAL DIAGNOSTIC--THE SECOND INPUT ARGUMENT ",
     1"TO THE SORT   SUBROUTINE HAS THE VALUE 1 *****")
   47 FORMAT("***** THE VALUE OF THE ARGUMENT IS ",I8   ," *****")
!
!-----START POINT-----------------------------------------------------
!
!     COPY THE VECTOR X INTO THE VECTOR Y
      DO100I=1,N
        Y(I)=X(I)
  100 CONTINUE
!
!     CHECK TO SEE IF THE INPUT VECTOR IS ALREADY SORTED
!
      NM1=N-1
      DO200I=1,NM1
        IP1=I+1
        IF(Y(I).LE.Y(IP1))GOTO200
        GOTO250
  200 CONTINUE
      RETURN
  250 M=1
      I=1
      J=N
  305 IF(I.GE.J)GOTO370
  310 K=I
      MID=(I+J)/2
      AMED=Y(MID)
      IF(Y(I).LE.AMED)GOTO320
      Y(MID)=Y(I)
      Y(I)=AMED
      AMED=Y(MID)
  320 L=J
      IF(Y(J).GE.AMED)GOTO340
      Y(MID)=Y(J)
      Y(J)=AMED
      AMED=Y(MID)
      IF(Y(I).LE.AMED)GOTO340
      Y(MID)=Y(I)
      Y(I)=AMED
      AMED=Y(MID)
      GOTO340
  330 Y(L)=Y(K)
      Y(K)=TT
  340 L=L-1
      IF(Y(L).GT.AMED)GOTO340
      TT=Y(L)
  350 K=K+1
      IF(Y(K).LT.AMED)GOTO350
      IF(K.LE.L)GOTO330
      LMI=L-I
      JMK=J-K
      IF(LMI.LE.JMK)GOTO360
      IL(M)=I
      IU(M)=L
      I=K
      M=M+1
      GOTO380
  360 IL(M)=K
      IU(M)=J
      J=L
      M=M+1
      GOTO380
  370 M=M-1
      IF(M.EQ.0)RETURN
      I=IL(M)
      J=IU(M)
  380 JMI=J-I
      IF(JMI.GE.11)GOTO310
      IF(I.EQ.1)GOTO305
      I=I-1
  390 I=I+1
      IF(I.EQ.J)GOTO370
      AMED=Y(I+1)
      IF(Y(I).LE.AMED)GOTO390
      K=I
  395 Y(K+1)=Y(K)
      K=K-1
      IF(AMED.LT.Y(K))GOTO395
      Y(K+1)=AMED
      GOTO390
      END






