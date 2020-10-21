c-------------------------------------------------------------------------
c
c	Subroutine calcmisfit - calculates a misfit value between
c				observed data and predicted data
c				plus model roughness
c
c	Note: Calls no other routines
c
c------------------------------------------------------------------------
c
      subroutine calcmisfit( predicted_data, observed_data, ndata,
     &       weight, nwaves, totmisfit, fname, amp_obs_real, 
     &       amp_obs_imag, write_best )
c
c
         include 'rfi_param.inc'
c
c
         complex*8	observed_data(maxdata,maxwave),
     &              predicted_data(maxdata,maxwave),
     &              observed_phase,
     &              predicted_phase

         real*4		weight(maxwave),
     &              abs_value_obs, abs_value_synt
         real*4		real_mft, imag_mft
         real*4		totmisfit, misfitval, aval_real, aval_imag, aval, 
     &              mft_abs, mft_fft
         real*4		mft_phase, mft_absolute, amplitude
         real*4		totmisfit_P,totmisfit_S,totmisfit_W
         integer		n_P, n_S, n_W
c
         integer         ndata(maxwave), nwaves
         character*40	fname(maxwave)
         character*40	instring, string_r, string1, string2, delim
         character	wavetype

         real*4		amp_obs_real(maxwave),
     &              amp_obs_imag(maxwave)
         real*4	 	a(15)
         real*4		m
         real*4		misfit_P_array(maxwave),misfit_S_array(maxwave),
     &              misfit_W_array(maxwave)
         real*4		median_P, median_S, median_W
         logical		write_best

         pi =4.D0*DATAN(1.D0)
c
c						misfit between observed and
c						predicted
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
c			----------------------------
c			Original code, do not edit
               real_mft = (real(predicted_data(i,iw)) - real(observed_data(i,iw)))**2
               imag_mft = (imag(predicted_data(i,iw)) - imag(observed_data(i,iw)))**2
c			----------------------------
               aval = aval + (real_mft + imag_mft)

            end do

            aval = weight(iw) * aval / (amp_obs_real(iw) + amp_obs_imag(iw))
            write(123,*) "STATION: ", fname(iw)(15:26) , " Misfit: ", aval

            if (write_best) then
               write(111,*) "STATION: ", fname(iw)(15:24) , " Misfit: ", aval
            end if

c		*** ORGINAL LINE ***
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

         write(*,*) "median P: ", median_P, " median S: ", median_S, " median W: ", median_W
         write(*,*) "n P: ", n_P, " n S: ", n_S, " n W: ", n_W

         return
      end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine reverse_string(instring, outstring)
         character*40	instring, outstring
         character	temp
         integer		i, length


         length = len_trim(instring) ! ignores trailing blanks.
         ! use len(string) to reverse those as well

         outstring = instring
         do i = 1, length/2
            temp = outstring(i:i)
            outstring(i:i) = outstring(length+1-i:length+1-i)
            outstring(length+1-i:length+1-i) = temp
         end do

      end subroutine reverse_string

      
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
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


!=========================================================================================
      SUBROUTINE MEDIAN(X,N,IWRITE,XMED)
C
C     PURPOSE--THIS SUBROUTINE COMPUTES THE
C              SAMPLE MEDIAN
C              OF THE DATA IN THE INPUT VECTOR X.
C              THE SAMPLE MEDIAN = THAT VALUE SUCH THAT HALF THE
C              DATA SET IS BELOW IT AND HALF ABOVE IT.
C     INPUT  ARGUMENTS--X      = THE SINGLE PRECISION VECTOR OF
C                                (UNSORTED OR SORTED) OBSERVATIONS.
C                     --N      = THE INTEGER NUMBER OF OBSERVATIONS
C                                IN THE VECTOR X.
C                     --IWRITE = AN INTEGER FLAG CODE WHICH
C                                (IF SET TO 0) WILL SUPPRESS
C                                THE PRINTING OF THE
C                                SAMPLE MEDIAN
C                                AS IT IS COMPUTED;
C                                OR (IF SET TO SOME INTEGER
C                                VALUE NOT EQUAL TO 0),
C                                LIKE, SAY, 1) WILL CAUSE
C                                THE PRINTING OF THE
C                                SAMPLE MEDIAN
C                                AT THE TIME IT IS COMPUTED.
C     OUTPUT ARGUMENTS--XMED   = THE SINGLE PRECISION VALUE OF THE
C                                COMPUTED SAMPLE MEDIAN.
C     OUTPUT--THE COMPUTED SINGLE PRECISION VALUE OF THE
C             SAMPLE MEDIAN.
C     PRINTING--NONE, UNLESS IWRITE HAS BEEN SET TO A NON-ZERO
C               INTEGER, OR UNLESS AN INPUT ARGUMENT ERROR
C               CONDITION EXISTS.
C     RESTRICTIONS--THE MAXIMUM ALLOWABLE VALUE OF N
C                   FOR THIS SUBROUTINE IS 15000.
C     OTHER DATAPAC   SUBROUTINES NEEDED--SORT.
C     FORTRAN LIBRARY SUBROUTINES NEEDED--NONE.
C     MODE OF INTERNAL OPERATIONS--SINGLE PRECISION.
C     LANGUAGE--ANSI FORTRAN.
C     REFERENCES--KENDALL AND STUART, THE ADVANCED THEORY OF
C                 STATISTICS, VOLUME 1, EDITION 2, 1963, PAGE 326.
C               --KENDALL AND STUART, THE ADVANCED THEORY OF
C                 STATISTICS, VOLUME 2, EDITION 1, 1961, PAGE 49.
C               --DAVID, ORDER STATISTICS, 1970, PAGE 139.
C               --SNEDECOR AND COCHRAN, STATISTICAL METHODS,
C                 EDITION 6, 1967, PAGE 123.
C               --DIXON AND MASSEY, INTRODUCTION TO STATISTICAL
C                 ANALYSIS, EDITION 2, 1957, PAGE 70.
C     WRITTEN BY--JAMES J. FILLIBEN
C                 STATISTICAL ENGINEERING LABORATORY (205.03)
C                 NATIONAL BUREAU OF STANDARDS
C                 WASHINGTON, D. C. 20234
C                 PHONE:  301-921-2315
C     ORIGINAL VERSION--JUNE      1972.
C     UPDATED         --SEPTEMBER 1975.
C     UPDATED         --NOVEMBER  1975.
C     UPDATED         --FEBRUARY  1976.
C
C---------------------------------------------------------------------
C
         DIMENSION X(1)
         DIMENSION Y(15000)
         COMMON /BLOCK2/ WS(15000)
         EQUIVALENCE (Y(1),WS(1))
C
         IPR=6
         IUPPER=15000
C
C     CHECK THE INPUT ARGUMENTS FOR ERRORS
C
         IF(N.LT.1.OR.N.GT.IUPPER)GOTO50
         IF(N.EQ.1)GOTO55
         HOLD=X(1)
         DO60I=2,N
            IF(X(I).NE.HOLD)GOTO90
  60     CONTINUE
         WRITE(IPR, 9)HOLD
         XMED=X(1)
         GOTO101
  50     WRITE(IPR,17)IUPPER
         WRITE(IPR,47)N
         RETURN
  55     WRITE(IPR,18)
         XMED=X(1)
         GOTO101
   90    CONTINUE
    9    FORMAT(1H ,109H***** NON-FATAL DIAGNOSTIC--THE FIRST  INPUT ARGUMENT (A VECTOR) TO THE MEDIAN SUBROUTINE HAS ALL ELEMENTS = ,E15.8,6H *****)
   17    FORMAT(1H , 98H***** FATAL ERROR--THE SECOND INPUT ARGUMENT TO THE MEDIAN SUBROUTINE IS OUTSIDE THE ALLOWABLE (1,,I6,16H) INTERVAL *****)
   18    FORMAT(1H ,100H***** NON-FATAL DIAGNOSTIC--THE SECOND INPUT ARGUMENT TO THE MEDIAN SUBROUTINE HAS THE VALUE 1 *****)
   47    FORMAT(1H , 35H***** THE VALUE OF THE ARGUMENT IS ,I8   ,6H *****)
C
C-----START POINT-----------------------------------------------------
C
         CALL SORT(X,N,Y)
         IFLAG=N-(N/2)*2
         NMID=N/2
         NMIDP1=NMID+1
         IF(IFLAG.EQ.0)XMED=(Y(NMID)+Y(NMIDP1))/2.0
         IF(IFLAG.EQ.1)XMED=Y(NMIDP1)
C
  101    IF(IWRITE.EQ.0)RETURN
         WRITE(IPR,999)
         WRITE(IPR,105)N,XMED
  105    FORMAT(1H ,25HTHE SAMPLE MEDIAN OF THE ,I6,17H OBSERVATIONS IS ,E5.8)
  999    FORMAT(1H )
         RETURN
      END



!===================================================================00
      SUBROUTINE SORT(X,N,Y)
C
C     PURPOSE--THIS SUBROUTINE SORTS (IN ASCENDING ORDER)
C              THE N ELEMENTS OF THE SINGLE PRECISION VECTOR X
C              AND PUTS THE RESULTING N SORTED VALUES INTO THE
C              SINGLE PRECISION VECTOR Y.
C     INPUT  ARGUMENTS--X      = THE SINGLE PRECISION VECTOR OF
C                                OBSERVATIONS TO BE SORTED.
C                     --N      = THE INTEGER NUMBER OF OBSERVATIONS
C                                IN THE VECTOR X.
C     OUTPUT ARGUMENTS--Y      = THE SINGLE PRECISION VECTOR
C                                INTO WHICH THE SORTED DATA VALUES
C                                FROM X WILL BE PLACED.
C     OUTPUT--THE SINGLE PRECISION VECTOR Y
C             CONTAINING THE SORTED
C             (IN ASCENDING ORDER) VALUES
C             OF THE SINGLE PRECISION VECTOR X.
C     PRINTING--NONE UNLESS AN INPUT ARGUMENT ERROR CONDITION EXISTS.
C     RESTRICTIONS--THE DIMENSIONS OF THE VECTORS IL AND IU
C                   (DEFINED AND USED INTERNALLY WITHIN
C                   THIS SUBROUTINE) DICTATE THE MAXIMUM
C                   ALLOWABLE VALUE OF N FOR THIS SUBROUTINE.
C                   IF IL AND IU EACH HAVE DIMENSION K,
C                   THEN N MAY NOT EXCEED 2**(K+1) - 1.
C                   FOR THIS SUBROUTINE AS WRITTEN, THE DIMENSIONS
C                   OF IL AND IU HAVE BEEN SET TO 36,
C                   THUS THE MAXIMUM ALLOWABLE VALUE OF N IS
C                   APPROXIMATELY 137 BILLION.
C                   SINCE THIS EXCEEDS THE MAXIMUM ALLOWABLE
C                   VALUE FOR AN INTEGER VARIABLE IN MANY COMPUTERS,
C                   AND SINCE A SORT OF 137 BILLION ELEMENTS
C                   IS PRESENTLY IMPRACTICAL AND UNLIKELY,
C                   THEN THERE IS NO PRACTICAL RESTRICTION
C                   ON THE MAXIMUM VALUE OF N FOR THIS SUBROUTINE.
C                   (IN LIGHT OF THE ABOVE, NO CHECK OF THE
C                   UPPER LIMIT OF N HAS BEEN INCORPORATED
C                   INTO THIS SUBROUTINE.)
C     OTHER DATAPAC   SUBROUTINES NEEDED--NONE.
C     FORTRAN LIBRARY SUBROUTINES NEEDED--NONE.
C     MODE OF INTERNAL OPERATIONS--SINGLE PRECISION.
C     LANGUAGE--ANSI FORTRAN.
C     COMMENT--THE SMALLEST ELEMENT OF THE VECTOR X
C              WILL BE PLACED IN THE FIRST POSITION
C              OF THE VECTOR Y,
C              THE SECOND SMALLEST ELEMENT IN THE VECTOR X
C              WILL BE PLACED IN THE SECOND POSITION
C              OF THE VECTOR Y, ETC.
C     COMMENT--THE INPUT VECTOR X REMAINS UNALTERED.
C     COMMENT--IF THE ANALYST DESIRES A SORT 'IN PLACE',
C              THIS IS DONE BY HAVING THE SAME
C              OUTPUT VECTOR AS INPUT VECTOR IN THE CALLING SEQUENCE.
C              THUS, FOR EXAMPLE, THE CALLING SEQUENCE
C              CALL SORT(X,N,X)
C              IS ALLOWABLE AND WILL RESULT IN
C              THE DESIRED 'IN-PLACE' SORT.
C     COMMENT--THE SORTING ALGORTHM USED HEREIN
C              IS THE BINARY SORT.
C              THIS ALGORTHIM IS EXTREMELY FAST AS THE
C              FOLLOWING TIME TRIALS INDICATE.
C              THESE TIME TRIALS WERE CARRIED OUT ON THE
C              UNIVAC 1108 EXEC 8 SYSTEM AT NBS
C              IN AUGUST OF 1974.
C              BY WAY OF COMPARISON, THE TIME TRIAL VALUES
C              FOR THE EASY-TO-PROGRAM BUT EXTREMELY
C              INEFFICIENT BUBBLE SORT ALGORITHM HAVE
C              ALSO BEEN INCLUDED--
C              NUMBER OF RANDOM        BINARY SORT       BUBBLE SORT
C               NUMBERS SORTED
C                N = 10                 .002 SEC          .002 SEC
C                N = 100                .011 SEC          .045 SEC
C                N = 1000               .141 SEC         4.332 SEC
C                N = 3000               .476 SEC        37.683 SEC
C                N = 10000             1.887 SEC      NOT COMPUTED
C     REFERENCES--CACM MARCH 1969, PAGE 186 (BINARY SORT ALGORITHM
C                 BY RICHARD C. SINGLETON).
C               --CACM JANUARY 1970, PAGE 54.
C               --CACM OCTOBER 1970, PAGE 624.
C               --JACM JANUARY 1961, PAGE 41.
C     WRITTEN BY--JAMES J. FILLIBEN
C                 STATISTICAL ENGINEERING LABORATORY (205.03)
C                 NATIONAL BUREAU OF STANDARDS
C                 WASHINGTON, D. C. 20234
C                 PHONE--301-921-2315
C     ORIGINAL VERSION--JUNE      1972.
C     UPDATED         --NOVEMBER  1975.
C
C---------------------------------------------------------------------
C
         DIMENSION X(1),Y(1)
         DIMENSION IU(36),IL(36)
C
         IPR=6
C
C     CHECK THE INPUT ARGUMENTS FOR ERRORS
C
         IF(N.LT.1)GOTO50
         IF(N.EQ.1)GOTO55
         HOLD=X(1)
         DO60I=2,N
            IF(X(I).NE.HOLD)GOTO90
   60    CONTINUE
         WRITE(IPR, 9)HOLD
         DO61I=1,N
            Y(I)=X(I)
   61    CONTINUE
         RETURN
   50    WRITE(IPR,15)
         WRITE(IPR,47)N
         RETURN
   55    WRITE(IPR,18)
         Y(1)=X(1)
         RETURN
   90    CONTINUE
    9    FORMAT(1H ,108H***** NON-FATAL DIAGNOSTIC--THE FIRST  INPUT ARGUMENT (A VECTOR) TO THE SORT   SUBROUTINE HAS ALL ELEMENTS = ,E15.8,6H *****)
   15    FORMAT(1H , 91H***** FATAL ERROR--THE SECOND INPUT ARGUMENT TO THE SORT   SUBROUTINE IS NON-POSITIVE *****)
   18    FORMAT(1H ,100H***** NON-FATAL DIAGNOSTIC--THE SECOND INPUT ARGUMENT TO THE SORT   SUBROUTINE HAS THE VALUE 1 *****)
   47    FORMAT(1H , 35H***** THE VALUE OF THE ARGUMENT IS ,I8   ,6H *****)
C
C-----START POINT-----------------------------------------------------
C
C     COPY THE VECTOR X INTO THE VECTOR Y
         DO100I=1,N
            Y(I)=X(I)
  100    CONTINUE
C
C     CHECK TO SEE IF THE INPUT VECTOR IS ALREADY SORTED
C
         NM1=N-1
         DO200I=1,NM1
            IP1=I+1
            IF(Y(I).LE.Y(IP1))GOTO200
            GOTO250
  200    CONTINUE
         RETURN
  250    M=1
         I=1
         J=N
  305    IF(I.GE.J)GOTO370
  310    K=I
         MID=(I+J)/2
         AMED=Y(MID)
         IF(Y(I).LE.AMED)GOTO320
         Y(MID)=Y(I)
         Y(I)=AMED
         AMED=Y(MID)
  320    L=J
         IF(Y(J).GE.AMED)GOTO340
         Y(MID)=Y(J)
         Y(J)=AMED
         AMED=Y(MID)
         IF(Y(I).LE.AMED)GOTO340
         Y(MID)=Y(I)
         Y(I)=AMED
         AMED=Y(MID)
         GOTO340
  330    Y(L)=Y(K)
         Y(K)=TT
  340    L=L-1
         IF(Y(L).GT.AMED)GOTO340
         TT=Y(L)
  350    K=K+1
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
  360    IL(M)=K
         IU(M)=J
         J=L
         M=M+1
         GOTO380
  370    M=M-1
         IF(M.EQ.0)RETURN
         I=IL(M)
         J=IU(M)
  380    JMI=J-I
         IF(JMI.GE.11)GOTO310
         IF(I.EQ.1)GOTO305
         I=I-1
  390    I=I+1
         IF(I.EQ.J)GOTO370
         AMED=Y(I+1)
         IF(Y(I).LE.AMED)GOTO390
         K=I
  395    Y(K+1)=Y(K)
         K=K-1
         IF(AMED.LT.Y(K))GOTO395
         Y(K+1)=AMED
         GOTO390
      END






