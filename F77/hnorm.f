C----------------------------------------------------------------------
      SUBROUTINE NORMA(COEFF)
C----------------------------------------------------------------------
C
C Effectue la normalisations des histogrammes 
C
C----------------------------------------------------------------------
      IMPLICIT NONE
      REAL*8 COEFF
      REAL*4 WGTPAW,ZERO
      PARAMETER (ZERO=0.E0)
      INCLUDE 'paw.i'

      WGTPAW=SNGL(COEFF)

      CALL HOPERA(110,'+',110,110,WGTPAW,ZERO)
      CALL HOPERA(120,'+',120,120,WGTPAW,ZERO)
      CALL HOPERA(130,'+',130,130,WGTPAW,ZERO)
      CALL HOPERA(140,'+',140,140,WGTPAW,ZERO)
      CALL HOPERA(150,'+',150,150,WGTPAW,ZERO)
      CALL HOPERA(160,'+',160,160,WGTPAW,ZERO)
      CALL HOPERA(170,'+',170,170,WGTPAW,ZERO)
      CALL HOPERA(180,'+',180,180,WGTPAW,ZERO)
      CALL HOPERA(190,'+',190,190,WGTPAW,ZERO)

      CALL HOPERA(210,'+',210,210,WGTPAW,ZERO)
      CALL HOPERA(220,'+',220,220,WGTPAW,ZERO)
      CALL HOPERA(230,'+',230,230,WGTPAW,ZERO)
      CALL HOPERA(240,'+',240,240,WGTPAW,ZERO)
      CALL HOPERA(250,'+',250,250,WGTPAW,ZERO)
      CALL HOPERA(260,'+',260,260,WGTPAW,ZERO)
      CALL HOPERA(270,'+',270,270,WGTPAW,ZERO)
      CALL HOPERA(280,'+',280,280,WGTPAW,ZERO)
      CALL HOPERA(290,'+',290,290,WGTPAW,ZERO)

      RETURN
      END

C----------------------------------------------------------------------
