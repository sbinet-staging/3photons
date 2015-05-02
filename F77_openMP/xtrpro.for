C
C   Common pour passer les probabilites supplementaires
C
      REAL*4 PRPLUS,PRMOINS
      double precision EE1,EE2
      COMMON /xtrpro/ PRPLUS,PRMOINS,EE1,EE2

C present dans mc.f (naturellement)
C present dans book.f
      save /xtrpro/
!$OMP THREADPRIVATE(/xtrpro/)
