with Text_IO;
with Precision;
with Ppp;
with Spinor;
with Scalar;
with Angle;
with Result;
with Resfin;
with Cutpar;
with Resfin;
with Param;
with Xtrpro;
with Hasa;
with Ada.Text_IO;
--with Ada.Sequential_IO;
--with Ada.Streams.Stream_IO;
with Ada.Calendar;
with GNAT.Calendar.Time_IO;
with CPUClock;
----------------------------------------------------------------------
procedure Mc is --program RIEMANN
----------------------------------------------------------------------
   use Text_IO;
   use Precision;
   use Specific_Elementary;
   use Ppp;
   use Spinor;
   use Scalar;
   use Angle;
   use Result;
   use Resfin;
   use Cutpar;
   use Resfin;
   use Param;
   use Xtrpro;
   package Real_IO   is new Text_IO.Float_IO (Real);
   use  Real_IO;
   package Entier_IO is new Text_IO.Integer_IO (Integer);
   use  Entier_IO;
   package Boolean_IO is new Text_IO.Enumeration_io (Boolean);
   use Boolean_IO;
   --http://www.irit.fr/PERSONNEL/Veronique.Gaildrat/PagesEnseignement/PagesConceptsProgrammation/Ada95/Ada95_Genericite.html
   use CPUClock;

   --use Ada.Streams.Stream_IO;
   --package Sequoia is new Ada.Text_IO;
   --package Sequoia is new Ada.Sequential_IO(Character);
   --use Sequoia;
   use Ada.Calendar;
   use GNAT.Calendar.Time_IO;

   ETOT, FLUX, PI, NORM, POIDS, WTEV, SIGMA, PREC, VARIANCE : Real;
   POIDS2, PROBA : Real;
   MFIN : T_VECTEUR_PARTICULE; --MFIN : array (1..100) of Real;
   ITOT, NTOT, NBIN : Integer;
   ALPHA, CONVERS, SIN2W, COS2W, FACTCOM, E2, ALPHAZ, E2Z, GZR : Real;
   BETAPLUS, BETAMOINS, BETAMIN : Real;
   SSP, SSM, INCSSP, INCSSM : Real;
   PAA, PAB, PBB : Real;
   CAA, CAB, CBB : Real;
   DZETA, PROPAG, ECARTPIC : Real;
   BREPEM : Real;
   PLOT : Boolean;
   ICYCLE : Integer := 0;
   A : String (1..80);
   UNL : constant := 10;
   UNE : constant := 20;
   Now : Time := Clock;
   NowSeconds : Integer := Integer (Seconds (Now));

   --real*4 etime, tarray(2) : ;
   TIMESTR : String (1..8) := Image (Now, "%T");
   TODAY   : String (1..9) := Image (Now, "%d-%b-%y");
   --TODAY   : String (1..9) := Image (Now, "%y-%b-%d");

   OCutpar : Cutpar.Cutpar;
   OParam  : Param.Param;
   OResfin : Resfin.Resfin;
   OXtrpro : Xtrpro.Xtrpro;

   O1Ppp    : Ppp.Ppp;
   O1Spinor : Spinor.Spinor;
   O1Scalar : Scalar.Scalar;
   O1Angle  : Angle.Angle;
   O1Result : Result.Result;

   --Fichier : Sequoia.File_Type;
   Fichier : File_Type;
   Last : Natural;
   --Fichier : Ada.Streams.Stream_IO.File_Type;
   --Last    : Stream_Element_Offset;
   TrialTime : CPUClock.CPUSecond;   -- CPU time for each trial
   TrialSysTime : CPUClock.CPUSecond;   -- CPU time for each trial
begin
   CPUClock.ResetCPUTime;
-- 02-Dec-07   21:02:24
--   Integer'Image (Year (Now));
--   Integer'Image (Month (Now));
--   Integer'Image (Day (Now));
--   Integer'Image (NowSeconds / 3600); -- heures
--   Integer'Image ((NowSeconds rem 3600) / 60); -- minutes
--   Integer'Image ((NowSeconds rem 3600) rem 60); -- secondes
-- date +'%y-%m-%d %H:%M:%S'
-- date +'%y-%b-%d %H:%M:%S'
-- date +'%y-%b-%d %T'
--package GNAT.Calendar.Time_IO is
--   function Image
--     (Date    : Ada.Calendar.Time;
--      Picture : Picture_String) return String;
--   --  Return Date as a string with format Picture. Raise Picture_Error if
--   --  picture string is wrong.

--   Put (Year (Now)); Put ("-");
--   Put (Month (Now)); Put ("-");
--   Put (Day (Now)); Put (" ");
--   Put (NowSeconds / 3600);  Put (":");-- heures
--   Put ((NowSeconds rem 3600) / 60);  Put (":");-- minutes
--   Put ((NowSeconds rem 3600) rem 60); -- secondes
--   Put_Line ("");
--   Put_Line (Image (Now, "%y-%b-%d %T"));
--   Put_Line ("");
   --   Lecture des parametres
   Open (Fichier, In_File, "valeurs", ""); --open (UNIT=UNL, FILE='valeurs', STATUS='OLD')
   Get (Fichier, ITOT);           Get_Line (Fichier, A, Last);
   Get (Fichier, ETOT);           Get_Line (Fichier, A, Last);
   Get (Fichier, oCutpar.ACUT);   Get_Line (Fichier, A, Last);
   Get (Fichier, oCutpar.BCUT);   Get_Line (Fichier, A, Last);
   Get (Fichier, oCutpar.EMIN);   Get_Line (Fichier, A, Last);
   Get (Fichier, oCutpar.SINCUT); Get_Line (Fichier, A, Last);
   Get (Fichier, ALPHA);          Get_Line (Fichier, A, Last);
   Get (Fichier, ALPHAZ);         Get_Line (Fichier, A, Last);
   Get (Fichier, CONVERS);        Get_Line (Fichier, A, Last);
   Get (Fichier, oParam.MZ0);     Get_Line (Fichier, A, Last);
   Get (Fichier, oParam.GZ0);     Get_Line (Fichier, A, Last);
   Get (Fichier, SIN2W);          Get_Line (Fichier, A, Last);
   Get (Fichier, BREPEM);         Get_Line (Fichier, A, Last);
   Get (Fichier, BETAPLUS);       Get_Line (Fichier, A, Last);
   Get (Fichier, BETAMOINS);      Get_Line (Fichier, A, Last);
   Get (Fichier, NBIN);           Get_Line (Fichier, A, Last);
   Get (Fichier, oParam.IMPR);    Get_Line (Fichier, A, Last);
   Get (Fichier, PLOT);           Get_Line (Fichier, A, Last);
   --Put (" oParam.IMPR "); Put_Line (Boolean'Image (oParam.IMPR));
   --Put (" PLOT        "); Put_Line (Boolean'Image (PLOT));
   --Put ("ITOT "); Put (ITOT);
   --Put ("ETOT "); Put (ETOT);

   --   Read (Fichier, ITOT, Last);           Read (Fichier, A, Last);
   --   Read (Fichier, ETOT, Last);           Read (Fichier, A, Last);
   --   Read (Fichier, oCutpar.ACUT, Last);   Read (Fichier, A, Last);
   --   Read (Fichier, oCutpar.BCUT, Last);   Read (Fichier, A, Last);
   --   Read (Fichier, oCutpar.EMIN, Last);   Read (Fichier, A, Last);
   --   Read (Fichier, oCutpar.SINCUT, Last); Read (Fichier, A, Last);
   --   Read (Fichier, ALPHA, Last);          Read (Fichier, A, Last);
   --   Read (Fichier, ALPHAZ, Last);         Read (Fichier, A, Last);
   --   Read (Fichier, CONVERS, Last);        Read (Fichier, A, Last);
   --   Read (Fichier, oParam.MZ0, Last);     Read (Fichier, A, Last);
   --   Read (Fichier, oParam.GZ0, Last);     Read (Fichier, A, Last);
   --   Read (Fichier, SIN2W, Last);          Read (Fichier, A, Last);
   --   Read (Fichier, BREPEM, Last);         Read (Fichier, A, Last);
   --   Read (Fichier, BETAPLUS, Last);       Read (Fichier, A, Last);
   --   Read (Fichier, BETAMOINS, Last);      Read (Fichier, A, Last);
   --   Read (Fichier, NBIN, Last);           Read (Fichier, A, Last);
   --   Read (Fichier, oParam.IMPR, Last);    Read (Fichier, A, Last);
   --   Read (Fichier, PLOT, Last);           Read (Fichier, A, Last);
   Close (Fichier);

   --ITOT           := 10_000;
   --ETOT           := 91.187;
   --oCutpar.ACUT   := 0.9;
   --oCutpar.BCUT   := 0.9396e0;
   --oCutpar.EMIN   := 4.559;
   --oCutpar.SINCUT := 0.0;
   --ALPHA          := 7.297353079644818e-3;
   --ALPHAZ         := 7.8125e-3;
   --CONVERS        := 0.38937966e9;
   --oParam.MZ0     := 91.187;
   --oParam.GZ0     := 2.490;
   --SIN2W          := 0.2319;
   --BREPEM         := 0.03367;
   --BETAPLUS       := 1.0;
   --BETAMOINS      := 1.0;
   --NBIN           := 200;
   --oParam.IMPR    := False;
   --PLOT           := False;

--   Initialisation de PAW
  --CALL INIBOOK(ETOT, NBIN)

--   Masses des particules dans l'etat final
  MFIN := (others => 0.0);
--   Calcul du facteur de flux (=1/2s pour 2 particules initiales sans masse)
  FLUX := 1.0/(2.0*ETOT**2);
--   Inclue la normalisation totale de l'espace des phases
  PI := 4.0 * Arctan (1.0);
  NORM := ((2.0*PI)**(4-3*INP))/Real (ITOT); --NORM := ((2.0*PI)**(4-3*INP))/FLOAT(ITOT);
--   Facteur commun non-moyenne sur les spins=1
--                   /(facteur de symetrie)
--                   *(facteur de conversion GeV^-2->pb)
--   Pour moyenner sur les spins, rajouter :
--                   /(nombre d'helicites incidentes)
  FACTCOM := 1.0 / 6.0 * CONVERS;
  E2      := 4.0 * PI  * ALPHA;
  E2Z     := 4.0 * PI  * ALPHAZ;
  COS2W   := 1.0 - SIN2W;
  GZR     := oParam.GZ0 / oParam.MZ0;
  --Put_Line (Real'Image (Hasa.Rn));

--   Couplages
  oParam.GA  := -Sqrt (E2)**3;
  oParam.GBP := -Sqrt (E2Z/(4.0 * COS2W * SIN2W))/oParam.MZ0**4;
  oParam.GBM := oParam.GBP;
--   Facteurs du a la somme sur les polarisations
  oParam.POLP  :=     - 2.0 * SIN2W;
  oParam.POLM  := 1.0 - 2.0 * SIN2W;
  oParam.POLP2 := oParam.POLP**2;
  oParam.POLM2 := oParam.POLM**2;
  PAA := 2.0;
  PAB := 1.0 - 4.0 * SIN2W;
  PBB := 1.0 - 4.0 * SIN2W + 8.0 * SIN2W**2;
--   Coefficient pour homogeneiser
  CAA := FACTCOM * PAA;
  CAB := FACTCOM * PAB / oParam.MZ0**2;
  CBB := FACTCOM * PBB / oParam.MZ0**4;
--   Poids d'un evenement d'espace des phases a 3 photons
  WTEV := PI**2/8.0*ETOT**2;
--   Passage en variable a-dimensionee...
  DZETA    := (ETOT/oParam.MZ0)**2;
  ECARTPIC := (DZETA-1.0) / GZR;
  PROPAG   := 1.0 / (ECARTPIC**2 + 1.0);
----   Initialisation des impulsions entrantes
  o1Ppp := oPpp (ETOT);

--   Initialisation des cumulants
  NTOT := 0;
  --for Sp in 1..2 loop
  --   for K in Result_Category loop
  --      oResfin.SPM2 (Sp, K) := 0.0;
  --      oResfin.VAR  (Sp, K) := 0.0;
  --   end loop;
  --end loop;
  oResfin.SPM2 := (others => (others => 0.0));
  oResfin.VAR  := (others => (others => 0.0));
  SIGMA    := 0.0;
  VARIANCE := 0.0;

--   Debut de la boucle d'integration
  for I in 1..ITOT loop

     RAMBO (INP, ETOT, MFIN, o1Ppp, WTEV);
     WTEV := WTEV*NORM;
--   Trie les photons sortants par energie
     TRI (o1Ppp);
--   Calcul des produits spinoriels, produits scalaires et angles
     o1Spinor := oSpinor (o1Ppp);
     o1Scalar := oScalar (o1Spinor);
     o1Angle  := oAngle  (o1Ppp, o1Scalar);
--   Calcule le poids total de l'evenement avec l'element de matrice
--   (coupures s'il y a lieu)
     if not CUT (o1Ppp, oCutpar, o1Angle) then
        o1Result := oResult (o1Spinor, oParam, ETOT);
        oResfin.SPM2DIF := (others => 0.0);
        for K in Result_Category loop
           --oResfin.SPM2DIF (K) := 0.0;
           for L1 in 1..2 loop
              for L2 in 1..2 loop
                 for L3 in 1..2 loop
                    oResfin.SPM2DIF (K) := oResfin.SPM2DIF (K) + o1Result.M2 (L1, L2, L3, K);
                 end loop;
              end loop;
           end loop;
           oResfin.SPM2 (1, K) := oResfin.SPM2 (1, K) + oResfin.SPM2DIF (K);
           oResfin.VAR  (1, K) := oResfin.VAR  (1, K) + oResfin.SPM2DIF (K)**2;
        end loop;
        PROBA :=
             CAA*oResfin.SPM2DIF (1)
             +CBB*(BETAPLUS**2*oResfin.SPM2DIF (2)
             +BETAMOINS**2*oResfin.SPM2DIF (3))
             /GZR**2*PROPAG
             +CAB * 2.0 * BETAPLUS * (ECARTPIC*oResfin.SPM2DIF (4)
             -oResfin.SPM2DIF (5))
             /GZR*PROPAG;
        POIDS    := PROBA * WTEV/4.0;
        SIGMA    := SIGMA + POIDS;
        VARIANCE := VARIANCE + POIDS**2;
--    Stocke dans un histogramme les parametres de l'evenement
        if PLOT then
           POIDS2:=CAA*oResfin.SPM2DIF (1) * WTEV / 4.0;
           oXtrpro.PRPLUS := Real (CBB * BETAPLUS**2  * oResfin.SPM2DIF (2)
                /GZR**2*PROPAG
                *WTEV/4.0);
           oXtrpro.PRMOINS:= Real (CBB * BETAMOINS**2 * oResfin.SPM2DIF (3)
                /GZR**2*PROPAG
                *WTEV/4.0);
           --BOOK(POIDS, POIDS2)
        end if;
        NTOT := NTOT + 1;
     else
        POIDS := 0.0;
     end if;

  end loop; --$OMP END PARALLEL DO

-- Fin de la boucle d'echantillonnage

-- Calcul des incertitudes relatives
  for K in Result_Category loop
     oResfin.VAR (1, K) := (oResfin.VAR (1, K)-oResfin.SPM2(1, K)**2 / Real (ITOT)) / Real (ITOT-1);
     oResfin.VAR (1, K) := Sqrt (oResfin.VAR (1, K) / Real (ITOT)) / abs (oResfin.SPM2 (1, K)
          / Real (ITOT));
  end loop;
-- Copie pour les spins opposes
  for K in Result_Category loop
     oResfin.SPM2 (2, K) := oResfin.SPM2 (1, K);
     oResfin.VAR (2, K)  := oResfin.VAR (1, K);
  end loop;
-- Polarisations
  for K in 2..3 loop
     oResfin.SPM2 (1, K) := oResfin.SPM2 (1, K) * oParam.POLM2;
     oResfin.SPM2 (2, K) := oResfin.SPM2 (2, K) * oParam.POLP2;
  end loop;
  for K in 4..5 loop
     oResfin.SPM2 (1, K) := oResfin.SPM2 (1, K) * oParam.POLM;
     oResfin.SPM2 (2, K) := oResfin.SPM2 (2, K) * oParam.POLP;
  end loop;
-- Coefficients physiques et Propagateur du Z0
  for Sp in 1..2 loop
     for K in Result_Category loop
        oResfin.SPM2 (Sp, K) := oResfin.SPM2 (Sp, K) * FACTCOM * FLUX * WTEV;
     end loop;
     oResfin.SPM2 (Sp, 1) := oResfin.SPM2 (Sp, 1);
     oResfin.SPM2 (Sp, 2) := oResfin.SPM2 (Sp, 2) / GZR**2 / oParam.MZ0**4*PROPAG;
     oResfin.SPM2 (Sp, 3) := oResfin.SPM2 (Sp, 3) / GZR**2 / oParam.MZ0**4*PROPAG;
     oResfin.SPM2 (Sp, 4) := oResfin.SPM2 (Sp, 4) / GZR / oParam.MZ0**2*PROPAG*ECARTPIC;
     oResfin.SPM2 (Sp, 5) := oResfin.SPM2 (Sp, 5) / GZR / oParam.MZ0**2*PROPAG;
  end loop;
  BETAMIN := Sqrt ((oResfin.SPM2 (1, 1) + oResfin.SPM2 (2, 1))/(oResfin.SPM2 (1, 2) + oResfin.SPM2 (2, 2)));
  SSP := (oResfin.SPM2 (1, 2) + oResfin.SPM2 (2, 2))/ Sqrt (oResfin.SPM2 (1, 1) + oResfin.SPM2 (2, 1))/2.0;
  SSM := (oResfin.SPM2 (1, 3) + oResfin.SPM2 (2, 3))/ Sqrt (oResfin.SPM2 (1, 1) + oResfin.SPM2 (2, 1))/2.0;
  INCSSP :=
       Sqrt ((oResfin.SPM2 (1, 2) * oResfin.VAR (1, 2))**2
       +(oResfin.SPM2 (2, 2) * oResfin.VAR (2, 2))**2)/abs (oResfin.SPM2 (1, 2) + oResfin.SPM2 (2, 2))
       +Sqrt ((oResfin.SPM2 (1, 1) * oResfin.VAR (1, 1))**2
       +(oResfin.SPM2 (2, 1) * oResfin.VAR (2, 1))**2)/abs (oResfin.SPM2 (1, 1) + oResfin.SPM2 (2, 1))/2.0;
  INCSSM :=
       Sqrt ((oResfin.SPM2 (1, 3) * oResfin.VAR (1, 3))**2
       +(oResfin.SPM2 (2, 3) * oResfin.VAR (2, 3))**2)/abs (oResfin.SPM2 (1, 3) + oResfin.SPM2 (2, 3))
       +Sqrt ((oResfin.SPM2 (1, 1) * oResfin.VAR (1, 1))**2
       +(oResfin.SPM2 (2, 1) * oResfin.VAR (2, 1))**2)/abs (oResfin.SPM2 (1, 1) + oResfin.SPM2 (2, 1))/2.0;

  VARIANCE := (VARIANCE-SIGMA**2/Real (ITOT)) / Real (ITOT-1);
  PREC := Sqrt (VARIANCE/Real (ITOT))/abs (SIGMA/Real (ITOT));
  SIGMA := SIGMA * FLUX;
--   Normalisation des histogrammes
--  call NORMA (FLUX*NBIN)
--  call NORMSUP (ETOT, 4./(oResfin.SPM2 (1, 1)+oResfin.SPM2 (2, 1)), &
--                    4./(oResfin.SPM2 (1, 2)+oResfin.SPM2 (2, 2)), &
--                    4./(oResfin.SPM2 (1, 3)+oResfin.SPM2 (2, 3)))
--   Stockage des resultats de PAW
--  call HROUT(0, ICYCLE, ' ')
--  call HREND('DON')
--   Stockage des resultats numeriques
--   DATE(TODAY)
-- TIME(TIMESTR)
  --Open (UNE, Out_File, "res.dat", ""); --open (UNIT=UNE, FILE='res.dat', STATUS='UNKNOWN')
  --Put_Line (""); TODAY, '   ', TIMESTR
  Put (" ");
  Put (TODAY);
  Put ("   ");
  Put (TIMESTR);
  Put_Line ("");
  Put_Line ("");
  Put (" Nombre d'evenements            : ");  Put_Line (Integer'Image (ITOT));
  Put (" ... apres coupure              : ");  Put_Line (Integer'Image (NTOT));
  Put (" energie dans le CdM      (GeV) : ");  Put_Line (Real'Image (ETOT));
  Put (" coupure / cos(photon,faisceau) : ");  Put_Line (Real'Image (oCutpar.ACUT));
  Put (" coupure / cos(photon,photon)   : ");  Put_Line (Real'Image (oCutpar.BCUT));
  Put (" coupure / sin(normale,faisceau): ");  Put_Line (Real'Image (oCutpar.SINCUT));
  Put (" coupure sur l'energie    (GeV) : ");  Put_Line (Real'Image (oCutpar.EMIN));
  Put (" 1/(constante de structure fine): ");  Put_Line (Real'Image (1.0/ALPHA));
  Put (" 1/(structure fine au pic)      : ");  Put_Line (Real'Image (1.0/ALPHAZ));
  Put (" facteur de conversion GeV-2/pb : ");  Put_Line (Real'Image (CONVERS));
--Put (" Volume d'espace des phases     : ");  Put_Line (Real'Image (VOLUME/NTOT));
  Put (" Masse du Z0              (GeV) : ");  Put_Line (Real'Image (oParam.MZ0));
  Put (" Largeur du Z0            (GeV) : ");  Put_Line (Real'Image (oParam.GZ0));
  Put (" Sinus^2 Theta Weinberg         : ");  Put_Line (Real'Image (SIN2W));
  Put (" Taux de branchement Z--->e+e-  : ");  Put_Line (Real'Image (BREPEM));
  Put (" Beta plus                      : ");  Put_Line (Real'Image (BETAPLUS));
  Put (" Beta moins                     : ");  Put_Line (Real'Image (BETAMOINS));
  Put_Line ("---------------------------------------------");
  Put (" Section Efficace          (pb) : ");  Put_Line (Real'Image (SIGMA));
  Put (" Ecart-Type                (pb) : ");  Put_Line (Real'Image (SIGMA*PREC));
  Put (" Precision Relative             : ");  Put_Line (Real'Image (PREC));
  Put_Line ("---------------------------------------------");
  Put (" Beta minimum                   : ");  Put_Line (Real'Image (BETAMIN));
  Put (" Stat. Significance  B+(pb-1/2) : ");  Put_Line (Real'Image (SSP));
  Put (" Incert. Stat. Sign. B+(pb-1/2) : ");  Put_Line (Real'Image (SSP*INCSSP));
  Put (" Stat. Significance  B-(pb-1/2) : ");  Put_Line (Real'Image (SSM));
  Put (" Incert. Stat. Sign. B+(pb-1/2) : ");  Put_Line (Real'Image (SSM*INCSSM));
  TrialTime := CPUClock.CPUTime;
  TrialSysTime := CPUClock.CPUSysTime;
  Put (" Temps ecoule                   : ");  Put_Line (CPUSecond'Image (TrialTime+TrialSysTime));
  Put (" Temps ecoule utilisateur       : ");  Put_Line (CPUSecond'Image (TrialTime));
  Put (" Temps ecoule systeme           : ");  Put_Line (CPUSecond'Image (TrialSysTime));
  Put (" Temps ecoule par evenement     : ");  Put_Line (Real'Image (Real (TrialTime) / Real (ITOT)));
--Put (" Temps ecoule                   : ");  Put_Line (Real'Image (Etime (Tarray)));
--Put (" Temps ecoule utilisateur       : ");  Put_Line (Real'Image (Tarray (1)));
--Put (" Temps ecoule systeme           : ");  Put_Line (Real'Image (Tarray (2)));
--Put (" Temps ecoule par evenement     : ");  Put_Line (Real'Image (Tarray (1) / Real (ITOT)));
  Put_Line ("");
  for Sp in 1..2 loop
     for K in Result_Category loop
        --write (UNE, 930) --930 format (2I3, 3E15.7)
        Put (Sp, 3);
        Put (K, 3);
        Put (oResfin.SPM2 (Sp, K), 3, 7, 3);
        Put (abs oResfin.SPM2 (Sp, K) * oResfin.VAR (Sp, K), 3, 7, 3);
        Put (oResfin.VAR (Sp, K), 3, 7, 3);
        Put_Line ("");
--        Put (Integer'Image (Sp)); Put (" ");
--        Put (Integer'Image (K)); Put (" ");
--        Put (Real'Image (oResfin.SPM2 (Sp, K))); Put (" ");
--        Put (Real'Image (abs oResfin.SPM2 (Sp, K) * oResfin.VAR (Sp, K))); Put (" ");
--        Put (Real'Image (oResfin.VAR (Sp, K)));
--        Put_Line ("");
     end loop;
     Put_Line ("");
  end loop;
  for K in Result_Category loop
     --write (UNE, 940) --940 format (3X, I3, 3E15.7)
     Put ("   ");
     Put (K, 3);
     Put ((oResfin.SPM2 (1, K)+oResfin.SPM2 (2, K))/4.0, 3, 7, 3);
     Put (Sqrt ((oResfin.SPM2 (1, K) * oResfin.VAR(1, K))**2
                           +(oResfin.SPM2 (2, K) * oResfin.VAR(2, K))**2)/4.0, 3, 7, 3);
     Put (Sqrt ((oResfin.SPM2 (1, K) * oResfin.VAR(1, K))**2
                           +(oResfin.SPM2 (2, K) * oResfin.VAR(2, K))**2)
                       /abs (oResfin.SPM2 (1, K)+oResfin.SPM2 (2, K)), 3, 7, 3);
--     Put (Integer'Image (K)); Put (" ");
--     Put (Real'Image ((oResfin.SPM2 (1, K)+oResfin.SPM2 (2, K))/4.0)); Put (" ");
--     Put (Real'Image (Sqrt ((oResfin.SPM2 (1, K) * oResfin.VAR(1, K))**2
--                           +(oResfin.SPM2 (2, K) * oResfin.VAR(2, K))**2)/4.0)); Put (" ");
--     Put (Real'Image (Sqrt ((oResfin.SPM2 (1, K) * oResfin.VAR(1, K))**2
--                           +(oResfin.SPM2 (2, K) * oResfin.VAR(2, K))**2)
--                       /abs (oResfin.SPM2 (1, K)+oResfin.SPM2 (2, K))));
     Put_Line ("");
  end loop;
  ERIC (UNE, BREPEM, CONVERS, PI, oParam, oResfin);
  FAWZI (UNE, BREPEM, CONVERS, PI, ETOT, oParam, oCutpar, oResfin);
  --close (UNIT=UNE)

  if False then
     --open (UNIT=UNE, ACCESS='APPEND', FILE='pil.mc', STATUS='unknown')
     Put (TODAY);
     Put ("   ");
     Put (TIMESTR);
     Put_Line ("");
     --write (UNE, 960)
     --960 format (G12.4, 6E15.7)
     Put (Real'Image (ETOT)); Put (" ");
     Put (Real'Image ((oResfin.SPM2 (1, 1)+oResfin.SPM2 (2, 1))/4.0)); Put (" ");
     Put (Real'Image ((oResfin.SPM2 (1, 2)+oResfin.SPM2 (2, 2))/4.0*BETAPLUS**2)); Put (" ");
     Put (Real'Image ((oResfin.SPM2 (1, 3)+oResfin.SPM2 (2, 3))/4.0*BETAMOINS**2)); Put (" ");
     Put (Real'Image ((oResfin.SPM2 (1, 4)+oResfin.SPM2 (2, 4))/4.0*BETAPLUS)); Put (" ");
     Put (Real'Image (((oResfin.SPM2 (1, 1)+oResfin.SPM2 (2, 1))
                      +(oResfin.SPM2 (1, 2)+oResfin.SPM2 (2, 2))*BETAPLUS**2
                      +(oResfin.SPM2 (1, 3)+oResfin.SPM2 (2, 3))*BETAMOINS**2
                  +2.0*(oResfin.SPM2 (1, 4)+oResfin.SPM2 (2, 4))*BETAPLUS)/4.0)); Put (" ");
     Put (Real'Image (SIGMA)); Put (" ");
     Put_Line ("");
     --close (UNIT=UNE)
  end if;

end Mc;
