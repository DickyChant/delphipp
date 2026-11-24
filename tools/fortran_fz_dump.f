C SPDX-FileCopyrightText: 2026 DELPHI++ Contributors
C SPDX-License-Identifier: GPL-3.0-or-later
C
C Standalone Fortran 77 / CERNLIB ZEBRA FZ reference dumper.
C Reads an uncompressed DELPHI .fadsim / .sdst FZ file via CFOPEN+FZFILE('DL')
C and prints event headers and top-level / daughter ZEBRA bank hierarchies
C for direct side-by-side comparison with delphipp_fz_inspect & delphipp_reco.
C
      PROGRAM FZDUMP
      IMPLICIT NONE
      INTEGER NW
      PARAMETER (NW = 500000)
      COMMON /QUEST/ IQUEST(100)
      INTEGER IQUEST
      COMMON // FENCE(16), LQ(NW)
      INTEGER FENCE, LQ
      INTEGER IQ(NW)
      REAL Q(NW)
      EQUIVALENCE (LQ(9), IQ(1), Q(1))

      COMMON /FZLNK/ LSUP
      INTEGER LSUP

      INTEGER IXSTOR, LUN, LREC, ISTAT, NUH, NEVT, NBANK, ND, NL, NS
      INTEGER L, J, LL, NRUN, IEVT, NARGS
      INTEGER IUHEAD(100)
      CHARACTER*256 FNAME
      CHARACTER*4 CHIDH, CHSUB

      NARGS = IARGC()
      IF (NARGS .LT. 1) THEN
         WRITE(*,'(A)') 'Usage: fortran_fz_dump <file.fadsim|file.sdst>'
         STOP 2
      ENDIF
      CALL GETARG(1, FNAME)

      CALL MZEBRA(-1)
      CALL MZSTOR(IXSTOR, '//', ' ', FENCE, LQ, LQ, LQ,
     +            LQ(2000), LQ(NW))
      CALL MZLINK(IXSTOR, '/FZLNK/', LSUP, LSUP, LSUP)

      LUN  = 11
      LREC = 900
      CALL CFOPEN(IQUEST(1), 0, LREC, 'r ', 0, FNAME, ISTAT)
      IF (ISTAT .NE. 0) THEN
         WRITE(*,'(A,I6,1X,A)') 'CFOPEN failed ISTAT=', ISTAT,
     +        FNAME(1:LEN_TRIM(FNAME))
         STOP 1
      ENDIF
      CALL FZFILE(LUN, LREC, 'DL')

      WRITE(*,'(A,A)') '=== FORTRAN CERNLIB FZ DUMP: ',
     +     FNAME(1:LEN_TRIM(FNAME))

      NEVT = 0
 10   CONTINUE
         NUH = 100
         LSUP = 0
         CALL FZIN(LUN, IXSTOR, LSUP, 1, ' ', NUH, IUHEAD)
         IF (IQUEST(1) .LT. 0) THEN
            WRITE(*,'(A,I6)') 'FZIN read error IQUEST(1)=', IQUEST(1)
            GOTO 90
         ELSEIF (IQUEST(1) .EQ. 1) THEN
            WRITE(*,'(A,I8)') '  [Start-of-Run] IQUEST(1)=1 hdr=',
     +           IUHEAD(1)
            GOTO 10
         ELSEIF (IQUEST(1) .EQ. 2) THEN
            WRITE(*,'(A)') '  [End-of-Run]'
            GOTO 10
         ELSEIF (IQUEST(1) .GE. 3) THEN
            GOTO 90
         ENDIF

         NEVT = NEVT + 1
         NRUN = 0
         IEVT = 0
         NBANK = 0
         IF (LSUP .NE. 0) THEN
            CALL UHTOC(IQ(LSUP - 4), 4, CHIDH, 4)
            ND = IQ(LSUP - 1)
            IF (CHIDH .EQ. 'HEAD' .AND. ND .GE. 5) THEN
               NRUN = IQ(LSUP + 4)
               IEVT = IQ(LSUP + 5)
            ELSEIF (NUH .GE. 2) THEN
               NRUN = IUHEAD(1)
               IEVT = IUHEAD(2)
            ENDIF
            WRITE(*,'(A,I4,A,I6,A,I6,A,A4,A,I5)')
     +           '  EVT #', NEVT, ' run=', NRUN, ' evt=', IEVT,
     +           ' top=', CHIDH, ' nd=', ND

            L = LSUP
 20         IF (L .NE. 0) THEN
               NBANK = NBANK + 1
               NL = IQ(L - 3)
               NS = IQ(L - 2)
               DO 30 J = 1, NS
                  LL = LQ(L - J)
                  IF (LL .NE. 0) THEN
                     NBANK = NBANK + 1
                     CALL UHTOC(IQ(LL - 4), 4, CHSUB, 4)
                     IF (NEVT .LE. 3) THEN
                        WRITE(*,'(A,I2,A,A4,A,I6,A,I4,A,I4)')
     +                       '    link[', J, ']=', CHSUB,
     +                       ' nd=', IQ(LL - 1), ' nl=', IQ(LL - 3),
     +                       ' ns=', IQ(LL - 2)
                     ENDIF
                  ENDIF
 30            CONTINUE
               L = LQ(L)
               GOTO 20
            ENDIF
         ENDIF
         CALL MZWIPE(IXSTOR)
         GOTO 10

 90   CONTINUE
      CALL FZENDI(LUN, 'TX')
      WRITE(*,'(A,I6)') 'Total events read by CERNLIB FZIN: ', NEVT
      IF (NEVT .LE. 0) THEN
         STOP 1
      ENDIF
      END
