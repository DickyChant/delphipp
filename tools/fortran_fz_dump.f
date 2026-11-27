C SPDX-FileCopyrightText: 2026 DELPHI++ Contributors
C SPDX-License-Identifier: GPL-3.0-or-later
C
C Standalone Fortran 77 / CERNLIB ZEBRA FZ reference dumper.
C Reads an uncompressed DELPHI .fadsim / .sdst FZ file via CFOPEN+FZFILE('DL')
C with LREC=3840 (DELPHI default) and automatic fallback to LREC=900.
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
      INTEGER L, J, LL, NRUN, IEVT, NARGS, NTRY, FLEN, HANDL
      INTEGER LPV, LPA, NMOD, IP1, NPA, NCHG, NNEU, LRAW, LPMB, LRD
      INTEGER NRAWD, NPMB
      REAL PX, PY, PZ, EN, CH, PT, PMAG, SUMP, SUME
      INTEGER IUHEAD(4000)
      CHARACTER*256 FNAME
      CHARACTER*32 CLREC
      CHARACTER*4 CHIDH, CHSUB

      NARGS = IARGC()
      IF (NARGS .LT. 1) THEN
         WRITE(*,'(A)') 'Usage: fortran_fz_dump <file> [lrec]'
         STOP 2
      ENDIF
      CALL GETARG(1, FNAME)
      FLEN = LEN_TRIM(FNAME)
      LREC = 3840
      IF (NARGS .GE. 2) THEN
         CALL GETARG(2, CLREC)
         READ(CLREC, *) LREC
      ENDIF

      CALL MZEBRA(-1)
      CALL MZSTOR(IXSTOR, '//', ' ', FENCE, LQ, LQ, LQ,
     +            LQ(2000), LQ(NW))
      CALL MZLINK(IXSTOR, '/FZLNK/', LSUP, LSUP, LSUP)

      WRITE(*,'(A,A)') '=== FORTRAN CERNLIB FZ DUMP: ',
     +     FNAME(1:FLEN)

      NTRY = 0
 5    CONTINUE
      NTRY = NTRY + 1
      LUN = 11
      IQUEST(1) = 0
      CALL CFOPEN(HANDL, 0, LREC, 'r ', 0, FNAME(1:FLEN), ISTAT)
      IF (ISTAT .NE. 0) THEN
         WRITE(*,'(A,I6,1X,A)') 'CFOPEN failed ISTAT=', ISTAT,
     +        FNAME(1:FLEN)
         STOP 1
      ENDIF
      IQUEST(1) = HANDL
      CALL FZFILE(LUN, LREC, 'DL')

      NEVT = 0
 10   CONTINUE
         NUH = 4000
         LSUP = 0
         CALL FZIN(LUN, 0, LSUP, 1, ' ', NUH, IUHEAD)
         IF (IQUEST(1) .LT. 0) THEN
            WRITE(*,'(A,I6)') 'FZIN read error IQUEST(1)=', IQUEST(1)
            GOTO 90
         ELSEIF (IQUEST(1) .EQ. 1) THEN
            WRITE(*,'(A,I8)') '  [Start-of-Run] hdr=', IUHEAD(1)
            GOTO 10
         ELSEIF (IQUEST(1) .EQ. 2) THEN
            WRITE(*,'(A)') '  [End-of-Run]'
            GOTO 10
         ELSEIF (IQUEST(1) .GE. 3) THEN
            GOTO 90
         ENDIF

         IF (LSUP .EQ. 0) THEN
            WRITE(*,'(A,I5)') '  [Pilot-only record] nuh=', NUH
            GOTO 10
         ENDIF

         NEVT = NEVT + 1
         NRUN = 0
         IEVT = 0
         NBANK = 0
         CALL UHTOC(IQ(LSUP - 4), 4, CHIDH, 4)
         ND = IQ(LSUP - 1)
         IF (NUH .GE. 9) THEN
            NRUN = IUHEAD(6)
            IEVT = IUHEAD(9)
         ELSEIF (CHIDH .EQ. 'HEAD' .AND. ND .GE. 5) THEN
            NRUN = IQ(LSUP + 4)
            IEVT = IQ(LSUP + 5)
         ELSEIF (NUH .GE. 2) THEN
            NRUN = IUHEAD(1)
            IEVT = IUHEAD(2)
         ENDIF

         L = LSUP
 20      IF (L .NE. 0) THEN
            NBANK = NBANK + 1
            NL = IQ(L - 3)
            NS = IQ(L - 2)
            DO 30 J = 1, NS
               LL = LQ(L - J)
               IF (LL .NE. 0) THEN
                  NBANK = NBANK + 1
                  CALL UHTOC(IQ(LL - 4), 4, CHSUB, 4)
                  IF (NEVT .LE. 3 .AND. J .LE. 12) THEN
                     WRITE(*,'(A,I2,A,A4,A,I7,A,I4,A,I4)')
     +                    '    link[', J, ']=', CHSUB,
     +                    ' nd=', IQ(LL - 1),
     +                    ' nl=', IQ(LL - 3),
     +                    ' ns=', IQ(LL - 2)
                  ENDIF
               ENDIF
 30         CONTINUE
            L = LQ(L)
            GOTO 20
         ENDIF

         WRITE(*,'(A,I5,A,I8,A,I8,A,A4,A,I7,A,I6)')
     +        '  EVT #', NEVT, ' run=', NRUN, ' evt=', IEVT,
     +        ' top=', CHIDH, ' nd=', ND, ' banks=', NBANK

         IF (CHIDH .EQ. 'DTOP') THEN
            NPA = 0
            NCHG = 0
            NNEU = 0
            SUMP = 0.0
            SUME = 0.0
            LPV = LQ(LSUP - 2)
 35         IF (LPV .GT. 0) THEN
               LPA = LQ(LPV - 1)
 36            IF (LPA .GT. 0) THEN
                  IF (IQ(LPA - 1) .GE. 20) THEN
                     NMOD = IQ(LPA + 4)
                     IF (NMOD .GE. 1 .AND. NMOD .LE. 32) THEN
                        IP1 = 4 + NMOD
                        IF (IP1 + 11 .LE. IQ(LPA - 1)) THEN
                           PX = Q(LPA + IP1 + 3)
                           PY = Q(LPA + IP1 + 4)
                           PZ = Q(LPA + IP1 + 5)
                           EN = Q(LPA + IP1 + 6)
                           CH = Q(LPA + IP1 + 8)
                           PT = SQRT(PX * PX + PY * PY)
                           PMAG = SQRT(PT * PT + PZ * PZ)
                           IF (EN .GE. 0.0 .AND. EN .LE. 200.0) THEN
                              IF (ABS(CH - 1.0) .LT. 0.1 .OR.
     +                            ABS(CH - 2.0) .LT. 0.1) THEN
                                 IF (PT .GE. 1.0E-3) THEN
                                    NPA = NPA + 1
                                    NCHG = NCHG + 1
                                    SUMP = SUMP + PMAG
                                    SUME = SUME + EN
                                 ENDIF
                              ELSEIF (PMAG .GT. 1.0E-4 .AND.
     +                                EN .GT. 0.0) THEN
                                 NPA = NPA + 1
                                 NNEU = NNEU + 1
                                 SUMP = SUMP + PMAG
                                 SUME = SUME + EN
                              ENDIF
                           ENDIF
                        ENDIF
                     ENDIF
                  ENDIF
                  LPA = LQ(LPA)
                  GOTO 36
               ENDIF
               LPV = LQ(LPV)
               GOTO 35
            ENDIF
            WRITE(*,'(A,I5,A,I5,A,I5,A,F10.3,A,F10.3)')
     +           '    BANK_STAT npa=', NPA, ' nchg=', NCHG,
     +           ' nneu=', NNEU, ' sum_e=', SUME, ' sum_p=', SUMP
         ELSEIF (CHIDH .EQ. 'TOP ') THEN
            LRAW = LQ(LSUP - 1)
            NRAWD = 0
            NPMB = 0
            IF (LRAW .GT. 0) THEN
               IF (IQ(LRAW - 2) .GE. 22) THEN
                  DO 40 J = 5, 22
                     LPMB = LQ(LRAW - J)
                     IF (LPMB .GT. 0) NPMB = NPMB + 1
 42                  IF (LPMB .GT. 0) THEN
                        IF (IQ(LPMB - 2) .GE. 5) THEN
                           LRD = LQ(LPMB - 5)
 44                        IF (LRD .GT. 0) THEN
                              IF (IQ(LRD - 1) .GE. 2) NRAWD = NRAWD + 1
                              LRD = LQ(LRD)
                              GOTO 44
                           ENDIF
                        ENDIF
                        LPMB = LQ(LPMB)
                        GOTO 42
                     ENDIF
 40               CONTINUE
               ENDIF
            ENDIF
            WRITE(*,'(A,I6,A,I4)') '    RAW_STAT nrawd=', NRAWD,
     +           ' npmb=', NPMB
         ENDIF

         CALL MZWIPE(0)
         GOTO 10

 90   CONTINUE
      CALL FZENDI(LUN, 'TX')
      IF (NEVT .EQ. 0 .AND. NTRY .EQ. 1 .AND. LREC .EQ. 3840) THEN
         WRITE(*,'(A)') '  Retrying with LREC=900...'
         LREC = 900
         GOTO 5
      ENDIF
      WRITE(*,'(A,I6,A,I5)') 'Total events read by CERNLIB FZIN: ',
     +     NEVT, ' LREC=', LREC
      IF (NEVT .LE. 0) THEN
         STOP 1
      ENDIF
      END
