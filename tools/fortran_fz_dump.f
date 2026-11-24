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
      INTEGER IUHEAD(100)
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
         NUH = 100
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
