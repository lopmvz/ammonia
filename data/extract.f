      PROGRAM EXTRACT
      IMPLICIT REAL*8 (A-H,O-Z) 
      PARAMETER(MAXA=500)
      DIMENSION CELLDM(6),HT(3,3),HTM1(3,3)
     X         ,X(MAXA),Y(MAXA),Z(MAXA)
     X         ,X2(MAXA),Y2(MAXA),Z2(MAXA)
     X         ,INDEXO(MAXA),INDEXQM(MAXA)

      CHARACTER*80 CUMMY,BASIS1,BASIS2,POTO,POTH
      CHARACTER*3 CATOM(MAXA)
C
C
 1    READ(5,*,END=10)  NATOM
      IF(NATOM.GT.MAXA) STOP 'Too many atoms'
      READ(5,*) CUMMY,A,B,C,ALPHA,BETA,GAMMA
C
C     select size in angstrom of cluster
C
C      OODIST=DSQRT(A*A+B*B+C*C)/2.D0
      OODIST=12.0D0
      BONDDIST=1.3D0
C
      PI = DACOS(-1.0D0)
      TODEGREE = 180.D0/PI
      CELLDM(1) = A
      CELLDM(2) = B/CELLDM(1)
      CELLDM(3) = C/CELLDM(1)
      CELLDM(4) = DCOS(ALPHA/TODEGREE)
      CELLDM(5) = DCOS(BETA/TODEGREE)
      CELLDM(6) = DCOS(GAMMA/TODEGREE)
C
C MAKE HT AND HTM1
      SINGAM=DSQRT(1.D0-CELLDM(6)**2) 
      TERM=DSQRT((1.D0+2.D0*CELLDM(4)*CELLDM(5)*CELLDM(6) 
     X     -CELLDM(4)**2-CELLDM(5)**2-CELLDM(6)**2)/(1.D0-CELLDM(6)**2)) 
      HT(1,1)=CELLDM(1)
      HT(1,2)=0.0D0
      HT(1,3)=0.0D0
      HT(2,1)=CELLDM(1)*CELLDM(2)*CELLDM(6) 
      HT(2,2)=CELLDM(1)*CELLDM(2)*SINGAM 
      HT(2,3)=0.0D0
      HT(3,1)=CELLDM(1)*CELLDM(3)*CELLDM(5) 
      HT(3,2)=CELLDM(1)*CELLDM(3)
     X     *(CELLDM(4)-CELLDM(5)*CELLDM(6))/SINGAM 
      HT(3,3)=CELLDM(1)*CELLDM(3)*TERM 
C     XYZ -> ABC
      HTM1(1,1)= 1.D0/HT(1,1)
      HTM1(1,2)= 0.D0
      HTM1(1,3)= 0.D0
      HTM1(2,1)=-HT(2,1)/(HT(1,1)*HT(2,2))
      HTM1(2,2)= 1.D0/HT(2,2)
      HTM1(2,3)= 0.D0
      HTM1(3,1)= (-HT(2,2)*HT(3,1)+HT(2,1)*HT(3,2))
     X     /(HT(1,1)*HT(2,2)*HT(3,3))
      HTM1(3,2)=-HT(3,2)/(HT(2,2)*HT(3,3))
      HTM1(3,3)=1.D0/HT(3,3)
C
      DO I = 1,NATOM
         INDEXO(I)=0
         READ(5,*) CATOM(I),X(I),Y(I),Z(I)
!         WRITE(*,'(3A4)') "XXXX",CATOM(I),"XXXX"
      END DO
C
C Find atN
C
      DO I = 1,NATOM
         IF(CATOM(I).EQ.'N  ') MOL1=I
      END DO
C
C Printing the molecules closest to MOL1
C
      DDCRIT=OODIST*OODIST
      DO I = 1,NATOM
         IF(CATOM(I).NE.'H  ') THEN
!         WRITE(*,'(3A4)') "XXXX",CATOM(I),"XXXX"
            DX=X(I)-X(MOL1)
            DY=Y(I)-Y(MOL1)
            DZ=Z(I)-Z(MOL1)
            CALL DISTANCE(DX,DY,DZ,HT,HTM1)
            X2(I)=DX
            Y2(I)=DY
            Z2(I)=DZ
            DD=DX*DX+DY*DY+DZ*DZ
            IF(DD.LE.DDCRIT) THEN
               INDEXO(I)=2
               END IF
         END IF
      END DO
C     
      DDCRIT=BONDDIST*BONDDIST
      DO I = 1,NATOM
         DO J = 1,NATOM
            IF(INDEXO(I).EQ.2) THEN
               IF(CATOM(J).EQ.'H  ') THEN
                  DX=X(J)-X(I)
                  DY=Y(J)-Y(I)
                  DZ=Z(J)-Z(I)
                  CALL DISTANCE(DX,DY,DZ,HT,HTM1)
                  DD=DX*DX+DY*DY+DZ*DZ
                  IF(DD.LE.DDCRIT) THEN
                  X2(J)=DX+X2(I)
                  Y2(J)=DY+Y2(I)
                  Z2(J)=DZ+Z2(I)
!                     WRITE(*,*) CATOM(I),CATOM(J),I,J,DX,DY,DZ
                     INDEXO(J)=1
                  END IF
               END IF
            END IF
         END DO
      END DO
      
C     
C     
C     
      II=0
      ICHARGE=0
      DO I = 1,NATOM
         IF(INDEXO(I).NE.0) THEN
            II=II+1
            IF(CATOM(I).EQ.'N') ICHARGE=ICHARGE-3
            IF(CATOM(I).EQ.'O') ICHARGE=ICHARGE-2
            IF(CATOM(I).EQ.'H') ICHARGE=ICHARGE+1
            IF(CATOM(I).EQ.'Cl') ICHARGE=ICHARGE-1
            END IF
!            write(*,*) I,INDEXO(I),"X",CATOM(I),"X",'Cl ',"X",ICHARGE
      END DO
      WRITE(*,'(I6)')   II
      WRITE(*,'(A20,6F8.4)') CUMMY,A,B,C,ALPHA,BETA,GAMMA
      WRITE(10,'(2I5)') ICHARGE,1
C
      DO I = 1,NATOM
         IF(INDEXO(I).NE.0) WRITE(*,'(A9,3F16.6,A1,I6,A1,I5,A1,I1)') 
     &        CATOM(I),X2(I),Y2(I),Z2(I)
     &        ,' ',INDEXO(I),' ',I,' ', 1
!     add molecule #, atominmol, type_of_atom
      END DO
      
C 
C Select the number of solvent molecules NM2 to include in the QM
C region and the number of atoms per solute molecule NATM2
C
C

C Set the number of solvent molecules in the QM region:
      NM2=1 !Number of solvent molecules in the QM region

c      NATM1 !Number of atoms of a solute molecule
c      NATM2 !Number of atoms of a solvent molecule

C Set the basis set:
      BASIS1="6-311++G**" !Solute 
      BASIS2="6-311++G**" !Solvent  
C   
C Set the potential parameters:
      POTO="-0.67444   5.73935"
      POTH=" 0.33722   2.30839"



C MOLECULAR.INP

      DO I=1,NATOM
         INDEXQM(I)=0
      END DO

      DDCRIT=BONDDIST*BONDDIST
      NATM1=0
      DO I=1,NATOM
         IF (CATOM(I).EQ.'N') THEN
            INDEXQM(I)=1
            NATM1=NATM1+1
            DO J=1,NATOM
               IF (CATOM(J).EQ.'H') THEN
                  DD=X2(J)*X2(J)+Y2(J)*Y2(J)+Z2(J)*Z2(J)
                  IF (DD.LE.DDCRIT) THEN
                     INDEXQM(J)=1
                     NATM1=NATM1+1
                  END IF
               END IF
            END DO
         END IF
      END DO

      DO I=1,NM2
         DD0=1000
         NATM2=0
         DO J=1,NATOM
            IF (CATOM(J).EQ.'O' .AND. INDEXQM(J).EQ.0) THEN
               DD=X2(J)*X2(J)+Y2(J)*Y2(J)+Z2(J)*Z2(J)
               IF (DD.LE.DD0) THEN
                  DD0=DD
                  MOL2=J
               END IF
             END IF
          END DO
          INDEXQM(MOL2)=2
          NATM2=NATM2+1
          DO J=1,NATOM
             IF (CATOM(J).EQ.'H' .AND. INDEXQM(J).EQ.0) THEN
               DX=X2(J)-X2(MOL2)
               DY=Y2(J)-Y2(MOL2)
               DZ=Z2(J)-Z2(MOL2)
               DD=DX*DX+DY*DY+DZ*DZ
               IF (DD.LE.DDCRIT) THEN
                  INDEXQM(J)=2 
                  NATM2=NATM2+1
               END IF
             END IF
          END DO
      END DO


      OPEN (UNIT=1,FILE="qm.mol")
      WRITE(1,'(A9)') "ATOMBASIS"
      WRITE(1,*) ""
      WRITE(1,*) ""
      WRITE(1,'(A11,A8,I1,A7,A9)') "Atomtypes=4","Charge=",ICHARGE
     &                  ,"Nosymm", "Angstrom"


      WRITE(1,'(A10,I5,2A10)') "7.00000", 1,"Basis=",BASIS1
      DO I=1,NATOM
         IF (CATOM(I).EQ.'N' .AND. INDEXQM(I).EQ.1) THEN
            WRITE(1,'(A3,3F16.6)') CATOM(I),X2(I),Y2(I),Z2(I)
         END IF
      END DO

      WRITE(1,'(A10,I5,2A10)') "1.00000",NATM1-1,"Basis=", BASIS1
      DO I=1,NATOM
         IF (CATOM(I).EQ.'H' .AND. INDEXQM(I).EQ.1) THEN
            WRITE(1,'(A3,3F16.6)') CATOM(I),X2(I),Y2(I),Z2(I)
        END IF
      END DO

      WRITE(1,'(A10,I5,2A10)') "8.00000",NM2,"Basis=",BASIS2
      DO I=1,NATOM
         IF (CATOM(I).EQ.'O' .AND. INDEXQM(I).EQ.2) THEN
            WRITE(1,'(A3,3F16.6)') CATOM(I),X2(I),Y2(I),Z2(I)
        END IF
      END DO

      WRITE(1,'(A10,I5,2A10)') "1.00000",NM2*2,"Basis=",BASIS2
      DO I=1,NATOM
         IF (CATOM(I).EQ.'H' .AND. INDEXQM(I).EQ.2) THEN
            WRITE(1,'(A3,3F16.6)') CATOM(I),X2(I),Y2(I),Z2(I)
        END IF
      END DO

C POTENTIAL.INP
      OPEN(UNIT=2,FILE='mm.pot')     
      WRITE(2,*) "AA"
      WRITE(2,'(I3,3A3)') II-1*NATM1-NM2*NATM2, '0', '1', '1'
      NMM=0
      DO I=1,NATOM
         IF (CATOM(I).EQ.'O' .AND. INDEXQM(I).EQ.0 .AND.
     & INDEXO(I).NE.0) THEN
            NMM=NMM+1
c            WRITE(2,'(A3,3F16.6)') CATOM(I),X2(I),Y2(I),Z2(I)
            WRITE(2,'(I3,3F16.6,A3,A20)') NMM,X2(I),Y2(I),Z2(I)," ",POTO
            DO J=1,NATOM
               IF (CATOM(J).EQ.'H') THEN
                  DX=X2(J)-X2(I)
                  DY=Y2(J)-Y2(I)
                  DZ=Z2(J)-Z2(I)
                  DD=DX*DX+DY*DY+DZ*DZ
                  IF (DD.LE.DDCRIT) THEN
c                    WRITE(2,'(A3,3F16.6)') CATOM(J),X2(J),Y2(J),Z2(J)
                     WRITE(2,'(I3,3F16.6,A3,A20)') NMM,X2(J),Y2(J)
     &,Z2(J)," ",POTH
                  END IF
               END IF
            END DO
         END IF
      END DO
 

      GOTO 1
 10   END
C====================================================================
      SUBROUTINE DISTANCE(DX,DY,DZ,HT,HTM1)
      IMPLICIT REAL*8 (A-H,O-Z) 
      DIMENSION R(3),S(3),HT(3,3),HTM1(3,3)
C     
      R(1)=DX 
      R(2)=DY 
      R(3)=DZ 
C     XYZ ->ABC ; A = X * M-1
      DO I=1,3
         S(I) = 0.D0
         DO J=1,3
            S(I) = S(I) + R(J)*HTM1(J,I)
         END DO
      END DO
C     PBC
      S(1)=S(1)-NINT(S(1))
      S(2)=S(2)-NINT(S(2))
      S(3)=S(3)-NINT(S(3))
C     ABC->XYZ ; A = X * M
      DO I=1,3
         R(I) = 0.D0
         DO J=1,3
            R(I) = R(I) + S(J)*HT(J,I)
         END DO
      END DO
C     
      DX=R(1) 
      DY=R(2) 
      DZ=R(3) 
C     
      RETURN
      END
C====================================================================
      SUBROUTINE XYZTOABC(DX,DY,DZ,HT,HTM1)
      IMPLICIT REAL*8 (A-H,O-Z) 
      DIMENSION R(3),S(3),HT(3,3),HTM1(3,3)
C     
      R(1)=DX 
      R(2)=DY 
      R(3)=DZ 
C     XYZ ->ABC ; A = X * M-1
      DO I=1,3
         S(I) = 0.D0
         DO J=1,3
            S(I) = S(I) + R(J)*HTM1(J,I)
         END DO
      END DO
C     
C      S(1)=S(1)-NINT(S(1))
C      S(2)=S(2)-NINT(S(2))
C      S(3)=S(3)-NINT(S(3))
      DX=S(1) 
      DY=S(2) 
      DZ=S(3) 
C     
      RETURN
      END
C====================================================================
      SUBROUTINE ABCTOXYZ(DX,DY,DZ,HT,HTM1)
      IMPLICIT REAL*8 (A-H,O-Z) 
      DIMENSION R(3),S(3),HT(3,3),HTM1(3,3)
C     
      S(1)=DX 
      S(2)=DY 
      S(3)=DZ 
C     ABC->XYZ ; A = X * M
      DO I=1,3
         R(I) = 0.D0
         DO J=1,3
            R(I) = R(I) + S(J)*HT(J,I)
         END DO
      END DO
C     
      DX=R(1) 
      DY=R(2) 
      DZ=R(3) 
C     
      RETURN
      END

