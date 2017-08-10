PROGRAM ROTMOL

IMPLICIT NONE
INTEGER NATOMS
REAL PHI, THETA, GAMA
REAL, ALLOCATABLE :: XYZ(:,:)
CHARACTER*3, ALLOCATABLE :: ATOM(:)

OPEN(UNIT=1, FILE="input")
OPEN(UNIT=2, FILE="output")

READ(1,*) NATOMS, PHI, THETA, GAMA
ALLOCATE(XYZ(1:NATOMS,1:3),ATOM(1:NATOMS))

CALL READMOL(NATOMS,ATOM,XYZ)
CALL WRITEMOL(NATOMS,ATOM,XYZ)

CLOSE(1)
CLOSE(2)

CONTAINS

SUBROUTINE READMOL(NATOMS,ATOM,XYZ)
IMPLICIT NONE
INTEGER i, j, NATOMS
REAL XYZ
CHARACTER*3, ATOM
DIMENSION XYZ(:,:), ATOM(:)

DO 10 i = 1, NATOMS
    READ(1,*) ATOM(i), (XYZ(i,j), j=1,3)
10 CONTINUE

END SUBROUTINE

SUBROUTINE WRITEMOL(NATOMS,ATOM,XYZ)
IMPLICIT NONE
INTEGER i, j, NATOMS
REAL XYZ
CHARACTER*3, ATOM
DIMENSION XYZ(:,:), ATOM(:)

WRITE(2,*) NATOMS
WRITE(2,*) "INPUT"
DO 20 i = 1, NATOMS
    WRITE(2,*) ATOM(i), (XYZ(i,j), j=1,3)
20 CONTINUE

END SUBROUTINE

END PROGRAM ROTMOL
