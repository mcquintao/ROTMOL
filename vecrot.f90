PROGRAM VECROT
    IMPLICIT NONE
    INTEGER i, j
    REAL XYZ, VEC
    DIMENSION XYZ(3,3), VEC(3)

    DO 1 i = 1, 3
        VEC(i) = 100
        DO 2 j = 1, 3
            XYZ(i,j) = i+j
        2 CONTINUE
    1 CONTINUE
    
    XYZ(1,:) = VEC(:)
    
    DO 3 i = 1,3
        PRINT *, (XYZ(i,j), j=1,3)
    3 CONTINUE

END PROGRAM
