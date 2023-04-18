! Copyright (C) 2022 TurboRVB group
!
! This program is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with this program. If not, see <http://www.gnu.org/licenses/>.

! SLATEC Common Mathematical Library, Version 4.1, July 1993
! a comprehensive software library containing over
! 1400 general purpose mathematical and statistical routines
! written in Fortran 77.

subroutine DSORT(X, Y, N, KFLAG)
    !***BEGIN PROLOGUE  DSORT
    !***DATE WRITTEN   761101   (YYMMDD)
    !***REVISION DATE  820801   (YYMMDD)
    !***CATEGORY NO.  N6A2B1
    !***KEYWORDS  QUICKSORT,SINGLETON QUICKSORT,SORT,SORTING
    !***AUTHOR  JONES, R. E., (SNLA)
    !           WISNIEWSKI, J. A., (SNLA)
    !***PURPOSE  SSORT sorts array X and optionally makes the same
    !            interchanges in array Y.  The array X may be sorted in
    !            increasing order or decreasing order.  A slightly modified
    !            QUICKSORT algorithm is used.
    !***DESCRIPTION
    !
    !     Written by Rondall E. Jones
    !     Modified by John A. Wisniewski to use the Singleton quicksort
    !     algorithm.  Date 18 November 1976.
    !
    !     Abstract
    !         SSORT sorts array X and optionally makes the same
    !         interchanges in array Y.  The array X may be sorted in
    !         increasing order or decreasing order.  A slightly modified
    !         quicksort algorithm is used.
    !
    !     Reference
    !         Singleton, R. C., Algorithm 347, An Efficient Algorithm for
    !         Sorting with Minimal Storage, CACM,12(3),1969,185-7.
    !
    !     Description of Parameters
    !         X - array of values to be sorted   (usually abscissas)
    !         Y - array to be (optionally) carried along
    !         N - number of values in array X to be sorted
    !         KFLAG - control parameter
    !             =2  means sort X in increasing order and carry Y along.
    !             =1  means sort X in increasing order (ignoring Y)
    !             =-1 means sort X in decreasing order (ignoring Y)
    !             =-2 means sort X in decreasing order and carry Y along.
    !***REFERENCES  SINGLETON,R.C., ALGORITHM 347, AN EFFICIENT ALGORITHM
    !                 FOR SORTING WITH MINIMAL STORAGE, CACM,12(3),1969,
    !                 185-7.
    !***ROUTINES CALLED  XERROR
    !***END PROLOGUE  SSORT
    real*8 x, y
    dimension X(N), Y(N), IL(21), IU(21)
    !***FIRST EXECUTABLE STATEMENT  SSORT
    NN = N
    if (NN .ge. 1) GO TO 10
    !      CALL XERROR ( 'SSORT- THE NUMBER OF VALUES TO BE SORTED WAS NOT PO
    !     1SITIVE.',58,1,1)
    write (*, *) 'SSORT- THE NUMBER OF VALUES TO BE SORTED WAS NOT >0'
    return
10  KK = IABS(KFLAG)
    if ((KK .eq. 1) .or. (KK .eq. 2)) GO TO 15
    write (*, *) 'SSORT- THE SORT !ONTROL PARAMETER, K, WAS NOT 2, 1,&
            & -1, OR -2.'
    !      CALL XERROR ( 'SSORT- THE SORT CONTROL PARAMETER, K, WAS NOT 2, 1,
    !     1 -1, OR -2.',62,2,1)
    return
    !
    ! ALTER ARRAY X TO GET DECREASING ORDER IF NEEDED
    !
15  if (KFLAG .ge. 1) GO TO 30
    do 20 I = 1, NN
20      X(I) = -X(I)
30      GO TO(100, 200), KK
        !
        ! SORT X ONLY
        !
100     continue
        M = 1
        I = 1
        J = NN
        R = .375
110     if (I .eq. J) GO TO 155
115     if (R .gt. .5898437) GO TO 120
        R = R + 3.90625e-2
        GO TO 125
120     R = R - .21875
125     K = I
        !                                  SELECT A CENTRAL ELEMENT OF THE
        !                                  ARRAY AND SAVE IT IN LOCATION T
        IJ = I + IFIX(FLOAT(J - I)*R)
        T = X(IJ)
        !                                  IF FIRST ELEMENT OF ARRAY IS GREATER
        !                                  THAN T, INTERCHANGE WITH T
        if (X(I) .le. T) GO TO 130
        X(IJ) = X(I)
        X(I) = T
        T = X(IJ)
130     L = J
        !                                  IF LAST ELEMENT OF ARRAY IS LESS THAN
        !                                  T, INTERCHANGE WITH T
        if (X(J) .ge. T) GO TO 140
        X(IJ) = X(J)
        X(J) = T
        T = X(IJ)
        !                                  IF FIRST ELEMENT OF ARRAY IS GREATER
        !                                  THAN T, INTERCHANGE WITH T
        if (X(I) .le. T) GO TO 140
        X(IJ) = X(I)
        X(I) = T
        T = X(IJ)
        GO TO 140
135     TT = X(L)
        X(L) = X(K)
        X(K) = TT
        !                                  FIND AN ELEMENT IN THE SECOND HALF OF
        !                                  THE ARRAY WHICH IS SMALLER THAN T
140     L = L - 1
        if (X(L) .gt. T) GO TO 140
        !                                  FIND AN ELEMENT IN THE FIRST HALF OF
        !                                  THE ARRAY WHICH IS GREATER THAN T
145     K = K + 1
        if (X(K) .lt. T) GO TO 145
        !                                  INTERCHANGE THESE ELEMENTS
        if (K .le. L) GO TO 135
        !                                  SAVE UPPER AND LOWER SUBSCRIPTS OF
        !                                  THE ARRAY YET TO BE SORTED
        if (L - I .le. J - K) GO TO 150
        IL(M) = I
        IU(M) = L
        I = K
        M = M + 1
        GO TO 160
150     IL(M) = K
        IU(M) = J
        J = L
        M = M + 1
        GO TO 160
        !                                  BEGIN AGAIN ON ANOTHER PORTION OF
        !                                  THE UNSORTED ARRAY
155     M = M - 1
        if (M .eq. 0) GO TO 300
        I = IL(M)
        J = IU(M)
160     if (J - I .ge. 1) GO TO 125
        if (I .eq. 1) GO TO 110
        I = I - 1
165     I = I + 1
        if (I .eq. J) GO TO 155
        T = X(I + 1)
        if (X(I) .le. T) GO TO 165
        K = I
170     X(K + 1) = X(K)
        K = K - 1
        if (T .lt. X(K)) GO TO 170
        X(K + 1) = T
        GO TO 165
        !
        ! SORT X AND CARRY Y ALONG
        !
200     continue
        M = 1
        I = 1
        J = NN
        R = .375
210     if (I .eq. J) GO TO 255
215     if (R .gt. .5898437) GO TO 220
        R = R + 3.90625e-2
        GO TO 225
220     R = R - .21875
225     K = I
        !                                  SELECT A CENTRAL ELEMENT OF THE
        !                                  ARRAY AND SAVE IT IN LOCATION T
        IJ = I + IFIX(FLOAT(J - I)*R)
        T = X(IJ)
        TY = Y(IJ)
        !                                  IF FIRST ELEMENT OF ARRAY IS GREATER
        !                                  THAN T, INTERCHANGE WITH T
        if (X(I) .le. T) GO TO 230
        X(IJ) = X(I)
        X(I) = T
        T = X(IJ)
        Y(IJ) = Y(I)
        Y(I) = TY
        TY = Y(IJ)
230     L = J
        !                                  IF LAST ELEMENT OF ARRAY IS LESS THAN
        !                                  T, INTERCHANGE WITH T
        if (X(J) .ge. T) GO TO 240
        X(IJ) = X(J)
        X(J) = T
        T = X(IJ)
        Y(IJ) = Y(J)
        Y(J) = TY
        TY = Y(IJ)
        !                                  IF FIRST ELEMENT OF ARRAY IS GREATER
        !                                  THAN T, INTERCHANGE WITH T
        if (X(I) .le. T) GO TO 240
        X(IJ) = X(I)
        X(I) = T
        T = X(IJ)
        Y(IJ) = Y(I)
        Y(I) = TY
        TY = Y(IJ)
        GO TO 240
235     TT = X(L)
        X(L) = X(K)
        X(K) = TT
        TTY = Y(L)
        Y(L) = Y(K)
        Y(K) = TTY
        !                                  FIND AN ELEMENT IN THE SECOND HALF OF
        !                                  THE ARRAY WHICH IS SMALLER THAN T
240     L = L - 1
        if (X(L) .gt. T) GO TO 240
        !                                  FIND AN ELEMENT IN THE FIRST HALF OF
        !                                  THE ARRAY WHICH IS GREATER THAN T
245     K = K + 1
        if (X(K) .lt. T) GO TO 245
        !                                  INTERCHANGE THESE ELEMENTS
        if (K .le. L) GO TO 235
        !                                  SAVE UPPER AND LOWER SUBSCRIPTS OF
        !                                  THE ARRAY YET TO BE SORTED
        if (L - I .le. J - K) GO TO 250
        IL(M) = I
        IU(M) = L
        I = K
        M = M + 1
        GO TO 260
250     IL(M) = K
        IU(M) = J
        J = L
        M = M + 1
        GO TO 260
        !                                  BEGIN AGAIN ON ANOTHER PORTION OF
        !                                  THE UNSORTED ARRAY
255     M = M - 1
        if (M .eq. 0) GO TO 300
        I = IL(M)
        J = IU(M)
260     if (J - I .ge. 1) GO TO 225
        if (I .eq. 1) GO TO 210
        I = I - 1
265     I = I + 1
        if (I .eq. J) GO TO 255
        T = X(I + 1)
        TY = Y(I + 1)
        if (X(I) .le. T) GO TO 265
        K = I
270     X(K + 1) = X(K)
        Y(K + 1) = Y(K)
        K = K - 1
        if (T .lt. X(K)) GO TO 270
        X(K + 1) = T
        Y(K + 1) = TY
        GO TO 265
        !
        ! CLEAN UP
        !
300     if (KFLAG .ge. 1) return
        do 310 I = 1, NN
310         X(I) = -X(I)
            return
        end
