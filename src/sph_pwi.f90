!! SPH_pwi: Routines related to point particle-wall interactions.
!! sph_pwi-0.9.0
!!
!! Copyright (C) 2015-2021  Chris D. Dritselis
!!
!! Laboratory of Fluid Mechanics & Turbomachinery
!! Mechanical Engineering Department
!! University of Thessaly
!! Pedion Areos, 38334 Volos, GR
!!
!! This program is free software; you can redistribute it and/or
!! modify it under the terms of the GNU General Public License
!! as published by the Free Software Foundation; either version 2
!! of the License, or (at your option) any later version.
!! This program is distributed in the hope that it will be useful,
!! but WITHOUT ANY WARRANTY; without even the implied warranty of
!! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!! GNU General Public License for more details.
!!
!! You should have received a copy of the GNU General Public License
!! along with this program; if not, write to the Free Software
!! Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!                                                          !!
!!                 MODULE SPHBasics_module                  !!
!!                                                          !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
MODULE SPHBasics_module
USE Channelflow_mod
  IMPLICIT NONE

  INTEGER, PARAMETER :: iWP_sph   = SELECTED_INT_KIND( 12 )            ! Double integer
  INTEGER, PARAMETER :: lWP_sph   = SELECTED_INT_KIND( 12 )            ! Double logical
  INTEGER, PARAMETER :: WP_sph    = SELECTED_REAL_KIND( 15, 307 )      ! Double real
!  INTEGER, PARAMETER :: WP_sph    = SELECTED_REAL_KIND( 25, 1000 )     ! Quad real
  INTEGER, PARAMETER :: WPinc_sph = SELECTED_REAL_KIND( 25, 1000 )     ! Quad real

  REAL ( KIND = WP_sph ), PARAMETER :: SMALL_NUMBER  = 0.000000000000001_WP_sph !0.00000001_WP_sph
  REAL ( KIND = WP_sph ), PARAMETER :: SMALL_NUMBER1 = 0.000000000000001_WP_sph

END MODULE SPHBasics_module



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!                                                          !!
!!                     MODULE SPH_pwi_module                !!
!!                                                          !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
MODULE SPH_pwi_module
  USE SPHBasics_module
  IMPLICIT NONE

  !!   basitype   !!

  TYPE, PUBLIC ::  basicType
     REAL ( KIND = WP_sph ) :: x, y, z
  END TYPE

  !!   point    !!

  TYPE, EXTENDS ( basicType ), PUBLIC ::  point

  CONTAINS
     PROCEDURE, PUBLIC :: assignP1toP0
     PROCEDURE, PUBLIC :: assignAtoP
     PROCEDURE, PUBLIC :: distancePointToPoint
     PROCEDURE, PUBLIC :: testEqualP0P1
     PROCEDURE, PUBLIC :: testNOTEqualP0P1
     PROCEDURE, PUBLIC :: translatePointByVector
     PROCEDURE, PUBLIC :: pointsDifference
     GENERIC, PUBLIC   :: assignment(=)      => assignP1toP0, assignAtoP
     GENERIC, PUBLIC   :: operator(.D.)      => distancePointToPoint
     GENERIC, PUBLIC   :: operator(.EQUAL.)  => testEqualP0P1
     GENERIC, PUBLIC   :: operator(.notEQUAL.) => testNOTEqualP0P1
     GENERIC, PUBLIC   :: operator(+)        => translatePointByVector
     GENERIC, PUBLIC   :: operator(-)        => pointsDifference

  END TYPE point

  !!   vector   !!

  TYPE, EXTENDS ( basicType ), PUBLIC ::  vector

  CONTAINS
     PROCEDURE, PUBLIC :: assignV1toV0
     PROCEDURE, PUBLIC :: assignAtoV
     PROCEDURE, PUBLIC :: testEqualV0V1
     PROCEDURE, PUBLIC :: add
     PROCEDURE, PUBLIC :: minus
     PROCEDURE, PUBLIC :: scaleVector
     PROCEDURE, PUBLIC, PASS(v1) :: scaleVectorByFactor
     PROCEDURE, PUBLIC :: cross
     PROCEDURE, PUBLIC :: dotProduct
     PROCEDURE, PUBLIC :: perpProduct
     PROCEDURE, PUBLIC :: norm
     PROCEDURE, PUBLIC :: distance
     GENERIC, PUBLIC   :: assignment(=)      => assignV1toV0, assignAtoV
     GENERIC, PUBLIC   :: operator(.EQUAL.)  => testEqualV0V1
     GENERIC, PUBLIC   :: operator(+)        => add
     GENERIC, PUBLIC   :: operator(-)        => minus
     GENERIC, PUBLIC   :: operator(*)        => scaleVector, scaleVectorByFactor, cross
     GENERIC, PUBLIC   :: operator(.X.)      => cross
     GENERIC, PUBLIC   :: operator(.DOT.)    => dotProduct
     GENERIC, PUBLIC   :: operator(.PERP.)   => perpProduct
     GENERIC, PUBLIC   :: operator(.D.)      => distance

  END TYPE vector

  !!   SPHLine   !!

  TYPE, PUBLIC :: SPHLine
     TYPE ( point ) :: p0, p1
  END TYPE SPHLine

  !!   SPHSegment   !!

   TYPE, EXTENDS( SPHLine ), PUBLIC :: SPHSegment
   END TYPE SPHSegment

  !!   SPHPlane   !!

  TYPE, PUBLIC :: SPHPlane
     TYPE ( point )  :: plp
     TYPE ( vector ) :: n
  END TYPE SPHPlane

  !!   SPHTrack   !!

  TYPE, PUBLIC :: SPHTrack
     TYPE ( point )  :: p0
     TYPE ( vector ) :: v
  END TYPE SPHTrack

  !!   SPHPolygon   !!

  TYPE, PUBLIC :: SPHPolygon
     INTEGER ( KIND = WP_sph )                             :: nVertex
     TYPE    ( point )        , ALLOCATABLE, DIMENSION (:) :: V
  END TYPE SPHPolygon



CONTAINS


!!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!!
!! < POINT type

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!                                                          !!
  !!                 (POINT)-assignP1toP0                     !!
  !!                                                          !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !                                                            !
  !     assignP1toP0: assign point P0 to P1                    !
  !     Input : points P0 and P1                               !
  !     Return: set P0 equal to P1                             !
  !                                                            !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE assignP1toP0( p0, p1 )
    CLASS ( point ) , INTENT ( OUT )  :: p0
    CLASS ( point ) , INTENT ( IN )   :: p1
    p0 % x = p1 % x
    p0 % y = p1 % y
    p0 % z = p1 % z
  END SUBROUTINE assignP1toP0



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!                                                          !!
  !!                   (POINT)-assignAtoP                     !!
  !!                                                          !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !                                                            !
  !     assignAtoP: assign array a to P0                       !
  !     Input : point P and array a                            !
  !     Return: set P equal to A                               !
  !                                                            !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE assignAtoP( p, a )
  CLASS ( point )        , INTENT ( OUT )                :: p
  REAL  ( KIND = WP_sph ), INTENT ( IN ) , DIMENSION (:) :: a
  p % x = a(1)
  p % y = a(2)
  p % z = a(3)
END SUBROUTINE assignAtoP



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!                                                          !!
  !!               (POINT)-distancePointToPoint               !!
  !!                                                          !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !                                                            !
  !     distancePointToPoint: distance between point P0 and P0 !
  !     Input : points P0 and P1                               !
  !     Return: d (distance)                                   !
  !                                                            !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  REAL ( KIND = WP_sph ) FUNCTION distancePointToPoint( p0, p1 ) RESULT ( d )
    CLASS ( point ), INTENT ( IN ) :: p0, p1
    d = SQRT( (p1 % x - p0 % x) * (p1 % x - p0 % x)     &
        +     (p1 % y - p0 % y) * (p1 % y - p0 % y)     &
        +     (p1 % z - p0 % z) * (p1 % z - p0 % z)  )
  END FUNCTION distancePointToPoint



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!                                                          !!
  !!                   (POINT)-testEqualP0P1                  !!
  !!                                                          !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !                                                            !
  !     testEqualP0P1: test equality for points P0 to P1       !
  !     Input : points P0 and P1                               !
  !     Return: TRUE when P0==P11                              !
  !             FALSE when P0/=P11                             !
  !                                                            !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  LOGICAL ( KIND = lWP_sph ) FUNCTION testEqualP0P1( p0, p1 ) RESULT ( test )
    CLASS ( point ), INTENT ( IN )  :: p0, p1
    IF (  (ABS(p0 % x - p1 % x) <= SMALL_NUMBER) .AND.              &
          (ABS(p0 % y - p1 % y) <= SMALL_NUMBER) .AND.              &
          (ABS(p0 % z - p1 % z) <= SMALL_NUMBER) ) test = .TRUE.
    IF (  (ABS(p0 % x - p1 % x) > SMALL_NUMBER) .OR.                &
          (ABS(p0 % y - p1 % y) > SMALL_NUMBER) .OR.                &
          (ABS(p0 % z - p1 % z) > SMALL_NUMBER) ) test = .FALSE.
  END FUNCTION testEqualP0P1


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!                                                          !!
!!                   (POINT)-testNOTEqualP0P1               !!
!!                                                          !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                            !
!     testEqualP0P1: test NOT equality for points P0 to P1   !
!     Input : points P0 and P1                               !
!     Return: TRUE when P0/=P11                              !
!             FALSE when P0==P11                             !
!                                                            !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
LOGICAL ( KIND = lWP_sph ) FUNCTION testNOTEqualP0P1( p0, p1 ) RESULT ( test )
CLASS ( point ), INTENT ( IN )  :: p0, p1
IF (  (ABS(p0 % x - p1 % x) <= SMALL_NUMBER) .AND.              &
(ABS(p0 % y - p1 % y) <= SMALL_NUMBER) .AND.              &
(ABS(p0 % z - p1 % z) <= SMALL_NUMBER) ) test = .FALSE.
IF (  (ABS(p0 % x - p1 % x) > SMALL_NUMBER) .OR.                &
(ABS(p0 % y - p1 % y) > SMALL_NUMBER) .OR.                &
(ABS(p0 % z - p1 % z) > SMALL_NUMBER) ) test = .TRUE.
END FUNCTION testNOTEqualP0P1

!-----------------------------------------------------------------------------------!
#if experimental_
TYPE ( point ) FUNCTION translatePointByVectorBasic( p, v ) RESULT ( pp )
CLASS ( point )    , INTENT ( IN )  :: p
CLASS ( basicType ), INTENT ( IN )  :: v
pp % x = p % x + v % x
pp % y = p % y + v % y
pp % z = p % z + v % z
END FUNCTION translatePointByVectorBasic
#endif
!-----------------------------------------------------------------------------------!


!! end of POINT type >
!!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!!



!!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!!
!! < VECTOR type


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!                                                          !!
  !!                   (VECTOR)-assignV1toV0                  !!
  !!                                                          !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !                                                            !
  !     assignV1toV0: assign vector V1 to vector V2            !
  !     Input : vectors V1 and V2                              !
  !     Return: set V1 equal to V2                             !
  !                                                            !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE assignV1toV0( v1, v2 )
    CLASS ( vector ), INTENT ( OUT ) :: v1
    CLASS ( vector ), INTENT ( IN )  :: v2
    v1 % x = v2 % x
    v1 % y = v2 % y
    v1 % z = v2 % z
  END SUBROUTINE assignV1toV0



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!                                                          !!
  !!                   (VECTOR)-assignAtoV                    !!
  !!                                                          !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !                                                            !
  !     assignAtoV: assign array A to vector V                 !
  !     Input : vector V and array a                           !
  !     Return: set V equal to A                               !
  !                                                            !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE assignAtoV( v, a )
    CLASS ( vector )       , INTENT ( OUT )               :: v
    REAL  ( KIND = WP_sph ), INTENT ( IN ), DIMENSION (:) :: a
    v % x = a( 1 )
    v % y = a( 2 )
    v % z = a( 3 )
  END SUBROUTINE assignAtoV



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!                                                          !!
  !!                   (VECTOR)-testEqualV0V1                 !!
  !!                                                          !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !                                                            !
  !     testEqualV0V1: test equality for vectors V0 to V1      !
  !     Input : vectors V0 and V1                              !
  !     Return: TRUE when V0==V11                              !
  !             FALSE when V0/=V11                             !
  !                                                            !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  LOGICAL ( KIND = lWP_sph ) FUNCTION testEqualV0V1( v0, v1 ) RESULT ( test )
    CLASS ( vector ), INTENT ( IN )  :: v0, v1
    IF (  (ABS(v0 % x - v1 % x) <= SMALL_NUMBER) .AND.              &
          (ABS(v0 % y - v1 % y) <= SMALL_NUMBER) .AND.              &
          (ABS(v0 % z - v1 % z) <= SMALL_NUMBER) ) test = .TRUE.
    IF (  (ABS(v0 % x - v1 % x) > SMALL_NUMBER) .OR.                &
          (ABS(v0 % y - v1 % y) > SMALL_NUMBER) .OR.                &
          (ABS(v0 % z - v1 % z) > SMALL_NUMBER) ) test = .FALSE.
    END FUNCTION testEqualV0V1



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!                                                          !!
  !!                      (VECTOR)-add                        !!
  !!                                                          !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !                                                            !
  !     add   : add vectors V1 and V2                          !
  !     Input : vectors V1 and V2                              !
  !     Return: addition of  V1 and V2                         !
  !                                                            !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  TYPE ( vector ) FUNCTION add( v1, v2 ) RESULT ( v3 )
    CLASS ( vector ), INTENT ( IN )  :: v1, v2
    v3 % x = v1 % x + v2 % x
    v3 % y = v1 % y + v2 % y
    v3 % z = v1 % z + v2 % z
  END FUNCTION add



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!                                                          !!
  !!                     (VECTOR)-minus                       !!
  !!                                                          !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !                                                            !
  !     minus : substract vectors V1 and V2                    !
  !     Input : vectors V1 and V2                              !
  !     Return: substraction of V1 and V2                      !
  !                                                            !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  TYPE ( vector ) FUNCTION minus( v1, v2 ) RESULT ( v3 )
   CLASS ( vector ), INTENT ( IN ) :: v1, v2
   v3 % x = v1 % x - v2 % x
   v3 % y = v1 % y - v2 % y
   v3 % z = v1 % z - v2 % z
  END FUNCTION minus



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!                                                          !!
  !!                  (VECTOR)-scaleVector                    !!
  !!                                                          !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !                                                            !
  !     scaleVector: scale vector V1 by scalar s               !
  !     Input : vector V1 and scalar S                         !
  !     Return: V1 * S                                         !
  !                                                            !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  TYPE ( vector ) FUNCTION scaleVector( v1, s ) RESULT ( v2 )
    CLASS ( vector )      , INTENT ( IN ) :: v1
    REAL ( KIND = WP_sph ), INTENT ( IN ) :: s
    v2 % x = s * v1 % x
    v2 % y = s * v1 % y
    v2 % z = s * v1 % z
  END FUNCTION scaleVector



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!                                                          !!
  !!                (VECTOR)-scaleVectorByFactor              !!
  !!                                                          !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !                                                            !
  !     scaleVectorByFactor: scale vector V1 by scalar s       !
  !     Input : vector V1 and scalar S                         !
  !     Return: S * V1                                         !
  !                                                            !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  TYPE ( vector ) FUNCTION scaleVectorByFactor( s, v1 ) RESULT ( v2 )
    CLASS ( vector )      , INTENT ( IN ) :: v1
    REAL ( KIND = WP_sph ), INTENT ( IN ) :: s
    v2 % x = s * v1 % x
    v2 % y = s * v1 % y
    v2 % z = s * v1 % z
  END FUNCTION scaleVectorByFactor



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!                                                          !!
  !!                      (VECTOR)-cross                      !!
  !!                                                          !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !                                                            !
  !     cross: 3D cross product between vectors V1 and V2      !
  !     Input : vectors V1 and V2                              !
  !     Return: V1 x V2                                        !
  !                                                            !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  TYPE ( vector ) FUNCTION cross( v1, v2 ) RESULT ( v3 )
    CLASS ( vector ), INTENT ( IN ) :: v1, v2
    v3 % x = v1 % y * v2 % z - v1 % z * v2 % y
    v3 % y = v1 % z * v2 % x - v1 % x * v2 % z
    v3 % z = v1 % x * v2 % y - v1 % y * v2 % x
  END FUNCTION cross



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!                                                          !!
  !!                      (VECTOR)-dotProduct                 !!
  !!                                                          !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !                                                            !
  !     dotProduct: 3D dot product between vectors V1 and V2   !
  !     Input : vectors V1 and V2                              !
  !     Return: V1 . V2                                        !
  !                                                            !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  REAL (KIND = WP_sph ) FUNCTION dotProduct( v1, v2 ) RESULT ( v3 )
    CLASS ( vector ), INTENT ( IN ) :: v1, v2
    v3 = v1 % x * v2 % x + v1 % y * v2 % y + v1 % z * v2 % z
  END FUNCTION dotProduct



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!                                                          !!
  !!                      (VECTOR)-perpProduct                !!
  !!                                                          !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !                                                            !
  !     perpProduct: 2D perpedicular product between           !
  !                  vectors V1 and V2.                        !
  !     Input : vectors V1 and V2                              !
  !     Return: V1 _|_ V2                                      !
  !                                                            !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  REAL (KIND = WP_sph ) FUNCTION perpProduct( v1, v2 ) RESULT ( v3 )
    CLASS ( vector ), INTENT ( IN ) :: v1, v2
    v3 = v1 % x * v2 % y - v1 % y * v2 % x
  END FUNCTION perpProduct



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!                                                          !!
  !!                        (VECTOR)-norm                     !!
  !!                                                          !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !                                                            !
  !     norm: norm of vector V (length of vector V)            !
  !     Input : vector V                                       !
  !     Return: norm of V                                      !
  !                                                            !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  REAL (KIND = WP_sph ) FUNCTION norm( v ) RESULT ( normV )
    CLASS ( vector ), INTENT ( IN ) :: v
    normV = SQRT( v % x * v % x + v % y * v % y + v % z * v % z )
  END FUNCTION norm



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!                                                          !!
  !!                      (VECTOR)-distance                   !!
  !!                                                          !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !                                                            !
  !     perpProduct: distance between vectors V1 and V2 =      !
  !                  norm(V1-V2) (length of vector V1-V2)      !
  !     Input : vectors V1 and V2                              !
  !     Return: distance d = norm(V1-V2)                       !
  !                                                            !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  REAL (KIND = WP_sph ) FUNCTION distance( v1, v2 ) RESULT ( d )
    CLASS ( vector ), INTENT ( IN ) :: v1, v2
    d = SQRT( (v1 % x - v2 % x) * (v1 % x - v2 % x)      &
        +     (v1 % y - v2 % y) * (v1 % y - v2 % y)      &
        +     (v1 % z - v2 % z) * (v1 % z - v2 % z) )
  END FUNCTION distance


!-----------------------------------------------------------------------------------!
#if experimental_
TYPE ( vector ) FUNCTION pointsDifferenceBasic( p1, p0 ) RESULT ( v )
CLASS ( basicType ), INTENT ( IN )  :: p1, p0
v % x = p1 % x - p0 % x
v % y = p1 % y - p0 % y
v % z = p1 % z - p0 % z
END FUNCTION pointsDifferenceBasic
#endif
!-----------------------------------------------------------------------------------!



!! end of VECTOR type >
!!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!!


!!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!!
!! < POINT/VECTOR mixed type


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!                                                          !!
  !!              (POINT/VECTOR)-translatePointByVector       !!
  !!                                                          !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !                                                            !
  !     translatePointByVector: translate point P by a vector  !
  !                             V to obtain a new point.       !
  !     Input : point P and vector V                           !
  !     Return: a new point                                    !
  !                                                            !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  TYPE ( point ) FUNCTION translatePointByVector( p, v ) RESULT ( pp )
    CLASS ( point ) , INTENT ( IN )  :: p
    CLASS ( vector ), INTENT ( IN )  :: v
    pp % x = p % x + v % x
    pp % y = p % y + v % y
    pp % z = p % z + v % z
  END FUNCTION translatePointByVector



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!                                                          !!
  !!                 (POINT/VECTOR)-pointsDifference          !!
  !!                                                          !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !                                                            !
  !     pointsDifference: the difference between two points    !
  !                       P0 and P1 may be considered vector   !
  !     Input : points P0 and P1                               !
  !     Return: a new vector                                   !
  !                                                            !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  TYPE ( vector ) FUNCTION pointsDifference( p1, p0 ) RESULT ( v )
    CLASS ( point ), INTENT ( IN )  :: p1, p0
    v % x = p1 % x - p0 % x
    v % y = p1 % y - p0 % y
    v % z = p1 % z - p0 % z
  END FUNCTION pointsDifference


!! end of POINT/VECTOR mixed type >
!!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!!


!!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!!
!! < LINE/PLANE/TRACK type

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!                                                          !!
  !!                    distance_PointToLine                  !!
  !!                                                          !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !                                                            !
  !     distance_PointToLine: get the distance of a point      !
  !                           to a line.                       !
  !     Input : a Point P and a Line L (in any dimension)      !
  !     Return: the shortest distance from P to L              !
  !                                                            !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  REAL ( KIND = WP_sph ) FUNCTION distance_PointToLine( P, L ) RESULT ( dist )
    CLASS ( point )        , INTENT ( IN ) :: P
    CLASS ( SPHLine )      , INTENT ( IN ) :: L
    TYPE  ( vector )                       :: v_, w_
    TYPE  ( point )                        :: p_
    REAL  ( KIND = WP_sph )                :: c1_, c2_, b_


    v_   = L % P1 - L % P0
    w_   = P - L % P0
    c1_  = dotProduct( w_, v_ )
    c2_  = dotProduct( v_, v_ )
    b_   = c1_ / c2_
    p_   = L % P0 + v_ * b_
    dist = distancePointToPoint( P, p_ )


  END FUNCTION distance_PointToLine



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!                                                          !!
  !!                   distance_PointToSegment                !!
  !!                                                          !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !                                                            !
  !     distance_PointToSegment: get the distance of a point   !
  !                           to a segment.                    !
  !     Input : a Point P and a Segment S (in any dimension)   !
  !     Return: the shortest distance from P to S              !
  !                                                            !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  REAL ( KIND = WP_sph ) FUNCTION distance_PointToSegment( P, S ) RESULT ( dist )
    CLASS ( point )        , INTENT ( IN ) :: P
    CLASS ( SPHSegment )   , INTENT ( IN ) :: S
    TYPE  ( vector )                       :: v_, w_
    TYPE  ( point )                        :: p_
    REAL  ( KIND = WP_sph )                :: c1_, c2_, b_


    v_   = S % P1 - S % P0
    w_   = P - S % P0
    c1_  = dotProduct( w_, v_ )
    IF ( c1_ <= 0.0_WP_sph ) THEN
       dist = distancePointToPoint( P, S % P0 )
       RETURN
    END IF
    c2_  = dotProduct( v_, v_ )
    IF ( c2_ <= c1_ ) THEN
       dist = distancePointToPoint( P, S % P1 )
       RETURN
    END IF

    b_   = c1_ / c2_
    p_   = S % P0 + v_ * b_
    dist = distancePointToPoint( P, p_ )


  END FUNCTION distance_PointToSegment



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!                                                          !!
  !!                   distance_PointToPlane                  !!
  !!                                                          !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !                                                            !
  !     distance_PointToPlane: get the distance of a point     !
  !                           to a plane  .                    !
  !     Input : a 3D Point P and a Plane PL                    !
  !     Output: PB base point on PL of perpedicular from P     !
  !     Return: the distance from P to the plane PL            !
  !                                                            !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE distance_PointToPlane( P, PL, PB, dist )
    CLASS ( point )        , INTENT ( IN ) :: P
    CLASS ( SPHPlane )     , INTENT ( IN ) :: PL
    CLASS ( point )        , INTENT ( OUT ):: PB
    REAL  ( KIND = WP_sph ), INTENT ( OUT ):: dist
    REAL  ( KIND = WP_sph )                :: sn_, sd_, sb_


    sn_ = - dotProduct( PL % n, ( P - PL % PLP ) )
    sd_ =   dotProduct( PL % n, PL % n )
    sb_ = sn_ / sd_
    PB = P + sb_ * PL % n
    dist = distancePointToPoint( P, PB )


  END SUBROUTINE distance_PointToPlane



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!                                                          !!
  !!                   distance3D_LineToLine                  !!
  !!                                                          !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !                                                            !
  !     distance3D_LineToLine: get the 3D minimum              !
  !                            between two lines.              !
  !     Input : a 3D Lines L1 and L2                           !
  !     Return: the shortest distance between L1 and L2        !
  !                                                            !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  REAL ( KIND = WP_sph ) FUNCTION distance3D_LineToLine( L1, L2 ) RESULT ( dist )
    CLASS ( SPHLine )      , INTENT ( IN ) :: L1, L2
    TYPE  ( vector )                       :: u_, v_, w_, dp_
    REAL  ( KIND = WP_sph )                :: a_, b_, c_, d_, e_, DD_, sc_, tc_


    u_ = L1 % P1 - L1 % P0
    v_ = L2 % P1 - L2 % P0
    w_ = L1 % P0 - L2 % P0

    a_ = dotProduct( u_, u_ )   !   always >=0
    b_ = dotProduct( u_, v_ )
    c_ = dotProduct( v_, v_ )   !   always >=0
    d_ = dotProduct( u_, w_ )
    e_ = dotProduct( v_, w_ )
    DD_ = a_ * c_ - b_ * b_     !   always >=0

    !!   compute the line parameters of the two closest points   !!

    IF ( DD_ < SMALL_NUMBER ) THEN       !   the lines are almost parallel

       sc_ = 0.0_WP_sph
       IF ( b_ > c_  ) tc_ = d_ / b_     !   use the largest demominator
       IF ( b_ <= c_ ) tc_ = e_ / c_

    ELSE

       sc_ = ( b_ * e_ - c_ * d_ ) / DD_
       tc_ = ( a_ * e_ - b_ * d_ ) / DD_

    END IF

    !!   get the difference of the two closest points   !!

    dp_ = w_ + ( sc_ * u_ ) - ( tc_ * v_ )  ! L1(sc) - L2(sc)
    dist = norm( dp_ )                      ! return the closest distance


  END FUNCTION distance3D_LineToLine



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!                                                          !!
  !!                distance3D_SegmentToSegment               !!
  !!                                                          !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !                                                            !
  !     distance3D_SegmentToSegment: get the 3D minimum        !
  !                                  between two segments.     !
  !     Input : a 3D line Segments S1 and S2                   !
  !     Return: the shortest distance between S1 and S2        !
  !                                                            !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  REAL ( KIND = WP_sph ) FUNCTION distance3D_SegmentToSegment( S1, S2 ) RESULT ( dist )
    CLASS ( SPHSegment )   , INTENT ( IN ) :: S1, S2
    TYPE ( vector )         :: u_, v_, w_, dp_
    REAL  ( KIND = WP_sph )                :: a_, b_, c_, d_, e_, DD_, sc_, tc_
    REAL  ( KIND = WP_sph )                :: sN_, sD_, tN_, tD_


    u_ = S1 % P1 - S1 % P0
    v_ = S2 % P1 - S2 % P0
    w_ = S1 % P0 - S2 % P0

    a_ = dotProduct( u_, u_ )   !   always >=0
    b_ = dotProduct( u_, v_ )
    c_ = dotProduct( v_, v_ )   !   always >=0
    d_ = dotProduct( u_, w_ )
    e_ = dotProduct( v_, w_ )
    DD_ = a_ * c_ - b_ * b_     !   always >=0

    sD_ = DD_   !   default sD = D >= 0
    tD_ = DD_   !   default tD = D >= 0

    !!   compute the line parameters of the two closest points   !!

    IF ( DD_ < SMALL_NUMBER1 ) THEN    !   the lines are almost parallel

       sN_ = 0.0_WP_sph        !   force using point P0 on segment S1
       sD_ = 1.0_WP_sph        !   to prevent possible division by 0.0 later
       tN_ = e_
       tD_ = c_

    ELSE      ! get the closest points on the infinite lines

       sN_ =  b_ * e_ - c_ * d_
       tN_ =  a_ * e_ - b_ * d_
       IF ( sN_ < 0.0_WP_sph ) THEN    ! sc < 0 => the s=0 edge is visible
          sN_ = 0.0_WP_sph
          tN_ = e_
          tD_ = c_
       ELSE IF ( sN_ > sD_ ) THEN     ! sc > 1  => the s=1 edge is visible
          sN_ = sD_
          tN_ = e_ + b_
          tD_ = c_
       END IF

    END IF


    IF ( tN_ < 0.0_WP_sph ) THEN    ! tc < 0 => the t=0 edge is visible

       tN_ = 0.0_WP_sph

      !!   recompute sc for this edge   !!

      IF ( -d_ < 0.0_WP_sph ) THEN
         sN_ = 0.0_WP_sph
      ELSE IF ( -d_ > a_ ) THEN
         sN_ = sD_
      ELSE
         sN_ = -d_
         sD_ = a_
      END IF

    ELSE IF ( tN_ > tD_ ) THEN   ! tc > 1  => the t=1 edge is visible

       tN_ = tD_

       !!   recompute sc for this edge   !!

       IF ( ( -d_ + b_ ) < 0.0_WP_sph ) THEN
          sN_ = 0.0_WP_sph
       ELSE IF ( ( -d_ + b_ ) > a_ ) THEN
          sN_ = sD_
       ELSE
          sN_ = -d_ + b_
          sD_ = a_
       END IF

    END IF


    !!   finally do the division to get sc and tc   !!

    IF ( ABS( sN_ ) < SMALL_NUMBER1 ) THEN
       sc_ = 0.0_WP_sph
    ELSE
       sc_ = sN_ / sD_
    END IF

    IF ( ABS( tN_ ) < SMALL_NUMBER1 ) THEN
       tc_ = 0.0_WP_sph
    ELSE
       tc_ = tN_ / tD_
    END IF

    !!   get the difference of the two closest points   !!

    dp_ = w_ + ( sc_ * u_ ) - ( tc_ * v_ )  ! S1(sc) - S2(sc)
    dist = norm( dp_ )                      ! return the closest distance


  END FUNCTION distance3D_SegmentToSegment



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!                                                          !!
  !!                         cpa_time                         !!
  !!                                                          !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !                                                            !
  !     cpa_time: compute the time of CPA                      !
  !               for two Tracks.                              !
  !     Input : two Tracks TR1 and TR2                         !
  !     Return: the time at which the two tracks are closest   !
  !                                                            !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  REAL ( KIND = WP_sph ) FUNCTION cpa_time( TR1, TR2 ) RESULT ( cpaTime )
    CLASS ( SPHTrack )   , INTENT ( IN ) :: TR1, TR2
    TYPE ( vector )         :: dv_, w0_
    REAL  ( KIND = WP_sph ) :: dv2_

    dv_ = TR1 % v - TR2 % v
    dv2_ = dotProduct( dv_, dv_ )

    IF ( dv2_ < SMALL_NUMBER ) THEN    ! the tracks are almost parallel
       cpaTime = 0.0_WP_sph                  ! any time is ok. use here time 0.
       RETURN
    END IF

    w0_ = TR1 % P0 - TR2 % P0
    cpaTime = -dotProduct( w0_, dv_ ) / dv2_   ! time of CPA

  END FUNCTION cpa_time



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!                                                          !!
  !!                distance3D_SegmentToSegment               !!
  !!                                                          !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !                                                            !
  !     cpa_distance: compute the distance at CPA              !
  !                   for two Tracks.                          !
  !     Input : two Tracks TR1 and TR2                         !
  !     Return: the distance for which the two tracks          !
  !             are closest                                    !
  !                                                            !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  REAL ( KIND = WP_sph ) FUNCTION cpa_distance( TR1, TR2 ) RESULT ( cpaDistance )
    CLASS ( SPHTrack )   , INTENT ( IN ) :: TR1, TR2
    REAL  ( KIND = WP_sph ) :: cpaTime_
    TYPE  ( point )         :: p1_, p2_

    cpaTime_ = cpa_time( TR1, TR2 )
    p1_ = TR1 % P0 + cpaTime_ * TR1 % v
    p2_ = TR2 % P0 + cpaTime_ * TR2 % v
    cpaDistance = distancePointToPoint( p1_, p2_ )


  END FUNCTION cpa_distance



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!                                                          !!
  !!                          inSegment                       !!
  !!                                                          !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !                                                            !
  !     inSegment: determine if a point is inside a segmanet   !
  !     Input : a point P and a collinear segment S            !
  !     Return: 0 = if P is inside S                           !
  !             1 = if P is not inside S                       !
  !                                                            !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  INTEGER ( KIND = iWP_sph ) FUNCTION inSegment( P, S ) RESULT ( status )
    TYPE ( POINT )     , INTENT ( IN ) :: P
    TYPE ( SPHSegment ), INTENT ( IN ) :: S


    IF ( S % P0 % x /= S % P1 % x ) THEN  ! S is not vertical

       IF ( S % P0 % x <= P % x .AND. P % x <= S % P1 % x ) THEN
          status = 1_iWP_sph
          RETURN
       END IF
       IF ( S % P0 % x >= P % x .AND. P % x >= S % P1 % x ) THEN
          status = 1_iWP_sph
          RETURN
       END IF

    ELSE  ! S is vertical, so test y coordinate

       IF ( S % P0 % y <= P % y .AND. P % y <= S % P1 % y ) THEN
          status = 1_iWP_sph
          RETURN
       END IF
       IF ( S % P0 % y >= P % y .AND. P % y >= S % P1 % y ) THEN
          status = 1_iWP_sph
         RETURN
       END IF

    END IF

    status = 0_iWP_sph

  END FUNCTION inSegment



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!                                                          !!
  !!                   intersection2D_2Segments               !!
  !!                                                          !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !                                                            !
  !     intersection2D_2Segments: find the 2D intersection of  !
  !                               two finite segments.         !
  !     Input : two finite segments S1 and S2                  !
  !     Output: IP0 = the intersect point (when it exists)     !
  !             IP1 = endpoint intersect segment[IP0, IP1]     !
  !                   (when it exists)                         !
  !     Return: 0 = disjoint (no intersection)                 !
  !             1 = intersection in the unique point IP0       !
  !             2 = overlap in segment from IP0 to IP1         !
  !                                                            !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE intersection2D_2Segments( S1, S2, IP0, IP1, status )
    CLASS   ( SPHSegment )    , INTENT ( IN )     :: S1, S2
    TYPE    ( point )         , INTENT ( IN OUT ) :: IP0, IP1
    INTEGER ( KIND = iWP_sph ), INTENT ( OUT )    :: status
    TYPE    ( vector )                            :: u_, v_, w_, w2_
    REAL    ( KIND = WP_sph )                     :: D_
    REAL    ( KIND = WP_sph )                     :: du_, dv_, t_, t0_, t1_, sI_, tI_


    u_ = S1 % P1 - S1 % P0
    v_ = S2 % P1 - S2 % P0
    w_ = S1 % P0 - S2 % P0
    D_ = perpProduct( u_, v_ )

    !!   test if  they are parallel (includes either being a point)   !!

    IF ( ABS (D_) < SMALL_NUMBER ) THEN   ! S1 and S2 are parallel

       IF ( perpProduct(u_,w_) /= 0.0_WP_sph .OR. perpProduct(v_,w_) /= 0.0_WP_sph ) THEN
          status = 0_iWP_sph    ! they are NOT collinear
          RETURN
       END IF

       !!   they are collinear or degenerate   !!
       !!   check if they are degenerate  points   !!

       du_ = dotProduct( u_, u_ )
       dv_ = dotProduct( v_, v_ )

       IF ( du_ == 0.0_WP_sph .AND. dv_ == 0.0_WP_sph ) THEN  ! both segments are points
          IF ( S1 % P0 .NOTEQUAL. S2 % P0 ) THEN   ! they are distinct points
             status = 0_iWP_sph
             RETURN
          END IF
          IP0 = S1 % P0
          status = 1_iWP_sph   ! they are the same point
          RETURN
       END IF

       IF ( du_ == 0.0_WP_sph ) THEN   ! S1 is a single point
          IF ( inSegment( S1 % P0, S2 ) == 0_iWP_sph ) THEN   ! but is not in S2
             status = 0_iWP_sph
             RETURN
          END IF
          IP0 = S1 % P0
          status = 1_iWP_sph
       END IF

       IF (dv_ == 0.0_WP_sph) THEN   ! S2 is a single point
          IF ( inSegment( S2 % P0, S1 ) == 0_iWP_sph ) THEN   ! but is not in S1
             status = 0_iWP_sph
             RETURN
          END IF
          IP0 = S2 % P0
          status = 1_iWP_sph
       END IF

       !!   they are collinear segments - get  overlap (or not)   !

       w2_ = S1 % P1 - S2 % P0
       IF ( v_ % x /= 0.0_WP_sph ) THEN
          t0_ = w_ % x / v_ % x
          t1_ = w2_ % x / v_ % x
       ELSE
          t0_ = w_ % y / v_ % y
          t1_ = w2_ % y / v_ % y
       END IF

       IF ( t0_ > t1_ ) THEN   ! must have t0 smaller than t1
          t_ = t0_; t0_ = t1_; t1_ = t_  ! swap if not
       END IF

       IF ( t0_ > 1.0_WP_sph .OR. t1_ < 0.0_WP_sph ) THEN   ! must have t0 smaller than t1
          status = 0_iWP_sph   ! No overlap
          RETURN
       END IF

       IF (t0_ < 0.0_WP_sph ) THEN
          t0_ = 0.0_WP_sph
       ELSE
          CONTINUE
       END IF
       IF (t1_ > 1.0_WP_sph  ) THEN
          t1_ = 1.0_WP_sph
       ELSE
          CONTINUE
       END IF

       !!   they overlap in a valid subsegment   !!

       IP0 = S2 % P0 + t0_ * v_
       IP1 = S2 % P0 + t1_ * v_
       status = 2_iWP_sph
       RETURN

    END IF

    !!   the segments are skew and may intersect in a point   !!
    !!   get the intersect parameter for S1   !!

    sI_ = perpProduct( v_, w_ ) / D_
    IF ( sI_ < 0.0_WP_sph .OR. sI_ > 1.0_WP_sph ) THEN   ! No intersection with S1
       status = 0_iWP_sph
       RETURN
    END IF

    !!   get the intersect parameter for S2   !!

    tI_ = perpProduct( u_, w_ ) / D_
    IF ( tI_ < 0.0_WP_sph .OR. tI_ > 1.0_WP_sph ) THEN   ! No intersection with S2
       status = 0_iWP_sph
       RETURN
    END IF

    IP0 = S1 % P0 + sI_ * u_
    status = 1_iWP_sph


  END SUBROUTINE intersection2D_2Segments



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!                                                          !!
  !!                 intersection3D_SegmentPlane              !!
  !!                                                          !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !                                                            !
  !     intersection3D_SegmentPlane: find the 3D intersection  !
  !                                  a segment and a plane.    !
  !     Input : S = a segment, and                             !
  !     PL = a plane = {Point PLP;  Vector n;}                 !
  !     Output: IP = the intersect point (when it exists)      !
  !     Return: 0 = disjoint (no intersection)                 !
  !             1 = intersection in the unique point IP        !
  !             2 = the  segment lies in the plane             !
  !                                                            !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE intersection3D_SegmentPlane( S, PL, IP, status )
    TYPE    ( point )         , INTENT ( OUT ) :: IP
    INTEGER ( KIND = iWP_sph ), INTENT ( OUT ) :: status
    CLASS   ( SPHSegment )    , INTENT ( IN )  :: S
    CLASS   ( SPHPlane )      , INTENT ( IN )  :: Pl
    TYPE    ( vector )                         :: u_, w_
    REAL    ( KIND = WP_sph )                  :: D_, N_
    REAL    ( KIND = WP_sph )                  :: sIP_
!    REAL    ( KIND = WP_sph ), PARAMETER       :: SMALL_NUMBER = 0.00000001_WP_sph


    u_   = S % P1 - S % P0
    w_   = S % P0 - PL % PLP
    D_ = dotProduct( PL % n, u_ )
    N_ = -dotProduct( PL % n, w_ )

    IF ( ABS( D_ ) < SMALL_NUMBER ) THEN      !  segment is parallel to plane

       IF ( N_ == 0.0_WP_sph ) THEN           !  segment lies in plane
          status = 2
          RETURN
       ELSE
          status = 0                          !  no intercection
          RETURN
       END IF

    END IF

    !   They are not parallel
    !   Compute intersection parameter

    sIP_ = N_ / D_
    IF ( sIP_ < 0.0_WP_sph .OR. sIP_ > 1.0_WP_sph ) THEN
       status = 0                              !  no intercection
       RETURN
    END IF

    IP = S % P0 + u_ * sIP_                 !  compute intersection point
    status = 1                              !  intercection


  END SUBROUTINE intersection3D_SegmentPlane



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!                                                          !!
  !!                          isLeft                          !!
  !!                                                          !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !                                                            !
  !     isLeft: tests if a point is Left|On|Right              !
  !             of an infinite line.                           !
  !     Input : three points P0, P1, and P2                    !
  !     Return: >0 for P2 left of the line through P0 and P1   !
  !             =1 for P2 on the line through P0 and P1        !
  !             <0 for P2 left of the line through P0 and P1   !
  !                                                            !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  REAL ( KIND = WP_sph ) FUNCTION isLeft( P0, P1, P2 ) RESULT ( status )
    CLASS ( point ), INTENT ( IN )  :: P0, P1, P2
    status = ( P1 % x - P0 % x ) * ( P2 % y - P0 % y )  &
        -    ( P2 % x - P0 % x ) * ( P1 % y - P0 % y )

  END FUNCTION isLeft



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!                                                          !!
  !!              crossingNumber_PointPolygon                 !!
  !!                                                          !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !                                                            !
  !     crossingNumber_PointPolygon: crossing number test      !
  !                                  for a point in a polygon. !
  !     Input : P = a point,                                   !
  !     Input : V() = vertex points of a polygon V(n+1)        !
  !             with V[n]=V[0]                                 !
  !     Return:  0 = outside, 1 = inside                       !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  INTEGER ( KIND = iWP_sph ) FUNCTION crossingNumber_PointPolygon( P, V, n ) RESULT ( cn )
    CLASS   ( point )         , INTENT ( IN )                       :: P
    CLASS   ( point )         , INTENT ( IN )    , DIMENSION ( : )  :: V
    INTEGER ( KIND = iWP_sph ), INTENT ( IN OUT )                   :: n
    INTEGER ( KIND = iWP_sph )                                      :: i_
    REAL    ( KIND = WP_sph )                                       :: vt_


    cn = 0_iWP_sph   ! the crossing number counter

    !!   Loop through all edges of the polygon   !!

    DO i_ = 0, n-1   !   edge from V(i) to V(i+1)

       IF (                                                            &
           ( (V(i_) % y <= P % y) .AND. (V(i_+1) % y > P % y) ) .OR.   &          !   an upward crossing
           ( (V(i_) % y > P % y)  .AND. (V(i_+1) % y <= P % y))        &          !   a downward crossing
           ) THEN
          vt_ = ( P % y - V(i_) % y ) / ( V(i_+1) % y - V(i_) % y )
          IF ( P % x < V(i_) % x + vt_ * (V(i_+1) % x - V(i_) % x) ) cn = cn + 1  ! a valid crossing of y=P%x right og P % x
       END IF

    END DO

    cn = MOD ( cn, 2_iWP_sph )   ! 0 if even (out), and 1 if odd (in)


  END FUNCTION crossingNumber_PointPolygon



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!                                                          !!
  !!               windingNumber_PointPolygon                 !!
  !!                                                          !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !                                                            !
  !     windingNumber_PointPolygon: winding number test        !
  !                                 for a point in a polygon.  !
  !     Input : P = a point,                                   !
  !     Input : V() = vertex points of a polygon V(n+1)        !
  !             with V[n]=V[0]                                 !
  !     Return: wn = the winding number (=0 only when P        !
  !             is outside)                                    !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  INTEGER ( KIND = iWP_sph ) FUNCTION windingNumber_PointPolygon( P, V, n ) RESULT ( wn )
    CLASS   ( point )         , INTENT ( IN )                       :: P
    CLASS   ( point )         , INTENT ( IN )    , DIMENSION ( : )  :: V
    INTEGER ( KIND = iWP_sph ), INTENT ( IN OUT )                   :: n
    INTEGER ( KIND = iWP_sph )                                      :: i_


    wn = 0_iWP_sph   ! the winding number counter

    !!   Loop through all edges of the polygon   !!

    DO i_ = 0, n-1   !   edge from V(i) to V(i+1)

       IF ( V(i_) % y <= P % y ) THEN              !   start y <= P % y
          IF ( V(i_+1) % y > P % y ) THEN          !   an upward crossing
             IF ( isLeft( V(i_), V(i_+1), P ) > 0.0_WP_sph ) wn = wn + 1
          END IF
       ELSE   ! start y > P % y (no test needed )
          IF ( V(i_+1) % y <= P % y ) THEN          !   a downward crossing
             IF ( isLeft( V(i_), V(i_+1), P ) < 0.0_WP_sph ) wn = wn - 1  ! P right of edge have a valid down intersection
          END IF
       END IF

    END DO


  END FUNCTION windingNumber_PointPolygon


#if add_

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!                                                          !!
!!             intersection2D_SegmentPolygon                !!
!!                                                          !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                            !
!     intersection2D_SegmentPolygon: intersect a 2D segment  !
!                                    with a convex polygon   !
!     Input: S = 2D segment to intersect with the            !
!            convex polygon V()                              !
!            n = number of 2D points in the polygon          !
!            V() = array of n+1 vertex points with           !
!            V(n) = V(0)                                     !
!     Note: The polygon MUST be convex and                   !
!           have vertices oriented counterclockwise (ccw).   !
!           This code does not check for and verify these    !
!           conditions.                                      !
!     Output: IS = the intersection segment (when it exists) !
!     Return: 0 = no intersection                            !
!             1 = a valid intersection segment exists        !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE intersection2D_SegmentPolygon( S, V, n, IS, status )
CLASS   ( Ssegment )      , INTENT ( IN )                       :: S
CLASS   ( Ssegment )      , INTENT ( IN OUT )                   :: IS
CLASS   ( point )         , INTENT ( IN )    , DIMENSION ( : )  :: V
INTEGER ( KIND = iWP_sph ), INTENT ( IN )                       :: n
INTEGER ( KIND = iWP_sph ), INTENT ( IN OUT )                   :: status
REAL ( KIND = WP_sph ) :: tE_, tL_, t_, N_, D_
TYPE ( vector ) :: dS_, e_
REAL    ( KIND = WP_sph ), PARAMETER       :: SMALL_NUMBER = 0.00000001_WP_sph
INTEGER ( KIND = iWP_sph )                                      :: i_

#if add_
IF ( S % P0 == S % P1 ) THEN   ! the segment is a single point
!!   test for inclusion of S % P0 in the polygon   !!
IS = S
!! ....
END IF
#endif

tE_ = 0.0_WP_sph  ! the maximum entering segment parameter
tL_ = 1.0_WP_sph  ! the minimum leaving segment parameter
dS_ = S % P1 - S % P0   ! the segment direction vector


DO i_ = 0, n-1

e_ = V(i_+1) - V(i_)   ! edge vector
N_ = perpProduct( e_, dS_ )  ! = -dot(ne, S%P0-V(i))
D_ =-perpProduct( e_, dS_ )  ! =dot(ne,dS)

IF ( ABS(D_) < SMALL_NUMBER ) THEN

IF (N_ < 0.0_WP_sph) THEN
status = 0
RETURN
ELSE
CONTINUE
END IF

END IF

#if add_
t_ = N_ / D_
IF ( D_ < 0.0_WP_sph) THEN
IF (t_>tE_)THEN
tE_ = t_
IF (tE_ > tL_)  THEN
status = 0
RETURN
END IF
END IF
END IF
#endif

END DO


! tE <= tL implies that there is a valid intersection subsegment

IS % P0 = S % P0 +  dS_ * tE_   ! = P(tE) = point where S enters polygon
IS % P1 = S % P0 +  dS_ * tL_   ! = P(tL) = point where S leaves polygon
status = 1

END SUBROUTINE intersection2D_SegmentPolygon

#endif



!! end of LINE/PLANE/TRACK type >
!!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!!


!! ----- !!
!! ----- !!


END MODULE SPH_pwi_module



!! --- !!
!! --- !!
!! --- !!


PROGRAM test
  USE SPH_pwi_module
  IMPLICIT NONE
  REAL    ( KIND = WP_sph )    , DIMENSION ( 3 ) :: a = (/ 1.0_WP_sph, 1.0_WP_sph, 1.0_WP_sph /)
  TYPE    ( point )                              :: p, p1, p2, p0
  TYPE    ( vector )                             :: v, v1, v2
  TYPE    ( SPHLine )                            :: L, LN1, LN2
  TYPE    ( SPHSegment )                         :: SG, SG1, SG2
  TYPE    ( SPHPlane )                           :: PLN
  TYPE    ( SPHPolygon )                         :: PLGN
  REAL    ( KIND = WP_sph )                      :: s
  INTEGER ( KIND = iWP_sph )                     :: i
  INTEGER ( KIND = iWP_sph )                     :: vec_test_
  INTEGER ( KIND = iWP_sph )                     :: pnt_test_
  INTEGER ( KIND = iWP_sph )                     :: ext_test_
  INTEGER ( KIND = iWP_sph )                     :: intersect_status_
  REAL    ( KIND =  WP_sph )                     :: dist_
!  LOGICAL ( KIND = lWP_sph )                     :: log_test_
  INTEGER ( KIND = iWP_sph )                     :: cn_
  INTEGER ( KIND = iWP_sph )                     :: wn_


  !! --- POINT ---  !!


  pnt_test_ = 0.0_iWP_sph

  PRINT '(/A)', ' *** Tests for POINT type *** '
  PRINT '(A)', 'Test #1 point constructor'
  p = point( 1.0_WP_sph, 1.0_WP_sph, 1.0_WP_sph )
  PRINT '(A,3(E15.8, 1x))', 'p= ', p
  IF ( ABS(p%x - 1.0_WP_sph) < SMALL_NUMBER .AND. ABS(p%y - 1.0_WP_sph) < SMALL_NUMBER .AND. &
       ABS(p%z - 1.0_WP_sph) < SMALL_NUMBER ) THEN
     PRINT '(A)', 'Test #1 passed ok <------------------------------------------'
     pnt_test_ = pnt_test_ + 1
  ELSE
     PRINT '(A)', 'Test #1 failed!'
  END IF

  !!

  PRINT '(/A)', 'Test #2 assignment = operator -> p1 = p'
  p1 = p
  PRINT '(A,3(E15.8, 1x))', 'p1= ', p1
  IF ( ABS(p1%x - 1.0_WP_sph) < SMALL_NUMBER .AND. ABS(p1%y - 1.0_WP_sph) < SMALL_NUMBER .AND. &
      ABS(p1%z - 1.0_WP_sph) < SMALL_NUMBER ) THEN
     PRINT '(A)', 'Test #2 passed ok <------------------------------------------'
     pnt_test_ = pnt_test_ + 1
  ELSE
     PRINT '(A)', 'Test #2 failed!'
  END IF

  !!

  PRINT '(/A)', 'Test #3 assignment = operator -> p2 = a'
  p2 = a
  PRINT '(A,3(E15.8, 1x))', 'p2= ', p2
  IF ( ABS(p2%x - 1.0_WP_sph) < SMALL_NUMBER .AND. ABS(p2%y - 1.0_WP_sph) < SMALL_NUMBER .AND. &
      ABS(p2%z - 1.0_WP_sph) < SMALL_NUMBER ) THEN
     PRINT '(A)', 'Test #3 passed ok <------------------------------------------'
     pnt_test_ = pnt_test_ + 1
  ELSE
     PRINT '(A)', 'Test #3 failed!'
  END IF

  !!

  PRINT '(/A)', 'Test #4 distance operator ->  p .D. p1'
  p  = point( 6.0_WP_sph, 6.0_WP_sph, 6.0_WP_sph )
  p1 = point( 4.0_WP_sph, 2.0_WP_sph, 1.0_WP_sph )
  s = p .D. p1
  PRINT '(A,1(E15.8, 1x))', 'p .D. p1= ', s
  IF ( ABS(s - 6.708203932499369_WP_sph) < SMALL_NUMBER  ) THEN
     PRINT '(A)', 'Test #4 passed ok <------------------------------------------'
     pnt_test_ = pnt_test_ + 1
  ELSE
     PRINT '(A)', 'Test #4 failed!'
  END IF

  !!

  PRINT '(/A)', 'Test #5 translate + operator ->  p + v'
  p = point( 1.0_WP_sph, 2.0_WP_sph, 3.0_WP_sph )
  v = vector( 1.0_WP_sph, 2.0_WP_sph, 3.0_WP_sph )
  p1 = p + v
  PRINT '(A,3(E15.8, 1x))', 'p + v= ', p1
  IF ( ABS(p1%x - 2.0_WP_sph) < SMALL_NUMBER .AND. ABS(p1%y - 4.0_WP_sph) < SMALL_NUMBER .AND. &
      ABS(p1%z - 6.0_WP_sph) < SMALL_NUMBER ) THEN
     PRINT '(A)', 'Test #5 passed ok <------------------------------------------'
     pnt_test_ = pnt_test_ + 1
  ELSE
     PRINT '(A)', 'Test #5 failed!'
  END IF

  !!

  PRINT '(/A)', 'Test #6 pointsDifference - operator ->  p1 - p2'
  p1 = point( 2.0_WP_sph, 3.0_WP_sph, 4.0_WP_sph )
  p2 = point( 1.0_WP_sph, 2.0_WP_sph, 3.0_WP_sph )
  v = p1 - p2
  PRINT '(A,3(E15.8, 1x))', 'p1 - p2= ', v
  IF ( ABS(v%x - 1.0_WP_sph) < SMALL_NUMBER .AND. ABS(v%y - 1.0_WP_sph) < SMALL_NUMBER .AND. &
      ABS(v%z - 1.0_WP_sph) < SMALL_NUMBER ) THEN
     PRINT '(A)', 'Test #6 passed ok <------------------------------------------'
     pnt_test_ = pnt_test_ + 1
  ELSE
     PRINT '(A)', 'Test #6 failed!'
  END IF

!p = a
!print*, p
!pt = point( 1.0_WP_sph, 1.0_WP_sph, 1.0_WP_sph )
!print*, p .EQUAL. pt, p .D. pt
!print*, p .EQUAL. pt, p .D. pt
!IF (pt .EQUAL. p) print*,'xxx'

  PRINT '(/A)', '*** SUMMARY of POINT tests ***'
  PRINT '(A,1x,I2,1x,A)', 'Number of tests passed: ', pnt_test_,   ' of 6 total.'
  PRINT '(A,1x,I2,1x,A)', 'Number of tests failed: ', -pnt_test_+6,' of 6 total.'


  !! --- VECTOR ---  !!


  vec_test_ = 0_iWP_sph

  PRINT '(/A)', ' *** Tests for VECTOR type *** '
  PRINT '(A)', 'Test #1 vector constructor'
  v = vector( 1.0_WP_sph, 1.0_WP_sph, 1.0_WP_sph )
  PRINT '(A,3(E15.8, 1x))', 'v= ', v
  IF ( ABS(v%x - 1.0_WP_sph) < SMALL_NUMBER .AND. ABS(v%y - 1.0_WP_sph) < SMALL_NUMBER .AND. &
       ABS(v%z - 1.0_WP_sph) < SMALL_NUMBER ) THEN
     PRINT '(A)', 'Test #1 passed ok <------------------------------------------'
     vec_test_ = vec_test_ + 1
  ELSE
     PRINT '(A)', 'Test #1 failed!'
  END IF

  !!

  PRINT '(/A)', 'Test #2 assignment = operator -> v1 = v'
  v1 = v
  PRINT '(A,3(E15.8, 1x))', 'v1= ', v1
  IF ( ABS(v1%x - 1.0_WP_sph) < SMALL_NUMBER .AND. ABS(v1%y - 1.0_WP_sph) < SMALL_NUMBER .AND. &
       ABS(v1%z - 1.0_WP_sph) < SMALL_NUMBER ) THEN
     PRINT '(A)', 'Test #2 passed ok <------------------------------------------'
     vec_test_ = vec_test_ + 1
  ELSE
     PRINT '(A)', 'Test #2 failed!'
  END IF

  !!

  PRINT '(/A)', 'Test #3 assignment = operator -> v2 = a'
  v2 = a
  PRINT '(A,3(E15.8, 1x))', 'v2= ', v2
  IF ( ABS(v2%x - 1.0_WP_sph) < SMALL_NUMBER .AND. ABS(v2%y - 1.0_WP_sph) < SMALL_NUMBER .AND. &
       ABS(v2%z - 1.0_WP_sph) < SMALL_NUMBER ) THEN
     PRINT '(A)', 'Test #3 passed ok <------------------------------------------'
     vec_test_ = vec_test_ + 1
  ELSE
     PRINT '(A)', 'Test #3 failed!'
  END IF

  !!

  PRINT '(/A)', 'Test #4 add + operator -> v = v1+v2'
  v = v1 + v2
  PRINT '(A,3(E15.8, 1x))', 'v= ', v
  IF ( ABS(v%x - 2.0_WP_sph) < SMALL_NUMBER .AND. ABS(v%y - 2.0_WP_sph) < SMALL_NUMBER .AND. &
       ABS(v%z - 2.0_WP_sph) < SMALL_NUMBER ) THEN
     PRINT '(A)', 'Test #4 passed ok <------------------------------------------'
     vec_test_ = vec_test_ + 1
  ELSE
     PRINT '(A)', 'Test #4 failed!'
  END IF

  !!

  PRINT '(/A)', 'Test #5 minus + operator -> v = v1-v2'
  v = v1 - v2
  PRINT '(A,3(E15.8, 1x))', 'v= ', v
  IF ( ABS(v%x - 0.0_WP_sph) < SMALL_NUMBER .AND. ABS(v%y - 0.0_WP_sph) < SMALL_NUMBER .AND. &
       ABS(v%z - 0.0_WP_sph) < SMALL_NUMBER ) THEN
     PRINT '(A)', 'Test #5 passed ok <------------------------------------------'
     vec_test_ = vec_test_ + 1
  ELSE
     PRINT '(A)', 'Test #5 failed!'
  END IF

  !!

  PRINT '(/A)', 'Test #6 scale * operator -> v = v1*s, (s=10.0)'
  s = 10.0_WP_sph
  v = v1 * s
  PRINT '(A,3(E15.8, 1x))', 'v= ', v
  IF ( ABS(v%x - 10.0_WP_sph) < SMALL_NUMBER .AND. ABS(v%y - 10.0_WP_sph) < SMALL_NUMBER .AND. &
       ABS(v%z - 10.0_WP_sph) < SMALL_NUMBER ) THEN
     PRINT '(A)', 'Test #6 passed ok <------------------------------------------'
     vec_test_ = vec_test_ + 1
  ELSE
     PRINT '(A)', 'Test #6 failed!'
  END IF

  !!

  PRINT '(/A)', 'Test #7 scale * operator -> v = s*v1, (s=20.0)'
  s = 20.0_WP_sph
  v = s * v1
  PRINT '(A,3(E15.8, 1x))', 'v= ', v
  IF ( ABS(v%x - 20.0_WP_sph) < SMALL_NUMBER .AND. ABS(v%y - 20.0_WP_sph) < SMALL_NUMBER .AND. &
       ABS(v%z - 20.0_WP_sph) < SMALL_NUMBER ) THEN
     PRINT '(A)', 'Test #7 passed ok <------------------------------------------'
     vec_test_ = vec_test_ + 1
  ELSE
     PRINT '(A)', 'Test #7 failed!'
  END IF

  !!

  PRINT '(/A)', 'Test #8 cross * operator -> v = v1xv2'
  v = v1 * v2
  PRINT '(A,3(E15.8, 1x))', 'v1xv2= ', v
  IF ( ABS(v%x - 0.0_WP_sph) < SMALL_NUMBER .AND. ABS(v%y - 0.0_WP_sph) < SMALL_NUMBER .AND. &
       ABS(v%z - 0.0_WP_sph) < SMALL_NUMBER ) THEN
     PRINT '(A)', 'Test #8 passed ok <------------------------------------------'
     vec_test_ = vec_test_ + 1
  ELSE
     PRINT '(A)', 'Test #8 failed!'
  END IF

!!

  PRINT '(/A)', 'Test #9 dotProduct operator -> v = v1.v2'
  s = v1 .DOT. v2
  PRINT '(A,1(E15.8, 1x))', 'v1.v2= ', s
  IF ( ABS( s-3.0_WP_sph ) < SMALL_NUMBER ) THEN
     PRINT '(A)', 'Test #9 passed ok <------------------------------------------'
     vec_test_ = vec_test_ + 1
  ELSE
     PRINT '(A)', 'Test #9 failed!'
  END IF

!!

!PRINT '(A)', 'Test perpProduct operator -> v = v1_|_v2'
!s = v1 .PERP. v2
!PRINT '(A,1(E15.8, 1x)/)', 'v1_|_v2= ', s

!!

  PRINT '(/A)', 'Test #10 dotProduct operator -> v = v1.v2'
  v1 = vector( 1.0_WP_sph, 3.0_WP_sph, -4.0_WP_sph )
  v2 = vector( 3.0_WP_sph, -5.0_WP_sph, 2.0_WP_sph )
  s = v1 .DOT. v2
  PRINT '(A,1(E15.8, 1x))', 'v1.v2= ', s
  IF ( ABS( s+20.0_WP_sph ) < SMALL_NUMBER ) THEN
     PRINT '(A)', 'Test #10 passed ok <------------------------------------------'
     vec_test_ = vec_test_ + 1
  ELSE
     PRINT '(A)', 'Test #10 failed!'
  END IF

  !!

  PRINT '(/A)', 'Test #11 cross * operator -> v = v1*v2'
  v1 = vector( 4.0_WP_sph, 9.0_WP_sph, 2.0_WP_sph )
  v2 = vector( 3.0_WP_sph, -3.0_WP_sph, 1.0_WP_sph )
  v = v1 * v2
  PRINT '(A,3(E15.8, 1x))', 'v1xv2= ', v
  IF ( ABS(v%x - 15.0_WP_sph) < SMALL_NUMBER .AND. ABS(v%y - 2.0_WP_sph) < SMALL_NUMBER .AND. &
      ABS(v%z +39.0_WP_sph) < SMALL_NUMBER ) THEN
     PRINT '(A)', 'Test #11 passed ok <------------------------------------------'
     vec_test_ = vec_test_ + 1
  ELSE
     PRINT '(A)', 'Test #11 failed!'
  END IF

  !!

  PRINT '(/A)', 'Test #12 norm operator -> norm(v)'
  v = vector( 3.0_WP_sph, 4.0_WP_sph, 0.0_WP_sph )
  s = norm( v )
  PRINT '(A,1(E15.8, 1x))', 'norm(v)= ', s
  IF ( ABS( s-5.0_WP_sph ) < SMALL_NUMBER ) THEN
     PRINT '(A)', 'Test #12 passed ok <------------------------------------------'
     vec_test_ = vec_test_ + 1
  ELSE
     PRINT '(A)', 'Test #12 failed!'
  END IF

  !!

  PRINT '(/A)', 'Test #13 distance operator -> distance(v1,v2)'
  v1 = vector( 6.0_WP_sph, 8.0_WP_sph, 0.0_WP_sph )
  v2 = vector( 3.0_WP_sph, 4.0_WP_sph, 0.0_WP_sph )
  s = v1 .D. v2
  PRINT '(A,1(E15.8, 1x))', 'distance(v1-v2)= ', s
  IF ( ABS( s-5.0_WP_sph ) < SMALL_NUMBER )  THEN
     PRINT '(A)', 'Test #13 passed ok <------------------------------------------'
     vec_test_ = vec_test_ + 1
  ELSE
     PRINT '(A)', 'Test #13 failed!'
  END IF

  PRINT '(/A)', '*** SUMMARY of VECTOR tests ***'
  PRINT '(A,1x,I2,1x,A)', 'Number of tests passed: ', vec_test_,   ' of 13 total.'
  PRINT '(A,1x,I2,1x,A)', 'Number of tests failed: ', -vec_test_+13,' of 13 total.'

!!

  ext_test_ = 0

  PRINT '(/A)', ' *** Tests for "distance_PointToLine *** '
  PRINT '(A)', 'Test #1 distance_PointToLine()'
  P  = point( 1.0_WP_sph, -4.0_WP_sph, 0.0_WP_sph )
  P0 = point( 0.0_WP_sph, -13.0_WP_sph / 8.0_WP_sph, 0.0_WP_sph )
  P1 = point( 1.0_WP_sph, -7.0_WP_sph / 8.0_WP_sph, 0.0_WP_sph )
  L %  P0 = point( 0.0_WP_sph, -13.0_WP_sph / 8.0_WP_sph, 0.0_WP_sph )
  L %  P1 = point( 1.0_WP_sph,  -7.0_WP_sph / 8.0_WP_sph, 0.0_WP_sph )
  s = distance_PointToLine( P, L )   !! 2.5
  PRINT '(A,1(E15.8, 1x))', 'distance_PointToLine= ', s
  IF ( ABS( s-2.5_WP_sph ) < SMALL_NUMBER )  THEN
     PRINT '(A)', 'Test #1 passed ok <------------------------------------------'
     ext_test_ = ext_test_ + 1
  ELSE
    PRINT '(A)', 'Test #1 failed!'
  END IF
!print*,s,distance_PointToSegment( P, SG )

  !!

!PRINT '(/A)', 'Test #2 distance_PointToLine'
!P  = point( -3.0_WP_sph, 7.0_WP_sph, 0.0_WP_sph )
!P0 = point( 0.0_WP_sph, 2.0_WP_sph, 0.0_WP_sph )
!P1 = point( 1.0_WP_sph, 16.0_WP_sph / 5.0_WP_sph, 0.0_WP_sph )
!L = SPHLine( P0, P1 )
!s = distance_PointToLine( P, L )   !! 5.506
!PRINT '(A,1(E15.8, 1x))', 'distance_PointToLine= ', s
!IF ( ABS( s-5.506_WP_sph ) < SMALL_NUMBER )  THEN
!PRINT '(A)', 'Test #2 passed ok <------------------------------------------'
!ext_test_ = ext_test_ + 1
!ELSE
!PRINT '(A)', 'Test #2 failed!'
!END IF
!print*,s,distance_PointToSegment( P, SG )

  !!

  PRINT '(/A)', 'Test #2 distance_PointToLine()'
  P  = point( 0.0_WP_sph, 2.0_WP_sph, 3.0_WP_sph )
  P0 = point( 5.0_WP_sph, 2.0_WP_sph, 1.0_WP_sph )
  P1 = point( 3.0_WP_sph, 1.0_WP_sph, -1.0_WP_sph )
  L = SPHLine( P0, P1 )
  s = distance_PointToLine( P, L )   !! 5
  PRINT '(A,1(E15.8, 1x))', 'distance_PointToLine= ', s
  IF ( ABS( s-5.0_WP_sph ) < SMALL_NUMBER )  THEN
     PRINT '(A)', 'Test #2 passed ok <------------------------------------------'
     ext_test_ = ext_test_ + 1
  ELSE
     PRINT '(A)', 'Test #2 failed!'
  END IF
!print*,s,distance_PointToSegment( P, SG )

  !!

  PRINT '(/A)', 'Test #3 distance_PointToLine()'
  P = point( 0.0_WP_sph, 2.0_WP_sph, 0.0_WP_sph )
  P0 = point( 4.0_WP_sph, 2.0_WP_sph, 1.0_WP_sph )
  P1 = point( 3.0_WP_sph, 2.0_WP_sph, 1.0_WP_sph )
  L = SPHLine( P0, P1 )
  s = distance_PointToLine( P, L )   !! 1
  PRINT '(A,1(E15.8, 1x))', 'distance_PointToLine= ', s
  IF ( ABS( s-1.0_WP_sph ) < SMALL_NUMBER )  THEN
     PRINT '(A)', 'Test #3 passed ok <------------------------------------------'
     ext_test_ = ext_test_ + 1
  ELSE
     PRINT '(A)', 'Test #3 failed!'
  END IF
!print*,s,distance_PointToSegment( P, SG )

  !!

  PRINT '(/A)', 'Test #4 distance_PointToLine()'
  P = point( -0.707_WP_sph, 0.707_WP_sph, 0.0_WP_sph )
  P0 = point( 0.0_WP_sph, 0.0_WP_sph, 0.0_WP_sph )
  P1 = point( 1.0_WP_sph, 0.0_WP_sph, 0.0_WP_sph )
  L = SPHLine( P0, P1 )
  s = distance_PointToLine( P, L )   !! 0.707/1.
  PRINT '(A,1(E15.8, 1x))', 'distance_PointToLine= ', s
  IF ( ABS( s-0.707_WP_sph ) < SMALL_NUMBER )  THEN
     PRINT '(A)', 'Test #4 passed ok <------------------------------------------'
     ext_test_ = ext_test_ + 1
  ELSE
     PRINT '(A)', 'Test #4 failed!'
  END IF
!print*,s,distance_PointToSegment( P, SG )

  !!

  PRINT '(/A)', 'Test #5 distance_PointToLine()'
  P = point( 0.707_WP_sph, 0.707_WP_sph, 0.0_WP_sph )
  P0 = point( 0.0_WP_sph, 0.0_WP_sph, 0.0_WP_sph )
  P1 = point( 1.0_WP_sph, 0.0_WP_sph, 0.0_WP_sph )
  L = SPHLine( P0, P1 )
  s = distance_PointToLine( P, L )   !! 0.707/0.707
  PRINT '(A,1(E15.8, 1x))', 'distance_PointToLine= ', s
  IF ( ABS( s-0.707_WP_sph ) < SMALL_NUMBER )  THEN
     PRINT '(A)', 'Test #5 passed ok <------------------------------------------'
    ext_test_ = ext_test_ + 1
  ELSE
    PRINT '(A)', 'Test #5 failed!'
  END IF
!print*,s,distance_PointToSegment( P, SG )

  !!

  PRINT '(/A)', 'Test #6 distance_PointToPlane()'
  P  = point( 3.0_WP_sph, 4.0_WP_sph, 3.0_WP_sph )
  P1 = point( 0.0_WP_sph, 0.0_WP_sph, 0.0_WP_sph )
  PLN % PLP = point(  0.0_WP_sph, 0.0_WP_sph, 2.0_WP_sph )
  PLN % N   = vector( 7.0_WP_sph, 5.0_WP_sph,-3.0_WP_sph )
  CALL distance_PointToPlane( P, PLN, P1, dist_ )
  PRINT '(A,1(E15.8, 1x))', 'distance_PointToPlane= ', dist_
  PRINT '(A,3(E15.8, 1x))', 'base point P         = ', P1
  IF ( ABS( dist_- 38.0_WP_sph/SQRT(83.0_WP_sph) ) < SMALL_NUMBER )  THEN
     PRINT '(A)', 'Test #6 passed ok <------------------------------------------'
     ext_test_ = ext_test_ + 1
  ELSE
     PRINT '(A)', 'Test #6 failed!'
  END IF
!print*, P1
!print*, dist_

  !!

  PRINT '(/A)', 'Test #7 distance_PointToPlane()'
  P  = point( 5.0_WP_sph, 1.0_WP_sph, 2.0_WP_sph )
  P1 = point( 0.0_WP_sph, 0.0_WP_sph, 0.0_WP_sph )
  PLN % PLP = point(  0.0_WP_sph, 0.0_WP_sph, 3.0_WP_sph )
  PLN % N   = vector( 2.0_WP_sph, -2.0_WP_sph,-1.0_WP_sph )
  CALL distance_PointToPlane( P, PLN, P1, dist_ )
  PRINT '(A,1(E15.8, 1x))', 'distance_PointToPlane= ', dist_
  PRINT '(A,3(E15.8, 1x))', 'base point P         = ', P1
  IF ( ABS( dist_- 3.0_WP_sph ) < SMALL_NUMBER )  THEN
     PRINT '(A)', 'Test #7 passed ok <------------------------------------------'
     ext_test_ = ext_test_ + 1
  ELSE
     PRINT '(A)', 'Test #8 failed!'
  END IF

  !!

  PRINT '(/A)', 'Test #8 distance_PointToPlane()'
  P  = point( 0.0_WP_sph, 0.0_WP_sph, 0.0_WP_sph )
  P1 = point( 0.0_WP_sph, 0.0_WP_sph, 0.0_WP_sph )
  PLN % PLP = point(  1.0_WP_sph, -3.0_WP_sph, 0.0_WP_sph )
  PLN % N   = vector( 3.0_WP_sph, -1.0_WP_sph,-4.0_WP_sph )
  CALL distance_PointToPlane( P, PLN, P1, dist_ )
  PRINT '(A,1(E15.8, 1x))', 'distance_PointToPlane= ', dist_
  PRINT '(A,3(E15.8, 1x))', 'base point P         = ', P1
  IF ( ABS( dist_ - 6.0_WP_sph / SQRT( 26.0_WP_sph ) ) < SMALL_NUMBER )  THEN
     PRINT '(A)', 'Test #8 passed ok <------------------------------------------'
     ext_test_ = ext_test_ + 1
  ELSE
     PRINT '(A)', 'Test #8 failed!'
  END IF

  !!

  PRINT '(/A)', 'Test #9 distance3D_LineToLine()'
  P0  = point( 0.0_WP_sph, 6.0_WP_sph, 0.0_WP_sph )
  P1  = point( 1.0_WP_sph, 7.0_WP_sph, 0.0_WP_sph )
  LN1 = SPHLine( P0, P1 )
  P0  = point( 0.0_WP_sph, -2.0_WP_sph, 0.0_WP_sph )
  P1  = point( 1.0_WP_sph, -1.0_WP_sph, 0.0_WP_sph )
  LN2 = SPHLine( P0, P1 )
  dist_ = distance3D_LineToLine( LN1, LN2 )
  PRINT '(A,1(E15.8, 1x))', 'distance3D_LineToLine= ', dist_
  IF ( ABS( dist_- 8.0_WP_sph / SQRT (2.0_WP_sph) ) < SMALL_NUMBER )  THEN
     PRINT '(A)', 'Test #9 passed ok <------------------------------------------'
     ext_test_ = ext_test_ + 1
  ELSE
     PRINT '(A)', 'Test #9 failed!'
  END IF
!print*, dist_

  !!

  PRINT '(/A)', 'Test #10 distance3D_LineToLine()'
  P0  = point( 0.0_WP_sph, 34.0_WP_sph/8.0_WP_sph, 0.0_WP_sph )
  P1  = point( 1.0_WP_sph, 19.0_WP_sph/8.0_WP_sph, 0.0_WP_sph )
  LN1 = SPHLine( P0, P1 )
  P0  = point( 0.0_WP_sph, -31.0_WP_sph/8.0_WP_sph, 0.0_WP_sph )
  P1  = point( 1.0_WP_sph, -46.0_WP_sph/8.0_WP_sph, 0.0_WP_sph )
  LN2 = SPHLine( P0, P1 )
  dist_ = distance3D_LineToLine( LN1, LN2 )
  PRINT '(A,1(E15.8, 1x))', 'distance3D_LineToLine= ', dist_
  IF ( ABS( dist_- 65.0_WP_sph / SQRT (289.0_WP_sph) ) < SMALL_NUMBER )  THEN
     PRINT '(A)', 'Test #10 passed ok <------------------------------------------'
     ext_test_ = ext_test_ + 1
  ELSE
     PRINT '(A)', 'Test #10 failed!'
  END IF

  !!

  PRINT '(/A)', 'Test #11 distance3D_SegmentToSegment()'
  P0  = point( 0.0_WP_sph, 0.0_WP_sph, 0.0_WP_sph )
  P1  = point( 1.0_WP_sph, 0.0_WP_sph, 0.0_WP_sph )
  SG1 = SPHSegment( P0, P1 )
  P0  = point( 0.0_WP_sph, 1.0_WP_sph, 0.0_WP_sph )
  P1  = point( 1.0_WP_sph, 1.0_WP_sph, 0.0_WP_sph )
  SG2 = SPHSegment( P0, P1 )
  dist_ = distance3D_SegmentToSegment( SG1, SG2 )
  PRINT '(A,1(E15.8, 1x))', 'distance3D_SegmentToSegment= ', dist_
  IF ( ABS( dist_- 1.0_WP_sph ) < SMALL_NUMBER )  THEN
     PRINT '(A)', 'Test #11 passed ok <------------------------------------------'
     ext_test_ = ext_test_ + 1
  ELSE
     PRINT '(A)', 'Test #11 failed!'
  END IF
!print*, dist_

  !!

  PRINT '(/A)', 'Test #12 distance3D_SegmentToSegment()'
  P0 = point( 0.0_WP_sph, 0.0_WP_sph, 0.0_WP_sph )
  P1 = point( 1.0_WP_sph, 0.0_WP_sph, 0.0_WP_sph )
  SG1 = SPHSegment( P0, P1 )
  P0 = point( -SQRT(0.25_WP_sph*0.0001_WP_sph), 0.00001_WP_sph + 0.25_WP_sph*0.0001_WP_sph, 0.0_WP_sph )
  P1 = point(  SQRT(0.25_WP_sph*0.0001_WP_sph), 0.00001_WP_sph - 0.25_WP_sph*0.0001_WP_sph, 0.0_WP_sph )
  SG2 = SPHSegment( P0, P1 )
  dist_ = distance3D_SegmentToSegment( SG1, SG2 )
  PRINT '(A,1(E15.8, 1x))', 'distance3D_SegmentToSegment= ', dist_
  IF ( ABS( dist_ - 0.0_WP_sph ) < SMALL_NUMBER )  THEN
     PRINT '(A)', 'Test #12 passed ok <------------------------------------------'
     ext_test_ = ext_test_ + 1
  ELSE
     PRINT '(A)', 'Test #12 failed!'
  END IF
!print*, dist_

  !!

  PRINT '(/A)', 'Test #13 distance3D_SegmentToSegment()'
  P0 = point( 0.77998990099877119_WP_sph, 0.61192502360790968_WP_sph, -0.22703111823648214_WP_sph )
  P1 = point( 0.53215344529598951_WP_sph, 0.85724585503339767_WP_sph, -0.10102437809109688_WP_sph )
  SG1 = SPHSegment( P0, P1 )
  P0 = point( -0.21277333982288837_WP_sph, 0.35091548087075353_WP_sph, -0.49557160679250956_WP_sph  )
  P1 = point( 0.11881479667499661_WP_sph, 0.022494725417345762_WP_sph, -0.66426620958372951_WP_sph )
  SG2 = SPHSegment( P0, P1 )
  dist_ = distance3D_SegmentToSegment( SG1, SG2 )
  PRINT '(A,1(E15.8, 1x))', 'distance3D_SegmentToSegment= ', dist_
  IF ( ABS( dist_- 0.98292397116488739_WP_sph ) < SMALL_NUMBER )  THEN
     PRINT '(A)', 'Test #13 passed ok <------------------------------------------'
     ext_test_ = ext_test_ + 1
  ELSE
     PRINT '(A)', 'Test #13 failed!'
  END IF
!print*, dist_

  !!

  PRINT '(/A)', 'Test #14 intersection3D_SegmentPlane()'
  SG % P0 = point( -1.0_WP_sph, 5.0_WP_sph, 1.0_WP_sph )
  SG % P1 = point( -3.0_WP_sph, 5.0_WP_sph, 2.0_WP_sph )
  PLN % PLP = point(  0.0_WP_sph, 0.0_WP_sph, 2.0_WP_sph )
  PLN % N   = vector( 4.0_WP_sph, 0.0_WP_sph, 1.0_WP_sph )
  CALL  intersection3D_SegmentPlane( SG, PLN, P, intersect_status_ )
  PRINT '(A,L2)', 'distance_PointToPlane= ', intersect_status_
  PRINT '(A,3(E15.8, 1x))', 'base point P         = ', P
  IF ( ABS( P % x - 3.0_WP_sph/7.0_WP_sph ) < SMALL_NUMBER .AND. &
       ABS( P % y - 5.0_WP_sph ) < SMALL_NUMBER .AND. &
       ABS( P % z - 2.0_WP_sph/7.0_WP_sph ) < SMALL_NUMBER )  THEN
     PRINT '(A)', 'Test #14 passed ok <------------------------------------------'
     ext_test_ = ext_test_ + 1
  ELSE
     PRINT '(A)', 'Test #14 failed!'
  END IF
!print*, 'P=', P
!print*,3.0_WP_sph/7.0_WP_sph, 5.0,2.0_WP_sph/7.0_WP_sph
!print*, 'status', intersect_status_

  !!

  PRINT '(/A)', 'Test #15 intersection3D_SegmentPlane()'
  SG % P0 = point( 0.0_WP_sph, -64.0_WP_sph, -64.0_WP_sph )
  SG % P1 = point( 0.0_WP_sph, -64.0_WP_sph,  0.0_WP_sph )
  PLN % PLP = point(  0.0_WP_sph, 0.0_WP_sph, -24.0_WP_sph )
  PLN % N   = vector( 0.0_WP_sph, 0.0_WP_sph, 1.0_WP_sph )
  CALL  intersection3D_SegmentPlane( SG, PLN, P, intersect_status_ )
  PRINT '(A,L2)', 'distance_PointToPlane= ', intersect_status_
  PRINT '(A,3(E15.8, 1x))', 'base point P         = ', P
  IF ( ABS( P % x - 0.0_WP_sph ) < SMALL_NUMBER .AND. &
       ABS( P % y + 64.0_WP_sph ) < SMALL_NUMBER .AND. &
       ABS( P % z + 24.0_WP_sph ) < SMALL_NUMBER )  THEN
     PRINT '(A)', 'Test #15 passed ok <------------------------------------------'
     ext_test_ = ext_test_ + 1
  ELSE
     PRINT '(A)', 'Test #15 failed!'
  END IF

  !!

  PRINT '(/A)', 'Test #16 isLeft()'
  P0 = point( 0.0_WP_sph, 0.0_WP_sph, 0.0_WP_sph )
  P1 = point( 0.0_WP_sph, 1.0_WP_sph, 0.0_WP_sph )
  P2 = point( -1.0_WP_sph, 0.0_WP_sph, 0.0_WP_sph )
  s = isLeft( P0, P1, P2 )
  PRINT '(A,E15.8,1x)', 'isLeft= ', s   ! the point is at the left of the line >0
  IF ( s > 0.0_WP_sph )  THEN
     PRINT '(A)', 'Test #16 passed ok <------------------------------------------'
     ext_test_ = ext_test_ + 1
  ELSE
     PRINT '(A)', 'Test #16 failed!'
  END IF

  !!

  PRINT '(/A)', 'Test #17 isLeft()'
  P0 = point( 0.0_WP_sph, 0.0_WP_sph, 0.0_WP_sph )
  P1 = point( 0.0_WP_sph, 1.0_WP_sph, 0.0_WP_sph )
  P2 = point( 1.0_WP_sph, 0.0_WP_sph, 0.0_WP_sph )
  s = isLeft( P0, P1, P2 )
  PRINT '(A,E15.8,1x)', 'isLeft= ', s
  IF ( s < 0.0_WP_sph )  THEN    ! the point is at the right of the line <0
     PRINT '(A)', 'Test #17 passed ok <------------------------------------------'
     ext_test_ = ext_test_ + 1
  ELSE
     PRINT '(A)', 'Test #17 failed!'
  END IF

  !!

  PRINT '(/A)', 'Test #18 isLeft()'
  P0 = point( 0.0_WP_sph, 0.0_WP_sph, 0.0_WP_sph )
  P1 = point( 0.0_WP_sph, 1.0_WP_sph, 0.0_WP_sph )
  P2 = point( 0.0_WP_sph, 0.5_WP_sph, 0.0_WP_sph )
  s = isLeft( P0, P1, P2 )
  PRINT '(A,E15.8,1x)', 'isLeft= ', s
  IF ( s == 0.0_WP_sph )  THEN   ! the point is on the line =0
     PRINT '(A)', 'Test #18 passed ok <------------------------------------------'
     ext_test_ = ext_test_ + 1
  ELSE
     PRINT '(A)', 'Test #18 failed!'
  END IF

  !!

  PRINT '(/A)', 'Test #19 crossingNumber_PointPolygon()'
  P = point( 0.0_WP_sph, 0.0_WP_sph, 0.0_WP_sph )
  PLGN % nVertex = 4_iWP_sph
  ALLOCATE ( PLGN % V(0:4) )
  PLGN % V(0) = point( -1.0_WP_sph, -1.0_WP_sph, 0.0_WP_sph )
  PLGN % V(1) = point(  1.0_WP_sph, -1.0_WP_sph, 0.0_WP_sph )
  PLGN % V(2) = point(  1.0_WP_sph,  1.0_WP_sph, 0.0_WP_sph )
  PLGN % V(3) = point( -1.0_WP_sph,  1.0_WP_sph, 0.0_WP_sph )
  PLGN % V(4) = PLGN % V(0)
  cn_ = crossingNumber_PointPolygon( P, PLGN % V, PLGN % nVertex )
  PRINT '(A,I3,1x)', 'crossing number cn= ', cn_  !! 0 outside, 1 inside
  IF ( cn_ == 1_iWP_sph )  THEN   ! the point is inside of the polygon
     PRINT '(A)', 'Test #19 passed ok <------------------------------------------'
     ext_test_ = ext_test_ + 1
  ELSE
     PRINT '(A)', 'Test #19 failed!'
  END IF

  !!

  PRINT '(/A)', 'Test #20 crossingNumber_PointPolygon()'
  P = point( 2.0_WP_sph, 0.0_WP_sph, 0.0_WP_sph )
  PLGN % nVertex = 4_iWP_sph
  PLGN % V(0) = point( -1.0_WP_sph, -1.0_WP_sph, 0.0_WP_sph )
  PLGN % V(1) = point(  1.0_WP_sph, -1.0_WP_sph, 0.0_WP_sph )
  PLGN % V(2) = point(  1.0_WP_sph,  1.0_WP_sph, 0.0_WP_sph )
  PLGN % V(3) = point( -1.0_WP_sph,  1.0_WP_sph, 0.0_WP_sph )
  PLGN % V(4) = PLGN % V(0)
  cn_ = crossingNumber_PointPolygon( P, PLGN % V, PLGN % nVertex )
  PRINT '(A,I3,1x)', 'crossing number cn= ', cn_  !! 0 outside, 1 inside
  IF ( cn_ == 0_iWP_sph )  THEN   ! the point is outside of the polygon
     PRINT '(A)', 'Test #20 passed ok <------------------------------------------'
     ext_test_ = ext_test_ + 1
  ELSE
     PRINT '(A)', 'Test #20 failed!'
  END IF
  DEALLOCATE ( PLGN % V )

  !!

  PRINT '(/A)', 'Test #21 windingNumber_PointPolygon()'
  P = point( 0.0_WP_sph, 0.0_WP_sph, 0.0_WP_sph )
  PLGN % nVertex = 4_iWP_sph
  ALLOCATE ( PLGN % V(0:4) )
  PLGN % V(0) = point( -1.0_WP_sph, -1.0_WP_sph, 0.0_WP_sph )
  PLGN % V(1) = point(  1.0_WP_sph, -1.0_WP_sph, 0.0_WP_sph )
  PLGN % V(2) = point(  1.0_WP_sph,  1.0_WP_sph, 0.0_WP_sph )
  PLGN % V(3) = point( -1.0_WP_sph,  1.0_WP_sph, 0.0_WP_sph )
  PLGN % V(4) = PLGN % V(0)
  wn_ = windingNumber_PointPolygon( P, PLGN % V, PLGN % nVertex )
  PRINT '(A,I3,1x)', 'winding number wn= ', wn_  !! 0 outside
  IF ( wn_ /= 0_iWP_sph )  THEN   ! the point is inside of the polygon
     PRINT '(A)', 'Test #21 passed ok <------------------------------------------'
     ext_test_ = ext_test_ + 1
  ELSE
     PRINT '(A)', 'Test #21 failed!'
  END IF

  !!

  PRINT '(/A)', 'Test #22 windingNumber_PointPolygon()'
  P = point( 2.0_WP_sph, 0.0_WP_sph, 0.0_WP_sph )
  PLGN % nVertex = 4_iWP_sph
  PLGN % V(0) = point( -1.0_WP_sph, -1.0_WP_sph, 0.0_WP_sph )
  PLGN % V(1) = point(  1.0_WP_sph, -1.0_WP_sph, 0.0_WP_sph )
  PLGN % V(2) = point(  1.0_WP_sph,  1.0_WP_sph, 0.0_WP_sph )
  PLGN % V(3) = point( -1.0_WP_sph,  1.0_WP_sph, 0.0_WP_sph )
  PLGN % V(4) = PLGN % V(0)
  wn_ = crossingNumber_PointPolygon( P, PLGN % V, PLGN % nVertex )
  PRINT '(A,I3,1x)', 'winding number wn= ', wn_  !! 0 outside
  IF ( wn_ == 0_iWP_sph )  THEN   ! the point is outside of the polygon
     PRINT '(A)', 'Test #22 passed ok <------------------------------------------'
     ext_test_ = ext_test_ + 1
  ELSE
     PRINT '(A)', 'Test #22 failed!'
  END IF
  DEALLOCATE ( PLGN % V )

  !!

  PRINT '(/A)', 'Test #23 intersection2D_2Segments()'
  P0  = point( 0.0_WP_sph, 0.0_WP_sph, 0.0_WP_sph )
  P1  = point( 1.0_WP_sph, 1.0_WP_sph, 0.0_WP_sph )
  SG1 = SPHSegment( P0, P1 )
  P0  = point( 1.0_WP_sph, 1.0_WP_sph, 0.0_WP_sph )
  P1  = point( 1.0_WP_sph, 0.0_WP_sph, 0.0_WP_sph )
  SG2 = SPHSegment( P0, P1 )
  P0 = point( -1000.0_WP_sph, -1000.0_WP_sph, -1000.0_WP_sph )
  P1 = P0
  CALL intersection2D_2Segments( SG1, SG2, P0, P1, i )
  PRINT '(A,I3,1x)', 'intersection2D_2Segments= ', i  !! 0-no intersect, 1=intersect, 2=overlap
  IF ( i <= 1 ) PRINT '(A, 3(E15.8, 1x))', 'intersection point P         = ', P0
  IF ( i == 2 ) THEN
     PRINT '(A, 3(E15.8, 1x))', 'start point P0 overlap       = ', P0
     PRINT '(A, 3(E15.8, 1x))', 'end point P1 overlap         = ', P1
  END IF
  IF ( ABS( P0 % x - 1.0_WP_sph ) < SMALL_NUMBER .AND. &
       ABS( P0 % y - 1.0_WP_sph ) < SMALL_NUMBER .AND. &
       ABS( P0 % z - 0.0_WP_sph ) < SMALL_NUMBER )  THEN
     PRINT '(A)', 'Test #23 passed ok <------------------------------------------'
     ext_test_ = ext_test_ + 1
  ELSE
     PRINT '(A)', 'Test #23 failed!'
  END IF

  !!

  PRINT '(/A)', 'Test #24 intersection2D_2Segments()'
  P0  = point( 1.8_WP_sph, 2.1_WP_sph, 0.0_WP_sph )
  P1  = point( 0.8_WP_sph, 1.1_WP_sph, 0.0_WP_sph )
  SG1 = SPHSegment( P0, P1 )
  P0  = point( 1.0_WP_sph, 1.25_WP_sph, 0.0_WP_sph )
  P1  = point( 0.0_WP_sph, 1.25_WP_sph, 0.0_WP_sph )
  SG2 = SPHSegment( P0, P1 )
  P0 = point( -1000.0_WP_sph, -1000.0_WP_sph, -1000.0_WP_sph )
  P1 = P0
  CALL intersection2D_2Segments( SG1, SG2, P0, P1, i )
  PRINT '(A,I3,1x)', 'intersection2D_2Segments= ', i  !! 0-no intersect, 1=intersect, 2=overlap
  IF ( i <= 1 ) PRINT '(A, 3(E15.8, 1x))', 'intersection point P         = ', P0
     IF ( i == 2 ) THEN
        PRINT '(A, 3(E15.8, 1x))', 'start point P0 overlap       = ', P0
        PRINT '(A, 3(E15.8, 1x))', 'end point P1 overlap         = ', P1
     END IF
  IF ( ABS( P0 % x - 0.95_WP_sph ) < SMALL_NUMBER .AND. &
       ABS( P0 % y - 1.25_WP_sph ) < SMALL_NUMBER .AND. &
       ABS( P0 % z - 0.0_WP_sph ) < SMALL_NUMBER )  THEN
     PRINT '(A)', 'Test #24 passed ok <------------------------------------------'
     ext_test_ = ext_test_ + 1
  ELSE
     PRINT '(A)', 'Test #24 failed!'
  END IF

  !!

  PRINT '(/A)', 'Test #25 intersection2D_2Segments()'
  P0  = point( -1.0_WP_sph, 1.0_WP_sph, 0.0_WP_sph )
  P1  = point(  1.0_WP_sph, 1.0_WP_sph, 0.0_WP_sph )
  SG1 = SPHSegment( P0, P1 )
  P0  = point( 0.0_WP_sph, 1.0_WP_sph, 0.0_WP_sph )
  P1  = point( 2.0_WP_sph, 1.0_WP_sph, 0.0_WP_sph )
  SG2 = SPHSegment( P0, P1 )
  P0 = point( -1000.0_WP_sph, -1000.0_WP_sph, -1000.0_WP_sph )
  P1 = P0
  CALL intersection2D_2Segments( SG1, SG2, P0, P1, i )
  PRINT '(A,I3,1x)', 'intersection2D_2Segments= ', i  !! 0-no intersect, 1=intersect, 2=overlap
  IF ( i <= 1 ) PRINT '(A, 3(E15.8, 1x))', 'intersection point P         = ', P0
  IF ( i == 2 ) THEN
     PRINT '(A, 3(E15.8, 1x))', 'start point P0 overlap       = ', P0
     PRINT '(A, 3(E15.8, 1x))', 'end point P1 overlap         = ', P1
  END IF
  IF ( ( i == 2 ) ) THEN

     IF ( ABS( P0 % x - 0.0_WP_sph ) < SMALL_NUMBER .AND.      &
          ABS( P0 % y - 1.0_WP_sph ) < SMALL_NUMBER .AND.      &
          ABS( P0 % z - 0.0_WP_sph ) < SMALL_NUMBER )  THEN
        IF ( ABS( P1 % x - 1.0_WP_sph ) < SMALL_NUMBER .AND.   &
             ABS( P1 % y - 1.0_WP_sph ) < SMALL_NUMBER .AND.   &
             ABS( P1 % z - 0.0_WP_sph ) < SMALL_NUMBER )  THEN
           PRINT '(A)', 'Test #25 passed ok <------------------------------------------'
           ext_test_ = ext_test_ + 1
        ELSE
           PRINT '(A)', 'Test #25 failed!'
        END IF
     ELSE
        PRINT '(A)', 'Test #25 failed!'
     END IF

  END IF


  PRINT '(/A)', '*** SUMMARY of EXTENDED tests ***'
  PRINT '(A,1x,I2,1x,A)', 'Number of tests passed: ', ext_test_,   ' of 25 total.'
  PRINT '(A,1x,I2,1x,A)', 'Number of tests failed: ', -ext_test_+25,' of 25 total.'


  PRINT '(/A)', 'SPHPWI test suite terminated!'


END PROGRAM test
