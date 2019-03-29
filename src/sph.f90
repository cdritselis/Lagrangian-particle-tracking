!! SPH_module: Routines related to sph/pp type.
!! sph-0.9.0
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


MODULE SPHBasics_module
  USE Channelflow_mod
  IMPLICIT NONE

  INTEGER, PARAMETER :: iWP_sph   = SELECTED_INT_KIND( 12 )            ! Double integer
  INTEGER, PARAMETER :: lWP_sph   = SELECTED_INT_KIND( 12 )            ! Double logical
  INTEGER, PARAMETER :: WP_sph    = SELECTED_REAL_KIND( 15, 307 )      ! Double real
  INTEGER, PARAMETER :: WPinc_sph = SELECTED_REAL_KIND( 25, 1000 )     ! Quad real
!!
INTEGER ( KIND = iWP_sph ), PARAMETER         :: NOT_INITIALIZED = 0_iWP_sph
INTEGER ( KIND = iWP_sph ), PARAMETER         :: FIRST_TIME_INITIALIZED = 1_iWP_sph
INTEGER ( KIND = iWP_sph ), PARAMETER         :: ALREADY_INITIALIZED = 2_iWP_sph
INTEGER ( KIND = iWP_sph ), PARAMETER         :: UNDEFINED = -1_iWP_sph

  !
END MODULE SPHBasics_module



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!                                                          !!
!!                     MODULE SPH_module                    !!
!!                                                          !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
MODULE SPH_module
  USE SPHBasics_module
  IMPLICIT NONE


  !!   ---   !!


  TYPE, PUBLIC ::  sphArgList
    INTEGER                                :: argc
    CHARACTER( LEN = 50 ), ALLOCATABLE, DIMENSION (:)  :: argv
    INTEGER   ( KIND = iWP_sph )  , ALLOCATABLE  :: status
    REAL      ( KIND = WP_sph )   , ALLOCATABLE  :: Version_
    CHARACTER ( LEN = : )         , ALLOCATABLE  :: Author_
    REAL      ( KIND = WP_sph )   , ALLOCATABLE  :: densityRatio_
    REAL      ( KIND = WP_sph )   , ALLOCATABLE  :: diameter_
    REAL      ( KIND = WP_sph )   , ALLOCATABLE  :: responseTime_
    REAL      ( KIND = WP_sph )   , ALLOCATABLE  :: gravity_
    CHARACTER ( LEN = : )         , ALLOCATABLE  :: sphFileName_
    REAL      ( KIND = WP_sph )   , ALLOCATABLE  :: sphDt_
    REAL      ( KIND = WP_sph )   , ALLOCATABLE  :: sphTime_
    INTEGER   ( KIND = iWP_sph )  , ALLOCATABLE  :: sphTimeSteps_
    INTEGER   ( KIND = iWP_sph )  , ALLOCATABLE  :: sphTotalParticles_

  CONTAINS

    PROCEDURE, PASS, PUBLIC :: sph_ReadCLI
    PROCEDURE, PASS, PUBLIC :: sph_GetInt
    PROCEDURE, PASS, PUBLIC :: sph_GetReal
    PROCEDURE, PASS, PUBLIC :: sph_GetStr

  END TYPE sphArgList


  !!   ---   !!


  !   Abstract particle type   !

  TYPE, ABSTRACT :: Particle_abstract
     LOGICAL ( KIND = lWP_sph ) :: sph_user_defined = .FALSE.
   CONTAINS
     PROCEDURE, PASS, PUBLIC :: sph_is_defined
     PROCEDURE, PASS, PUBLIC :: sph_mark_as_defined
     PROCEDURE, PASS, PUBLIC :: sph_mark_as_undefined
  END TYPE Particle_abstract


  !   Basic particle properties   !


  TYPE, EXTENDS ( Particle_abstract ) :: Particle_properties
     REAL ( KIND = WP_sph ) :: densityRatio = 1000.0_WP
     REAL ( KIND = WP_sph ) :: diameter     = 0.01_WP
     REAL ( KIND = WP_sph ) :: mass         = 1.0_WP
     REAL ( KIND = WP_sph ) :: responseTime = 1.0_WP
     REAL ( KIND = WP_sph ) :: gravity      = -1.0_WP
  END TYPE Particle_properties


  !   Particle parameters   !


  TYPE, EXTENDS ( Particle_properties ) :: Particle_parameters
     REAL ( KIND = WP_sph ) :: Reynolds               = 1.0_WP
     REAL ( KIND = WP_sph ) :: normalizedResponseTime = 1.0_WP
     REAL ( KIND = WP_sph ) :: dragCoefficient        = 1.0_WP
  END TYPE Particle_parameters


  !   Simulation parameters   !


  TYPE, EXTENDS ( Particle_parameters ) :: Simulation_parameters
     INTEGER   ( KIND =  iWP_sph )  :: sphID             = 1_iWP     !  ID of the particle
     CHARACTER ( LEN = 2 )          :: sphFileName       = " "       !  Name of the output file
     INTEGER   ( KIND = iWP_sph )   :: sphTotalParticles = 1_iWP     !  Total particles
     REAL      ( KIND = WP_sph )    :: sphDt             = 0.01_WP   !  Time step
     REAL      ( KIND = WP_sph )    :: sphTime           = 1.0_iWP   !  Total time
     INTEGER   ( KIND = iWP_sph )   :: sphTimeSteps      = 1_iWP     !  Total time steps
  END TYPE Simulation_parameters


  !   Simulation state   !


  TYPE, EXTENDS( Simulation_parameters ) :: Simulation_state
     REAL ( KIND = WP_sph ), DIMENSION ( 3 )  :: sphPosition      = (/ 0.0_WP, 0.0_WP, 0.0_WP /)   ! Particle position
     REAL ( KIND = WP_sph ), DIMENSION ( 3 )  :: sphVelocity      = (/ 0.0_WP, 0.0_WP, 0.0_WP /)   ! Particle velocity
     REAL ( KIND = WP_sph ), DIMENSION ( 3 )  :: sphFluidVelocity = (/ 0.0_WP, 0.0_WP, 0.0_WP /)   ! Fluid velocity at particle position
     REAL ( KIND = WP_sph ), DIMENSION ( 3 )  :: sphForce         = (/ 0.0_WP, 0.0_WP, 0.0_WP /)   ! Particle forces
     REAL ( KIND = WP_sph ), DIMENSION ( 3 )  :: sph              = (/ 0.0_WP, 0.0_WP, 0.0_WP /)   ! temporary sph

   CONTAINS

     PROCEDURE, PASS, PUBLIC :: sph_placeParticle_point
     PROCEDURE, PASS, PUBLIC :: sph_MoveParticle
     PROCEDURE, PASS, PUBLIC :: sph_setParticleVelocity
     PROCEDURE, PASS, PUBLIC :: sph_setFluidVelocityAtParticle
     PROCEDURE, PASS, PUBLIC :: sph_setParticleID
     PROCEDURE, PASS, PUBLIC :: sph_saveParticleTrajectory
     PROCEDURE, PASS, PUBLIC :: sph_checkState
     PROCEDURE, PASS, PUBLIC :: sph_dampReflect
     PROCEDURE, PASS, PUBLIC :: sph_cyclic
     PROCEDURE, PASS, PUBLIC :: sph_TransferCLI
  END TYPE Simulation_state


CONTAINS


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!                                                          !!
  !!                        sph_ReadCLI                       !!
  !!                                                          !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE sph_ReadCLI( this )
    IMPLICIT NONE
    CLASS   ( sphArgList ), INTENT ( IN OUT ) :: this
    INTEGER                                   :: arg


    this % status = NOT_INITIALIZED
    this % argc   = command_argument_count()
    IF ( this % argc > 0 ) THEN
       IF ( .NOT. ALLOCATED ( this % argv ) ) THEN
          ALLOCATE ( this % argv( this % argc ) )
          DO arg = 1, this % argc
             CALL get_command_argument( arg, this % argv(arg) )
          END DO
          this % status = FIRST_TIME_INITIALIZED
       ELSE
          this % status = ALREADY_INITIALIZED
       END IF

    ELSE

       this % status = UNDEFINED

    END IF


  END SUBROUTINE sph_ReadCLI



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!                                                          !!
  !!                       Argumentlist_str                   !!
  !!                                                          !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  FUNCTION Argumentlist_str( this, label_argv )  RESULT ( rs )
    IMPLICIT NONE
    CLASS     ( sphArgList )       , INTENT ( IN )  :: this
    CHARACTER ( LEN = 50_iWP_sph  ), INTENT ( IN )  :: label_argv
    INTEGER                                         :: arg_, arg_test_
    CHARACTER ( LEN = 50_iWP_sph  )                 :: rs


    rs        =  "  "
    arg_      = 0
    arg_test_ = 0
    DO WHILE ( arg_test_ < this % argc )
       arg_test_ = arg_test_ + 1
       IF ( TRIM( label_argv ) == TRIM( this % argv( arg_test_ ) ) ) THEN
          arg_ = arg_test_
          EXIT
       END IF
    END DO


    argumentlist: SELECT CASE ( label_argv )

    CASE ( '-au', '--author' )
       PRINT*, 'Author: ', sph_GetAuthor() // ' (University of Thessaly)'
       STOP " sph terminated ..."

    CASE ( '-v', '--version' )
       PRINT*, 'Version: ', sph_GetVersion()
       STOP " sph terminated ..."

    CASE ( '-h', '--help' )
       CALL sph_printHelp()
       STOP  " sph terminated ..."

    CASE ( '-f', '--file' )
       arg_ = arg_ + 1
       rs = this % argv( arg_ )

    CASE DEFAULT
       PRINT '(A,A,/)', 'Unrecognized command-line option: ',  this % argv( arg_ )
       CALL sph_printHelp()
       STOP

  END SELECT argumentlist


END FUNCTION Argumentlist_str



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!                                                          !!
  !!                       Argumentlist_int                   !!
  !!                                                          !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  FUNCTION Argumentlist_int( this, label_argv )  RESULT ( ri )
    IMPLICIT NONE
    CLASS     ( sphArgList )       , INTENT ( IN )  :: this
    CHARACTER ( LEN = 50_iWP_sph  ), INTENT ( IN )  :: label_argv
    INTEGER                                         :: arg_, arg_test_
    INTEGER   ( KIND = iWP_sph  )                   :: ri


    ri        = 0_iWP_sph
    arg_      = 0
    arg_test_ = 0
    DO WHILE ( arg_test_ < this % argc )
       arg_test_ = arg_test_ + 1
       IF ( TRIM( label_argv ) == TRIM( this % argv( arg_test_ ) ) ) THEN
          arg_ = arg_test_
          EXIT
       END IF
    END DO


    argumentlist: SELECT CASE ( label_argv )

    CASE ( '-Nsteps', '--totalTimeStep', '-tp', '--totalParticles' )
       arg_ = arg_ + 1
       ri = sph_Str2int( this % argv( arg_ ) )

    CASE DEFAULT
       PRINT '(A,A,/)', 'Unrecognized command-line option: ',  this % argv( arg_ )
       CALL sph_printHelp()
       STOP

    END SELECT argumentlist


  END FUNCTION Argumentlist_int



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!                                                          !!
  !!                       Argumentlist_real                  !!
  !!                                                          !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  FUNCTION Argumentlist_real( this, label_argv )  RESULT ( rr )
    IMPLICIT NONE
    CLASS     ( sphArgList )       , INTENT ( IN )  :: this
    CHARACTER ( LEN = 50_iWP_sph  ), INTENT ( IN )  :: label_argv
    INTEGER                                         :: arg_, arg_test_
    REAL      ( KIND = WP_sph  )                    :: rr


    rr        = 0.0_iWP_sph
    arg_      = 0
    arg_test_ = 0
    DO WHILE ( arg_test_ < this % argc )
       arg_test_ = arg_test_ + 1
       IF ( TRIM( label_argv ) == TRIM( this % argv( arg_test_ ) ) ) THEN
          arg_ = arg_test_
          EXIT
       END IF
    END DO


    argumentlist: SELECT CASE ( label_argv )

    CASE ( '-d', '--diameter', '-S', '--densityRatio', '-dt', '--timeStep', '-t', '--totalTime' )
       arg_ = arg_ + 1
       rr = sph_Str2real( this % argv( arg_ ) )

    CASE DEFAULT
       PRINT '(A,A,/)', 'Unrecognized command-line option: ',  this % argv( arg_ )
       CALL sph_printHelp()
       STOP

    END SELECT argumentlist


  END FUNCTION Argumentlist_real



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!                                                          !!
  !!                        sph_GetReal                       !!
  !!                                                          !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  FUNCTION sph_GetReal( this, string_alias, string_whole, default_value, label) RESULT ( return_real )
    IMPLICIT NONE
    CLASS     ( sphArgList )   , INTENT ( IN )              :: this
    REAL      ( KIND = WP )    , INTENT ( IN )              :: default_value
    CHARACTER ( LEN = * )      , INTENT ( IN )              :: string_alias, string_whole
    CHARACTER ( LEN = * )      , INTENT ( IN ), OPTIONAL    :: label
    CHARACTER ( LEN = 50_iWP  )                             :: label_argv
    REAL      ( KIND = WP )                   , ALLOCATABLE :: return_real
    INTEGER                                                 :: arg_tmp_


    IF ( ALLOCATED ( return_real )  ) DEALLOCATE ( return_real )
    ALLOCATE ( return_real )
    return_real = default_value

    label_argv =  "NULL"
    arg_tmp_   = 0
    DO WHILE ( arg_tmp_ < this % argc )

       arg_tmp_ = arg_tmp_ + 1
       IF ( TRIM(string_alias) == TRIM(this % argv( arg_tmp_ ) ) ) THEN
          label_argv =  string_alias
          return_real = Argumentlist_real( this, label_argv )
          EXIT
       END IF

       IF ( string_whole == this%argv( arg_tmp_ ) ) THEN
          label_argv =  string_whole
          return_real = Argumentlist_real( this, label_argv )
          EXIT
       END IF

    END DO

    IF ( PRESENT (label) ) &
         PRINT*, label, return_real


  END FUNCTION sph_GetReal



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!                                                          !!
  !!                         sph_GetInt                       !!
  !!                                                          !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  FUNCTION sph_GetInt( this, string_alias, string_whole, default_value, label) RESULT ( return_int )
    IMPLICIT NONE
    CLASS     ( sphArgList )   , INTENT ( IN )                 :: this
    INTEGER   ( KIND = iWP )   , INTENT ( IN )                 :: default_value
    CHARACTER ( LEN = * )      , INTENT ( IN )                 :: string_alias, string_whole
    CHARACTER ( LEN = * )      , INTENT ( IN ) , OPTIONAL      :: label
    CHARACTER ( LEN = 50_iWP  )                                :: label_argv
    INTEGER   ( KIND = iWP )                   , ALLOCATABLE   :: return_int
    INTEGER                                                    :: arg_tmp_


    IF ( ALLOCATED ( return_int )  ) DEALLOCATE ( return_int )
    ALLOCATE ( return_int )
    return_int = default_value

    label_argv =  "NULL"
    arg_tmp_   = 0
    DO WHILE ( arg_tmp_ < this % argc )

       arg_tmp_ = arg_tmp_ + 1
       IF ( string_alias == this % argv( arg_tmp_ ) ) THEN
          label_argv = string_alias
          return_int = Argumentlist_int( this, label_argv )
          EXIT
       END IF

       IF ( string_whole == this % argv( arg_tmp_ ) ) THEN
          label_argv = string_whole
          return_int = Argumentlist_int( this, label_argv )
          EXIT
       END IF

    END DO

    IF ( PRESENT (label) ) &
         PRINT*, label, return_int


  END FUNCTION sph_GetInt



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!                                                          !!
  !!                        sph_GetStr                        !!
  !!                                                          !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  FUNCTION sph_GetStr( this, string_alias, string_whole, default_value, label) RESULT ( return_str )
    IMPLICIT NONE
    CLASS     ( sphArgList )                             :: this
    CHARACTER ( LEN = * )     , INTENT ( IN )            :: default_value
    CHARACTER ( LEN = * )     , INTENT ( IN )            :: string_alias, string_whole
    CHARACTER ( LEN = * )     , INTENT ( IN ), OPTIONAL  :: label
    CHARACTER ( LEN = 50_iWP )                           :: label_argv
    CHARACTER ( LEN = : )     , ALLOCATABLE              :: return_str
    INTEGER                                              :: arg_tmp_


    return_str = default_value

    label_argv =  "NULL"
    arg_tmp_   = 0
    DO WHILE ( arg_tmp_ < this%argc )

       arg_tmp_ = arg_tmp_ + 1
       IF ( string_alias == this%argv( arg_tmp_ ) ) THEN
          label_argv = string_alias
          return_str = Argumentlist_str( this, label_argv )
          EXIT
       END IF

       IF ( string_whole == this % argv( arg_tmp_ ) ) THEN
          label_argv = string_whole
          return_str = Argumentlist_str( this, label_argv )
          EXIT
       END IF

    END DO

    IF ( PRESENT (label) ) &
         PRINT*, label, return_str


  END FUNCTION sph_GetStr



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!                                                          !!
  !!                      sph_TransferCLI                     !!
  !!                                                          !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE sph_TransferCLI( thisss, this )
    IMPLICIT NONE
    CLASS ( sphArgList )      , INTENT ( IN OUT ) :: this
    CLASS ( Simulation_state ), INTENT ( IN OUT ) :: thisss


    IF ( ALLOCATED ( this % densityRatio_ ) ) THEN
       thisss % densityRatio = this % densityRatio_
       DEALLOCATE ( this % densityRatio_ )
    END IF

    IF ( ALLOCATED ( this % diameter_ ) ) THEN
       thisss % diameter = this % diameter_
       DEALLOCATE ( this % diameter_ )
    END IF

    IF ( ALLOCATED ( this % responseTime_ ) ) THEN
       thisss % responseTime = this % responseTime_
       DEALLOCATE ( this % responseTime_ )
    END IF

    IF ( ALLOCATED ( this % gravity_ ) ) THEN
       thisss % gravity = this % gravity_
       DEALLOCATE ( this % gravity_ )
    END IF

    IF ( ALLOCATED ( this % sphTotalParticles_ ) ) THEN
       thisss % sphTotalParticles = this % sphTotalParticles_
       DEALLOCATE ( this % sphTotalParticles_ )
    END IF

    IF ( ALLOCATED ( this % sphTimeSteps_ ) ) THEN
       thisss % sphTimeSteps = this % sphTimeSteps_
       DEALLOCATE ( this % sphTimeSteps_ )
    END IF

    IF ( ALLOCATED ( this % sphDt_ ) ) THEN
       thisss % sphDt = this % sphDt_
       DEALLOCATE ( this % sphDt_ )
    END IF

    IF ( ALLOCATED ( this % sphTime_ ) ) THEN
       thisss % sphTime = this % sphTime_
       DEALLOCATE ( this % sphTime_ )
    END IF


  END SUBROUTINE sph_TransferCLI


  !! --- !!


  !   +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++   !
  PURE SUBROUTINE sph_defaultParticleProperties( this )
    IMPLICIT NONE
    CLASS   ( Simulation_state ), INTENT ( IN OUT ) :: this

    this % densityRatio = 1000.0_WP_sph
    this % diameter     = 0.001_WP_sph
    this % mass         = this % densityRatio * 3.14159_WP_sph / 6.0_WP_sph * this % diameter**3
    this % responseTime = this % densityRatio * this % diameter**2 / 18.0_WP_sph
    this % gravity      = -1.0_WP
  END SUBROUTINE sph_defaultParticleProperties
  !   +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++   !
  PURE SUBROUTINE sph_defaultSimulationParameters( this )
    IMPLICIT NONE
    CLASS   ( Simulation_state ), INTENT ( IN OUT ) :: this

    this % sphFileName       = "pt"
    this % sphTotalParticles = 1_iWP_sph
    this % sphDt             = 0.01_WP_sph
    this % sphTime           = 1.0_WP
    this % sphTimeSteps      = INT ( this % sphTime / this % sphDt, KIND = iWP_sph )
  END SUBROUTINE sph_defaultSimulationParameters
  !   +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++   !


  !! --- !!


  SUBROUTINE sph_printHelp()
    IMPLICIT NONE


    PRINT '(A)', 'Usage: sph [--version] [--author] [--help] ...'
    PRINT '(A)', ''
    PRINT '(A)', 'Argument list for sph program.'
    PRINT '(A)', ''
    PRINT '(A)', 'Options:'
    PRINT '(A)', ''
    PRINT '(A)', '  -au,     --author               print author information and exit'
    PRINT '(A)', '  -v,      --version              print version information and exit'
    PRINT '(A)', '  -h,      --help                 print usage information and exit'
    PRINT '(A)', '  -d,      --diameter             set particle diameter'
    PRINT '(A)', '  -S,      --densityRatio         set density ratio'
    PRINT '(A)', '  -tp,     --response             set particle response time'
    PRINT '(A)', '  -f,      --file                 set name of output file'
    PRINT '(A)', '  -TP,     --totalParticles       set total particles'
    PRINT '(A)', '  -dt,     --timestep             set time step'
    PRINT '(A)', '  -t,      --totalTime            set total time'
    PRINT '(A)', '  -Nsteps, --totalTimeStep        set total time steps'


  END SUBROUTINE sph_printHelp


  !   +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++   !

  ! Get_version

  FUNCTION sph_GetVersion()  RESULT ( return_value )
    IMPLICIT NONE
    CHARACTER ( LEN = 6_iWP_sph )      , PARAMETER :: VERSION_STRING = "0.9"
    CHARACTER ( LEN( VERSION_STRING ) )            :: return_value
    return_value = VERSION_STRING
  END FUNCTION sph_GetVersion
  !   +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++   !

  ! Get_author

  FUNCTION sph_GetAuthor()  RESULT ( return_value )
    IMPLICIT NONE
    CHARACTER ( LEN( 'C.D. Dritselis' ) )  :: return_value
    return_value = 'C.D. Dritselis'
  END FUNCTION sph_GetAuthor
  !   +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++   !

  ! str2real

  FUNCTION sph_str2real( real_string )  RESULT ( return_r )
    IMPLICIT NONE
    CHARACTER ( LEN = * )      , INTENT ( IN )  :: real_string
    REAL      ( KIND = WP_sph )                 :: return_r
    READ ( real_string, * ) return_r
  END FUNCTION sph_str2real

  !   +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++   !

  ! str2int

  FUNCTION sph_str2int( int_string )  RESULT ( return_i )
    IMPLICIT NONE
    CHARACTER ( LEN = * )       , INTENT ( IN )  :: int_string
    INTEGER   ( KIND = iWP_sph )                 :: return_i
    READ ( int_string, * ) return_i
   END FUNCTION sph_str2int
  !   +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++   !


  !! ========= !!

  LOGICAL FUNCTION sph_is_defined( this )
    CLASS ( Particle_abstract ), INTENT ( IN OUT ) :: this
    sph_is_defined = this%sph_user_defined
  END FUNCTION sph_is_defined
  SUBROUTINE sph_mark_as_defined( this )
    CLASS ( Particle_abstract ), INTENT ( IN OUT ) :: this
    this%sph_user_defined = .TRUE.
  END SUBROUTINE sph_mark_as_defined
  SUBROUTINE sph_mark_as_undefined( this )
    CLASS ( Particle_abstract ), INTENT ( IN OUT ) :: this
    this%sph_user_defined = .FALSE.
  END SUBROUTINE sph_mark_as_undefined


  !! ========= !!


  !   +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++   !
  PURE INTEGER ( KIND = iWP_sph ) FUNCTION sph_particleID( this )
    IMPLICIT NONE
    CLASS   ( Simulation_state ), INTENT ( IN ) :: this
    sph_particleID = this % sphID
  END FUNCTION sph_particleID
  !   +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++   !
  PURE REAL ( KIND = iWP_sph ) FUNCTION sph_particleReynolds( this )
    IMPLICIT NONE
    CLASS   ( Simulation_state ), INTENT ( IN ) :: this
    sph_particleReynolds = SQRT(                                                       &
        ( this % sphVelocity(1) - this % sphFluidVelocity(1) )**2 +                    &
        ( this % sphVelocity(2) - this % sphFluidVelocity(2) )**2 +                    &
        ( this % sphVelocity(3) - this % sphFluidVelocity(3) )**2 ) * this % diameter
  END FUNCTION sph_particleReynolds
  !   +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++   !
  PURE REAL ( KIND = iWP_sph ) FUNCTION sph_particleResponseTime( this )
    IMPLICIT NONE
    CLASS   ( Simulation_state ), INTENT ( IN ) :: this
    sph_particleResponseTime = this % densityRatio * this % diameter**2 / 18.0_WP_sph
  END FUNCTION sph_particleResponseTime
  !   +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++   !
  PURE REAL ( KIND = iWP_sph ) FUNCTION sph_dragCoefficient( this )
    IMPLICIT NONE
    CLASS   ( Simulation_state ), INTENT ( IN ) :: this
    sph_dragCoefficient = 1.0_WP_sph + 0.15_WP_sph * sph_particleReynolds( this ) ** 0.687_WP_sph
  END FUNCTION sph_dragCoefficient
  !   +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++   !
  PURE REAL ( KIND = iWP_sph ) FUNCTION sph_particleNormalizedResponseTime( this )
    IMPLICIT NONE
    CLASS   ( Simulation_state ), INTENT ( IN ) :: this
    sph_particleNormalizedResponseTime = this % ResponseTime / sph_dragCoefficient ( this )
  END FUNCTION sph_particleNormalizedResponseTime
  !   +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++   !


  !! --- !!


  LOGICAL FUNCTION sph_setParticleID( this, id )
    IMPLICIT NONE
    CLASS   ( Simulation_state ), INTENT ( IN OUT ) :: this
    INTEGER ( KIND = iWP_sph )  , INTENT ( IN )      :: id

    this % sphID = id
    sph_setParticleID = .TRUE.

  END FUNCTION sph_setParticleID



  LOGICAL FUNCTION sph_placeParticle_point( this, point )
    IMPLICIT NONE
    CLASS   ( Simulation_state ), INTENT ( IN OUT )                :: this
    REAL    ( KIND = WP_sph )   , INTENT ( IN )    , DIMENSION (:) :: point

    this % sphPosition ( 1 ) = point( 1 )
    this % sphPosition ( 2 ) = point( 2 )
    IF ( SIZE(point) == 3 ) this % sphPosition ( 3 ) = point( 3 )
    sph_placeParticle_point = .TRUE.

  END FUNCTION sph_placeParticle_point

  !   +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++   !

  ELEMENTAL PURE REAL ( KIND = WP_sph ) FUNCTION sph_point( point )
    IMPLICIT NONE
    REAL ( KIND = WP_sph ), INTENT ( IN ) :: point

    sph_point = point

  END FUNCTION sph_point

  !   +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++   !

  REAL ( KIND = WP_sph ) FUNCTION sph_random( seed )
    IMPLICIT NONE
    REAL ( KIND = WP_sph ), INTENT ( IN ), OPTIONAL :: seed
    REAL ( KIND = WP_sph )                          :: rand_

    IF ( PRESENT( seed ) ) CALL RANDOM_SEED( )
    CALL RANDOM_NUMBER( rand_ )
    sph_random = rand_

  END FUNCTION sph_random


  !   +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++   !


  !LOGICAL FUNCTION sph_placeParticle_random( this, distance )
  !IMPLICIT NONE
  !CLASS   ( Simulation_state ), INTENT ( IN OUT )                :: this
  !REAL    ( KIND = WP_sph )       , INTENT ( IN )    , DIMENSION (:) :: distance
  !REAL ( KIND = WP_sph ) :: position1_, position2_, position3_


  !CALL RANDOM_SEED( )
  !CALL RANDOM_NUMBER( position1_ )
  !CALL RANDOM_NUMBER( position2_ )
  !CALL RANDOM_NUMBER( position3_ )
  !this % sphPosition ( 1 ) = position1_ * distance
  !this % sphPosition ( 2 ) = position2_
  !IF ( SIZE(point) == 3 ) this % sphPosition ( 3 ) = position3_

  !sph_placeParticle_random = .TRUE.
  !
  !END FUNCTION sph_placeParticle_random


  !! --- !!


  ELEMENTAL PURE REAL ( KIND = WP_sph ) FUNCTION sph_shift ( value, shift )
    REAL ( KIND = WP_sph ), INTENT ( IN ) :: value, shift
    sph_shift = value + shift
  END FUNCTION sph_shift

  ELEMENTAL PURE REAL ( KIND = WP_sph ) FUNCTION sph_scale ( value, scale )
    REAL ( KIND = WP_sph ), INTENT ( IN ) :: value, scale
    sph_scale = value * scale
  END FUNCTION sph_scale

  ELEMENTAL PURE REAL ( KIND = WP_sph ) FUNCTION sph_shiftAndScale ( value, shift, scale )
    REAL ( KIND = WP_sph ), INTENT ( IN ) :: value, shift, scale
    sph_shiftAndScale = (value + shift) * scale
  END FUNCTION sph_shiftAndScale

  ELEMENTAL PURE REAL ( KIND = WP_sph ) FUNCTION sph_scaleAndShift ( value, scale, shift )
    REAL ( KIND = WP_sph ), INTENT ( IN ) :: value, shift, scale
    sph_scaleAndShift = value * scale + shift
  END FUNCTION sph_scaleAndShift


  !! --- !!


  ELEMENTAL PURE REAL ( KIND = WP_sph ) FUNCTION sph_assign ( value )
    REAL ( KIND = WP_sph ), INTENT ( IN ) :: value
    sph_assign = value
  END FUNCTION sph_assign

  ELEMENTAL PURE REAL ( KIND = WP_sph ) FUNCTION sph_add ( value, add )
    REAL ( KIND = WP_sph ), INTENT ( IN ) :: value, add
    sph_add = value + add
  END FUNCTION sph_add

  ELEMENTAL PURE REAL ( KIND = WP_sph ) FUNCTION sph_multiply ( value, myltiply )
    REAL ( KIND = WP_sph ), INTENT ( IN ) :: value, myltiply
    sph_multiply = value * myltiply
  END FUNCTION sph_multiply


  !! --- !!


  LOGICAL FUNCTION sph_setParticleVelocity( this, velocity )
    IMPLICIT NONE
    CLASS   ( Simulation_state ), INTENT ( IN OUT )                :: this
    REAL    ( KIND = WP_sph )       , INTENT ( IN )    , DIMENSION (:) :: velocity

    this % sphVelocity ( 1 ) = velocity( 1 )
    this % sphVelocity ( 2 ) = velocity( 2 )
    IF ( SIZE(velocity) == 3 ) this % sphVelocity ( 3 ) = velocity( 3 )

    sph_setParticleVelocity = .TRUE.

  END FUNCTION sph_setParticleVelocity



  LOGICAL FUNCTION sph_setFluidVelocityAtParticle( this, fluidVelocity )
    IMPLICIT NONE
    CLASS   ( Simulation_state ), INTENT ( IN OUT )                :: this
    REAL    ( KIND = WP_sph )       , INTENT ( IN )    , DIMENSION (:) :: fluidVelocity

    this % sphFluidVelocity ( 1 ) = fluidVelocity( 1 )
    this % sphFluidVelocity ( 2 ) = fluidVelocity( 2 )
    IF ( SIZE(fluidVelocity) == 3 ) this % sphFluidVelocity ( 3 ) = fluidVelocity( 3 )

    sph_setFluidVelocityAtParticle = .TRUE.

  END FUNCTION sph_setFluidVelocityAtParticle



  !! ========= !!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!                                                          !!
  !!                      sph_MoveParticle                    !!
  !!                                                          !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  PURE SUBROUTINE sph_MoveParticle( this, whichGravity )
    IMPLICIT NONE
    CLASS   ( Simulation_state ), INTENT ( IN OUT )                 :: this
    INTEGER ( KIND = iWP_sph )  , INTENT ( IN )    , DIMENSION ( :) :: whichGravity


    !!   Calculate prerequisites   !!

    this % normalizedResponseTime = sph_particleNormalizedResponseTime( this )
    this % sph                    = this % sphVelocity

    !!   Advance velocity   !!

    this % sphVelocity = this % sphFluidVelocity +                                                                   &
         ( this % sphVelocity - this % sphFluidVelocity ) * EXP ( -this % sphDt / this % normalizedResponseTime ) +  &
         whichGravity * this % normalizedResponseTime * this % gravity *                                             &
         ( 1.0_WP_sph - EXP ( -this % sphDt / this % normalizedResponseTime ) )

    !!   Advance position  !!

    this % sphPosition = this % sphPosition + 0.5_WP_sph * ( this % sphVelocity + this % sph ) * this % sphDt


  END SUBROUTINE sph_MoveParticle



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!                                                          !!
  !!                       sph_dampReflect                    !!
  !!                                                          !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  LOGICAL FUNCTION sph_dampReflect ( this, damp, which, barrier )
    IMPLICIT NONE
    CLASS   ( Simulation_state ), INTENT ( IN OUT )  :: this
    REAL    ( KIND = WP_sph )   , INTENT ( IN )      :: damp, barrier
    INTEGER ( KIND = iWP_sph )  , INTENT ( IN )      :: which
    REAL    ( KIND = WP_sph )                        :: tBounce


    !   Scale back the distance travelled based on the time from collision

    tBounce = ( this % sphPosition ( which ) - barrier ) / this % sphVelocity ( which )
    this % sphPosition  =  this % sphPosition - ( 1.0_WP_sph - damp ) * tBounce * this % sphVelocity

    !   Reflect the position and the velocity

    this % sphPosition ( which ) = 2.0_WP_sph * barrier - this % sphPosition ( which )
    this % sphVelocity ( which ) = - this % sphVelocity ( which )

    !   Damp the velocities

    this % sphVelocity = damp * this % sphVelocity

    sph_dampReflect = .TRUE.


  END FUNCTION sph_dampReflect



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!                                                          !!
  !!                         sph_cyclic                       !!
  !!                                                          !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  LOGICAL FUNCTION sph_cyclic ( this, which, barrier )
    IMPLICIT NONE
    CLASS   ( Simulation_state ), INTENT ( IN OUT )  :: this
    REAL    ( KIND = WP_sph )   , INTENT ( IN )      :: barrier
    INTEGER ( KIND = iWP_sph )  , INTENT ( IN )      :: which


    !   Apply cyclic / periodic conditions

    this % sphPosition ( which ) = MOD( this % sphPosition( which ), barrier )

    sph_cyclic = .TRUE.


  END FUNCTION sph_cyclic



  PURE LOGICAL FUNCTION sph_checkState( this, xMinMax, yMinMax, zMinMax )
    IMPLICIT NONE
    CLASS   ( Simulation_state ), INTENT ( IN )                 :: this
    REAL    ( KIND = WP_sph )   , INTENT ( IN ), DIMENSION (:)  :: xMinMax, yMinMax, zMinMax

    sph_checkState = .TRUE.
    IF ( this % sphPosition ( 1 ) > xMinMax(2) .OR. this % sphPosition ( 1 ) < xMinMax(1) ) sph_checkState = .FALSE.
    IF ( this % sphPosition ( 2 ) > yMinMax(2) .OR. this % sphPosition ( 2 ) < yMinMax(1) ) sph_checkState = .FALSE.
    IF ( this % sphPosition ( 3 ) > zMinMax(2) .OR. this % sphPosition ( 3 ) < zMinMax(2) ) sph_checkState = .FALSE.

  END FUNCTION sph_checkState



  ELEMENTAL PURE LOGICAL FUNCTION sph_crossDomainBoundary( sphPosition, Min, Max )
    IMPLICIT NONE
    REAL      ( KIND = WP_sph ), INTENT ( IN ) :: Min, Max
    REAL      ( KIND = WP_sph ), INTENT ( IN ) :: sphPosition

    sph_crossDomainBoundary = .FALSE.
    IF ( sphPosition > Max .OR. sphPosition < Min ) sph_crossDomainBoundary = .TRUE.

  END FUNCTION sph_crossDomainBoundary


  !! ========= !!


  SUBROUTINE sph_saveParticleTrajectory( this )
    IMPLICIT NONE
    CLASS     ( Simulation_state ), INTENT ( IN ) :: this
    CHARACTER ( LEN = 8 )                         :: tmp_id_

    WRITE ( tmp_id_, '(i1)' ) this % sphID
    OPEN  ( 8888, FILE = TRIM( 'p' //  tmp_id_ )  // '.trajectory', FORM = 'formatted', STATUS = 'unknown', POSITION = 'append' )
    WRITE ( 8888,* ) this % sphPosition ( 1 ), this % sphPosition ( 2 ), this % sphPosition ( 3 )
    CLOSE ( 8888 )

  END SUBROUTINE sph_saveParticleTrajectory


  !! ========= !!
  !! ========= !!


END MODULE SPH_module






PROGRAM testSPH
  USE SPH_module
  IMPLICIT NONE
  TYPE ( Simulation_state ), DIMENSION ( 1 ) :: p
  TYPE ( sphArgList )                        :: aList
  INTEGER ( KIND = iWP_sph ) :: step_
  LOGICAL ( KIND = lWP_sph ) :: particleID_, placeParticle_, velocityParticle_, fluidVelocityAtParticle_
  LOGICAL ( KIND = lWP_sph ) :: particle_reflect_, particle_cyclic_


  !!   Read input for SPH   !!

  CALL aList % sph_ReadCLI()
  ALLOCATE ( aList % densityRatio_ )
  aList % densityRatio_      = aList % sph_GetReal( '-S'     , '--densityRatio'  , 1000.0_WP_sph , "Density ratio" )
  ALLOCATE ( aList % diameter_ )
  aList % diameter_          = aList % sph_GetReal( '-d'     , '--diameter'      , 0.01_WP_sph   , "Diameter" )
  ALLOCATE ( aList % ResponseTime_ )
  aList % ResponseTime_      =  &
      aList % sph_GetReal( '-tp'     , '--response'     , &
      aList % densityRatio_ * aList % diameter_**2 / 18.0_WP_sph , "Particle response time" )
  ALLOCATE ( aList % sphTotalParticles_ )
  aList % sphTotalParticles_ = aList % sph_GetInt ( '-TP'    , '--totalParticles', 1_iWP_sph     , "Total particles" )
  ALLOCATE ( aList % sphTimeSteps_ )
  aList % sphTimeSteps_      = aList % sph_GetInt ( '-Nsteps', '--totalTimeStep' , 10000_iWP_sph , "Total time steps" )
  ALLOCATE ( aList % sphDt_ )
  aList % sphDt_             = aList % sph_GetReal( '-dt'    , '--timestep'      , 0.01_WP_sph   , "Time step" )
  ALLOCATE ( aList % gravity_ )
  aList % gravity_           = aList % sph_GetReal( '-g'     , '--gravity'       , -1.0_WP_sph   , "Gravity" )
  CALL p(1) % sph_TransferCLI( aList )

  !!   Initialize SPH   !!

  particleID_              = p(1) % sph_setParticleID ( 1_iWP_sph )
  placeParticle_           = p(1) % sph_placeParticle_point ( [0.0_WP_sph, 1.0_WP_sph, 0.0_WP_sph] )
  velocityParticle_        = p(1) % sph_setParticleVelocity ( [1.0_WP_sph, 0.0_WP_sph, 0.0_WP_sph] )
  fluidVelocityAtParticle_ = p(1) % sph_setFluidVelocityAtParticle ( [1.0_WP_sph, 0.0_WP_sph, 0.0_WP_sph] )
#if SPH_WARN_ALL_ || SPH_DEBUG_ALL_
  PRINT '(A,3(E12.5,1X))', 'Particle position: ', p(1) % sphPosition
  PRINT '(A,3(E12.5,1X))', 'Particle velocity: ', p(1) % sphVelocity
#endif

  !!   Advance   !!

trajectory:  DO step_ = 1, p(1) % sphTimeSteps
     fluidVelocityAtParticle_ = p(1) % sph_setFluidVelocityAtParticle ( p(1) % sphVelocity )
     CALL p(1) % sph_MoveParticle( whichGravity = [0_iWP_sph, 1_iWP_sph, 0_iWP_sph] )
     PRINT '(A,3(E12.5,1X))', 'Particle position: ', p(1) % sphPosition
     PRINT '(A,3(E12.5,1X))', 'Particle velocity: ', p(1) % sphVelocity

     !!   Particle-wall   !!


     IF ( p(1) % sphPosition(2) < 0.0_WP_sph + 0.5_WP_sph * p(1) % Diameter ) THEN
        particle_reflect_ = &
            p(1) % sph_dampReflect ( damp = 1.0_WP_sph, which = 2_iWP_sph, barrier = 0.0_WP_sph + 0.5_WP_sph * p(1) % Diameter )
#if SPH_WARN_ALL_ || SPH_DEBUG_ALL_
        PRINT '(A,3(E12.5,1X))', 'Particle position after collision: ', p(1) % sphPosition
        PRINT '(A,3(E12.5,1X))', 'Particle velocity after collision: ', p(1) % sphVelocity
#endif
     END IF
     IF ( p(1) % sphPosition(2) > 1.0_WP_sph - 0.5_WP_sph * p(1) % Diameter ) then
        particle_reflect_ = &
            p(1) % sph_dampReflect ( damp = 1.0_WP_sph, which = 2_iWP_sph, barrier = 1.0_WP_sph - 0.5_WP_sph * p(1) % Diameter )
#if SPH_WARN_ALL_ || SPH_DEBUG_ALL_
        PRINT '(A,3(E12.5,1X))', 'Particle position after collision: ', p(1) % sphPosition
        PRINT '(A,3(E12.5,1X))', 'Particle velocity after collision: ', p(1) % sphVelocity
#endif
     END IF

     !!   Particle-periodic boundary   !!

     IF ( p(1) % sphPosition(1) < 0.0_WP_sph ) THEN
        particle_cyclic_ = p(1) % sph_cyclic( which = 1_iWP_sph, barrier = 0.0_WP_sph )
#if SPH_WARN_ALL_ || SPH_DEBUG_ALL_
        PRINT '(A,1(E12.5,1X))', 'Position after crossing periodic boundary: ', p(1) % sphPosition( 1 )
#endif
     END IF
     IF ( p(1) % sphPosition(1) > 5.0_WP_sph ) THEN
        particle_cyclic_ = p(1) % sph_cyclic( which = 1_iWP_sph, barrier = 5.0_WP_sph )
#if SPH_WARN_ALL_ || SPH_DEBUG_ALL_
        PRINT '(A,1(E12.5,1X))', 'Position after crossing periodic boundary: ', p(1) % sphPosition ( 1 )
#endif
     END IF
     IF ( p(1) % sphPosition(3) < 0.0_WP_sph ) THEN
        particle_cyclic_ = p(1) % sph_cyclic( which = 3_iWP_sph, barrier = 0.0_WP_sph )
#if SPH_WARN_ALL_ || SPH_DEBUG_ALL_
        PRINT '(A,1(E12.5,1X))', 'Position after crossing periodic boundary: ', p(1) % sphPosition ( 3 )
#endif
     END IF
     IF ( p(1) % sphPosition(3) > 5.0_WP_sph ) THEN
        particle_cyclic_ = p(1) % sph_cyclic( which = 3_iWP_sph, barrier = 5.0_WP_sph )
#if SPH_WARN_ALL_ || SPH_DEBUG_ALL_
        PRINT '(A,1(E12.5,1X))', 'Position after crossing periodic boundary: ', p(1) % sphPosition ( 3 )
#endif
     END IF

     !!   Save trajectory   !!

     CALL p(1) % sph_saveParticleTrajectory

  END DO trajectory


END PROGRAM testSPH
