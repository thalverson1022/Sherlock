!-------------------------------------------------------------------------------------------
Module Constants
  Implicit NONE
  Double Precision, Parameter :: PI = 3.1415926535897932385d0
  Double Precision, Parameter :: CONV = 219474.6313705d0
  Integer,          Parameter :: SIZE = 20                  !Parameter for the large arrays
  Integer,          Parameter :: Max_Len = 10000000              !Maximum length for MandN
  Double Precision, Parameter :: hbar = 1.0d0
  Integer,          Parameter :: Max_IT = 10000000 !Max number of interations
End Module Constants
!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------
Module UserInputs
  Implicit NONE  
  Integer          :: DMax         !Problem Dimensionality
  Integer          :: coeffs_len   !Total number of coulping/anharmonic constants  
  Double Precision :: EMax         !Truncation parameter 
  Double Precision :: shift        !m grid shift
  Integer		   :: mtype        !Switch for (1) half or (2) whole int grid
  Double Precision :: R_max        !Maximum Radius used for truncation 
  
  Character(len = 128) :: Out_FN       !Output file: eigenvalues
  Character(len = 128) :: PE_FN        !Input file: potential energy
  Character(len = 128) :: masses_FN    !Input file: masses
  Character(len = 128) :: UserIns_FN   !Input file: Emax and output filename
  Character(len = 128) :: maxes_FN     !Input file: Max grid values
  Character(len = 128) :: start_FN     !Input file: Max grid values
  
End Module UserInputs
!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------
Module ParameterArrays
  
  Implicit NONE
  
  Double Precision, Allocatable, Dimension(:)   :: Coeffs
  Integer,          Allocatable, Dimension(:,:) :: Powers
  
  Double Precision, Allocatable, Dimension(:)   :: Masses  

  Double Precision, Allocatable, Dimension(:,:) :: Xpoints !Array for zero power coordinates

  Integer,          Allocatable, Dimension(:,:) :: maxes     ! array for setting hard limits for each dof

  Double Precision, Allocatable, Dimension(:)   :: start_point 

  !Double Precision, Allocatable, Dimension(:,:) :: Qendpts   ! array used to record the final max
                                                             ! and min in each dof
  !Double Precision, Allocatable, Dimension(:)   :: Pendpts   ! array used to record the final max
                                                             ! and min in each dof
                                                               
End Module
!-------------------------------------------------------------------------------------------
Module Outputs
  Implicit NONE
  Integer										:: iMax         !Numbe of valid phase space points
  Double Precision, Allocatable, Dimension(:,:) :: Points   		!Array for Phase Space Tuncation 
End Module Outputs
!-------------------------------------------------------------------------------------------
                                                             