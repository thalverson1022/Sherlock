!===================================
!    N-D Search 2.0
!   Created by: Tom  Halverson
!   Date: 10-14-2014
!
!   - Improved Algorithm ordering as 
!     per Dr. P's suggestions
!
!===================================


    !+++++++++++++++++++++++++++++++++++++++++
    !++++++++++  Main Program Loop +++++++++++
    !+++++++++++++++++++++++++++++++++++++++++   
Program Main
  
  Implicit NONE

  Double Precision :: t1, t2, tot_time
  Character(len = 5) :: units

  Call cpu_time(t1)

  Print *
  Print *,'========================================================'
  Print *,'|         Sherlock v2.0                                |'
  Print *,'|         Created by: Thomas Halverson                 |'
  Print *,'|         Date: 10-14-14                                |'
  Print *,'|                                                      |'
  Print *,'|                                                      |'
  Print *,'========================================================'
  Print *

  Call ReadUserInputs()
  Call MakeXArrays()
  Call PST()
  Call WriteFile()

  Call cpu_time(t2)

  If(t2-t1 > 18000.0d0) then
    tot_time = (t2-t1)/(60.0d0*60.0d0)
    units = 'hrs.'
  Else if (t2-t1 > 300) then
    tot_time = (t2-t1)/(60.0d0)
    units = 'mins.'
  Else 
    tot_time = (t2-t1)
    units = 'sec.'
  End If

  Print 9999, tot_time, units
  Print *,'*****************************************************************'
  Print *
  Print *

  9999 Format ('  Total time taken: ', f10.5, 2x, a)
   
End Program Main



    !************************************************************!
    !********************    SUBROUTINES   **********************!
    !************************************************************!

  Subroutine ReadUserInputs()
    !----Modules----
    Use Constants
    Use UserInputs
    Use ParameterArrays
    !---------------
    
    Implicit NONE
    Integer :: i,j
    Character(len = 40) :: input_file
 
    Call GETARG(1, input_file)
    input_file = trim(input_file)
    
    Print *,'*****************************************************************'
    Print *,'        STEP 1:'
    Print '(1x,a,2x,a20)',"Reading data from: ", input_file
    Print *
    
!   ---- Read in files names (and directories) of user inputs ----  
    Open(UNIT = 15, FILE = input_file, STATUS='OLD',  ACTION ='Read')
    
      Read (15,'(A)') PE_FN
        PE_FN = Trim(PE_FN)
      Read (15,'(A)') masses_FN
        masses_FN = Trim(masses_FN)
      Read (15,'(A)') maxes_FN
        maxes_FN = Trim(maxes_FN)
      Read (15,'(A)') start_FN
        start_FN = Trim(start_FN)
      Read (15,'(A)') UserIns_FN
        UserIns_FN = Trim(UserIns_FN)

    Close(15)
!   -----------------------------------------
  
  
    Print *, 'Use requested files: '
    Print 1002, "PE:                ", PE_FN
    Print 1002, "Masses:            ", masses_FN
    Print 1002, "Coordinate limits: ", maxes_FN
    Print 1002, "User input file:   ", UserIns_FN
    Print *
    Print *, 'Reading data from suggested files'
    
   
!   ---- Read in user input system paramters from file ----  
    Open(UNIT = 15, FILE = UserIns_FN, STATUS='OLD',  ACTION ='Read')
    
      Read (15,*) EMax
      Read (15,*) mtype
      Read (15,'(A)') Out_FN

      If(mtype==1) then
        shift = 0.5d0
      Else if (mtype==2) then
        shift = 1.0d0
      Else
         Print *, "ERROR 1: INVALID SWITCH VALUE"
      End If
    
    Close(15)
!   -----------------------------------------



!   -----Read in PE info from file-----
    Open(UNIT = 15, FILE = PE_FN, STATUS='OLD',  ACTION ='Read')
    
      Read(15,*) Dmax, coeffs_len

      
      
      Allocate(Coeffs(coeffs_len))
      Allocate(Powers(coeffs_len, Dmax))
     
      Do i = 1, coeffs_len
        Read (15,*) (Powers(i,j), j=1, Dmax), Coeffs(i)
        Coeffs(i) = Coeffs(i)
      End Do 
    
    Close(15)
!   -----------------------------

!   -----Read in PE info from file-----
    Open(UNIT = 15, FILE = start_FN, STATUS='OLD',  ACTION ='Read')
    
      Allocate(start_point(2*Dmax))

      Do i = 1, 2*Dmax
         Read(15,*) start_point(i)
      End Do
    
    Close(15)
!   -----------------------------
 

!   ----Read maxes from file-----
    Open(15, file = maxes_FN, STATUS='OLD',  ACTION ='Read')
    
      Allocate(maxes(2,Dmax))

      Read(15,*) R_max
      
      Do i = 1, Dmax
        Read(15,*) (maxes(j,i),j=1,2)
      End Do

      !Pendpts      = -20.0d0 !initialize max with small value
      !Qendpts(1,:) =  20.0d0 !initialize mins with large value
      !Qendpts(2,:) = -20.0d0 !initilaize maxes with small value     
      
    Close(15)
!   -----------------------------


!   ----Read masses from file-----
    Open(15, file = masses_FN, STATUS='OLD',  ACTION ='Read')
    
      Allocate(Masses(DMax))
     
      Do i = 1, Dmax
        Read(15,*) Masses(i)
      End Do
      
    Close(15)
!   ------------------------------


    
    Print *, "Data read: SUCCESS"
    Print *
    Print *, "    ==============================================================="
    Print 1000, "Degrees of Freedom:    ", Dmax
    Print 1000, "PS Dimensionality:     ", 2*Dmax
    Print 1000, "Number of parameters:  ", coeffs_len
    Print 1003, "Energy cuttoff:        ", EMax
    Print 1002, "Output filename:       ", Out_FN
    If (mtype == 1) then
      Print '(10x,a)',  " Lattice spacing:     half integer"
    Else if (mtype ==2) then
      Print '(10x,a)',  " Lattice spacing:     whole integer"
    Else
      Print *, "ERROR 1: INVALID SWITCH VALUE"
    End If
     Print *
    Print *, '            User Defined Search Maximums'
    Write(*,'(6x)',Advance='NO')
    Do i = 1,Dmax*7+2
      Write(*,'(a)',Advance='NO') "-"
    End Do
    Print *
    Do i = 1, 2
      Write(*,'(8x)',Advance='NO')
      Do j = 1, Dmax
        Write(*,'(f5.1, 2x)', Advance='NO') DBLE(maxes(i,j))-DBLE(Size)-DBLE(shift)
      End Do
      Print *
    End Do
    Write(*,'(6x)',Advance='NO')
    Do i = 1,Dmax*7+2
      Write(*,'(a)',Advance='NO') "-"
    End Do
    Print *
    Print *
    Print 1003, "R_max = ", R_max 
    Print *
    Print *, "    ==============================================================="
    Print *
    Print *,'*****************************************************************'
    Print *
      
      1000 Format(10x, a, i13)
      1002 Format(7x, a27, a)
      1003 Format(10x, a, f13.10)
   
    
  End Subroutine ReadUserInputs
!-----------------------------------------------------------------------------


!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
! Subroutine MakeXArrays()
!   - MakeXArrays subroutine populates the array for the specific values
!     of x^0, x^1, x^2, x^3, and x^4 at each point on the grid. This array is 
!     used to more efficiently calculate the potential energy.  
!-----------------------------------------------------------------------------
  Subroutine MakeXArrays()
  !---------------
  !    Modules 
    Use Constants 
    Use ParameterArrays
    Use UserInputs
  !---------------

    Implicit NONE

    Integer :: i, j
    Double Precision :: m, Delta

    Delta = Sqrt(PI)

!    ====  Half Integer Grid  ====
    If(mtype==1) then
      Allocate(Xpoints(5,2*Size))

      Do j = 1, (2*Size)
        Xpoints(1,j) = 1.0d0
      End Do

      Do i = 2, 5
        Do j = 1, (2*Size)
          m = 1.0d0*j - 1.0d0*SIZE - .50d0
          Xpoints(i,j) = (m*Delta)**(i-1)
        End Do
      End Do

!    ====  Whole Integer Grid ====
    Else If (mtype==2) then
      Allocate(Xpoints(5,2*Size + 1))

      Do j = 1, (2*Size)+1
        Xpoints(1,j) = 1.0d0
      End Do

      Do i = 2, 5
        Do j = 1, (2*Size)+1
          m = 1.0d0*j - 1.0d0*SIZE - 1.0d0
          Xpoints(i,j) = (m*Delta)**(i-1)
        End Do
      End Do
      
    Else
      Print *, "ERROR 3: INVALID SWITCH VALUE"
    End If

  End Subroutine MakeXArrays
!-----------------------------------------------------------------------------

  
  
!-----------------------------------------------------------------------------  
Subroutine PST()

  !---------------
  !    Modules 
    Use Constants 
    Use ParameterArrays
    Use UserInputs
    Use Outputs
  !---------------


  Implicit None
 
  Integer :: i,j
  
  Integer :: k !last value of the k-stack
  Integer :: f
  Integer :: points_len !length of the points array
  Integer :: stack_len !Length of the k-stack
  Integer :: zero_pos
  Integer :: point_index
  Integer :: Counter
  Integer :: max_counts
  Integer :: count_ticker
  
  Integer :: edge_counter
  Integer :: r_counter
  
  Double Precision :: t1_search, t2_search, t_search
  Double Precision :: t1_tot, t2_tot, t_tot
  Double Precision :: t1_H, t2_H, t_H
  Double Precision :: t_dummy
  Double Precision :: t1_test, t2_test, t_test
  Character(len = 5) :: units
  Character(len = 5) :: units_H
  Character(len = 5) :: units_search
  Character(len = 5) :: units_tot

  Logical :: k_flag
  Logical :: member_flag
  Logical :: p_flag
  logical :: q_flag

  Double Precision :: H, Radius
  Double Precision :: eM, eN
  Double Precision :: Delta
  Double Precision :: KE_tot, PE_tot
  
  Double Precision :: Maker_Specific

  Integer,          Allocatable, Dimension(:,:) :: visits
  Integer,          Allocatable, Dimension(:)   :: stack
  Integer,          Allocatable, Dimension(:)   :: Yous
  Double Precision, Allocatable, Dimension(:)   :: test_point        
  
  Print *,'*****************************************************************'
  Print *,'        STEP 2:'
  Print *,'        Searching Phase Space'
  Print *
 
  f = 2*Dmax
  points_len = 1
  stack_len = 1
  point_index = 1
  k_flag = .False.
  member_flag = .FALSE.
  Delta = SQRT(PI)
 
  t_H = 0.0d0
  t_Search = 0.0d0

  !fail safe for now, so we dont loop infinitely
  Counter = 0
  max_counts = 1000000*((2*f)+1)
  count_ticker = 20000*((2*f)+1)
  edge_counter = 0
  r_counter = 0

  Allocate(points(f,Max_len))
  Allocate(visits(2*f,Max_len))
  Allocate(stack(Max_len))
  Allocate(Yous(f))
  Allocate(test_point(f))

  !Inilizie visits to have one point that contains no visits (0,0,0,0)
  k = 1
  visits = 0
  points = -100.0d0
  stack = -100
  stack(1) = 1

  Call cpu_time(t1_search)

! ------------------------------     Check the starting point    -------------------------
  Print *, '         Initial Starting Point'
  j = 1
  Do i = 1, f
    If (Mod(i,2)==0) then
      Write(*,'(a,i2,a,3x,f10.4,3x)') "p_",j,"=",start_point(i)*Delta
      j = j + 1
    Else 
      Write(*,'(5x,a,i2,a,3x,f10.4,3x)', Advance ='NO') "q_",j,"=",start_point(i)*Delta
    End If
  End Do
  Print *
  
  KE_tot = 0.0d0
   
  Do i = 1, Dmax
    Yous(i) = INT((start_point((2*i)-1)) + SHIFT + SIZE)
    eN = start_point(2*i)
    KE_tot = KE_tot + (hbar*hbar)*(((eN)*(eN)*Delta*Delta)/(2.0d0*masses(i)))
  End Do

  PE_tot = Maker_specific(Yous, Powers, Coeffs, Xpoints, coeffs_len, Dmax, SIZE)

  H = KE_tot + PE_tot

  If(H.GT.Emax) then
    Print '(5x,f7.6,2x,a1,2x,f7.6)',H,'>',Emax
    Print *, "       !!!!!!  INVALID STARTING POINT  !!!!!!"
    Print *, "                       EXITING"
    Print *
    stack_len = -100
  Else 
    Do i = 1, f
      points(i,k) = start_point(i)
    End Do
  End If
!------------------------------------------------------------------------------------------

  Print *, '    ------------------------------------------'
  Do While (stack_len > 0 .AND. Counter < max_counts)
    If(Mod(Counter, count_ticker) == 0 .AND. counter.NE.0) Then
      Print *, '       Current basis size:   ', points_len
      Print *, '       Current stack length: ', stack_len
      Print *
    End If
  
  !Get the test point from the stack
    k = stack(stack_len)
    Do i = 1, f
      test_point(i) = points(i,k)
    End Do

    !Do i = 1, f
    !  Write(*,'(f4.1,2x)',Advance='NO') test_point(i)
    !End Do
    !Print '(f8.7)', PE_tot
    !Print *, k


   !Look at visits and figure out where to go
    zero_pos = -1
    i = 1
    Do While (zero_pos == -1 .AND. i < 2*f + 1) 
      If(visits(i,k) == 0) Then
        visits(i,k) = 1
        zero_pos = i
      End If
      i = i + 1 
    End Do

    If (zero_pos == -1) Then !point is finished
     
      Call cpu_time(t1_test)
      Call StackIt(stack, stack_len, Max_len, point_index, -1)
      Call cpu_time(t2_test)
      t_test = t_test + (t2_test-t1_test)
      point_index = k
      
    Else  !move to next test point
      
      If (zero_pos .LE. f) Then
        !+1
        test_point( Mod(zero_pos,f)+1 ) =  test_point( Mod(zero_pos,f)+1 ) + 1.0d0
      Else
        !-1
        test_point( Mod(zero_pos,f)+1 ) =  test_point( Mod(zero_pos,f)+1 ) - 1.0d0
      End If 


    !---------------------------   Check if it H<Emax    ---------------------------------

      
      !Convert grid points to VN lattice Points, X -> m's, n's
      Call cpu_time(t1_H)
      KE_tot = 0.0d0
      Radius = 0.0d0
      p_Flag = .True.
      q_Flag = .True.
      Do i = 1, Dmax
        Yous(i) = INT((test_point((2*i)-1)) + SHIFT + SIZE)
        eM = test_point((2*i)-1)
        eN = test_point(2*i)
        Radius = Radius + (eM*Delta*eM*Delta)
        If(eN<0.0d0)then
          p_Flag = .False.
        End If
        If(yous(i) < maxes(1,i) .OR. yous(i) > maxes(2,i)) then
          edge_counter = edge_counter + 1
          q_flag = .False.
        End If
        KE_tot = KE_tot + (hbar*hbar)*(((eN)*(eN)*Delta*Delta)/(2.0d0*masses(i)))
      End Do

      Radius = Sqrt(Radius)
      If(Radius > R_max) then
        r_counter = r_counter + 1
        q_flag = .False.
      End If

      PE_tot = Maker_specific(Yous, Powers, Coeffs, Xpoints, coeffs_len, Dmax, SIZE)
    
      H = KE_tot + PE_tot
      Call cpu_time(t2_H)
      t_H = t_H + (t2_H-t1_H)
! -------------------------------------------------------------------------------------------

      If(H.LE.Emax .AND. p_flag .AND. q_flag) Then
     
        Call cpu_time(t1_search)
        Call SearchList(member_flag, test_point, f, points, points_len, max_len)
        Call cpu_time(t2_search)
        t_search = t_search + (t2_search-t1_search)
        
        !Add it to the list
        If(.NOT.member_flag) then
          
          points_len = points_len + 1
          point_index = points_len

          Call StackIt(stack, stack_len, Max_len, point_index, 1)

          Do i = 1, f
            points(i,point_index) = test_point(i)
          End Do
          
        End if

      End If

    End If

  Counter = counter + 1
  
  End Do  

  Call cpu_time(t2_tot)
  t_tot = t2_tot - t1_tot
  iMax = points_len

  Deallocate(visits)
  Deallocate(stack)

  t_dummy = t_search
  If(t_dummy > 18000.0d0) then
    t_dummy = (t_dummy)/(60.0d0*60.0d0)
    units = 'hrs.'
  Else if (t_dummy > 300.0d0) then
    t_dummy = (t_dummy)/(60.0d0)
    units = 'mins.'
  Else 
    units= 'sec.'
  End If
  t_search = t_dummy
  units_search = units

  t_dummy = t_H
  If(t_dummy > 18000.0d0) then
    t_dummy = (t_dummy)/(60.0d0*60.0d0)
    units = 'hrs.'
  Else if (t_dummy > 300.0d0) then
    t_dummy = (t_dummy)/(60.0d0)
    units = 'mins.'
  Else 
    units= 'sec.'
  End If
  t_H = t_dummy
  units_H = units

  t_dummy = t_tot
  If(t_dummy > 18000.0d0) then
    t_dummy = (t_dummy)/(60.0d0*60.0d0)
    units = 'hrs.'
  Else if (t_dummy > 300.0d0) then
    t_dummy = (t_dummy)/(60.0d0)
    units = 'mins.'
  Else 
    units= 'sec.'
  End If
  t_tot = t_dummy
  units_tot = units

  Print *, '    ------------------------------------------'
  Print *
  Print *
  Print *, '    ========================================='
  Print *, '            Search Complete'
  Print *
  Print *, '      q_max reached: ', edge_counter
  Print *, '      r_max reached: ', r_counter
  Print *, '      total iterations: ', counter
  Print *  
  Print *, "      Total Points: ", points_len
  Print *, '    ========================================='
  Print *
  Print *, '    ========================================='
  Print *, '              Timing Summery'
  Print '(a,f10.5,2x,a)', '    Total Time:  ', t_tot, units_tot
  Print '(a,f10.5,2x,a)', '    Search time: ', t_search, units_search
  Print '(a,f10.5,2x,a)', '    H calc time: ', t_H, units_H
  !Print '(a,f10.5)', 'test time: ', t_test
  Print *, '    ========================================='
  Print *,'*****************************************************************'
   
End Subroutine PST
!-----------------------------------------------------------------------------



!-----------------------------------------------------------------------------
Subroutine StackIt(stack, stack_len, Max_len, point_index, x)

  Implicit NONE
  Integer, Intent(INOUT) :: stack_len
  Integer, Intent(IN)    :: Max_len
  Integer, Intent(INOUT) :: point_index
  Integer, Intent(IN)    :: x

  integer :: i

  Integer, Intent(INOUT), Dimension(Max_Len) :: stack
  Integer,                Dimension(Max_Len) :: dummy_stack

  dummy_stack = -100
  
  If (x == 1) Then !increase

   !====================================
   !--    Increase stack routine   -----
    Do i = 1, stack_len
      dummy_stack(i) = stack(i)
    End Do
    stack_len = stack_len + 1
    dummy_stack(stack_len) = point_index
   
    stack = -100
        
    Do i = 1, stack_len
      stack(i) = dummy_stack(i)
    End Do
   !====================================

  Else If (x == -1) Then
    !====================================
    !--    Decrease stack routine   -----
    stack_len = stack_len - 1
    Do i = 1, stack_len
      dummy_stack(i) = stack(i)
    End Do      
    
     stack = -100
    
    Do i = 1, stack_len
      stack(i) = dummy_stack(i)
    End Do
   
   !====================================
 End If
   
End Subroutine StackIt
!-----------------------------------------------------------------------------




!-----------------------------------------------------------------------------
Subroutine SearchList(member_flag, test_point, f, points, points_len, Max_len)

  Implicit None
  Logical, Intent(INOUT) :: member_flag
  Integer, Intent(IN)    :: f
  Integer, Intent(IN)    :: points_len
  Integer, Intent(IN)    :: Max_len

  Logical :: search_flag
  Integer :: search_counter
  Integer :: coord_counter

  Double Precision, Intent(IN), Dimension(f) :: test_point
  Double Precision, Intent(IN), Dimension(f,max_len) :: points
     
      

  member_flag = .FALSE.
  search_counter = 1
  Do While(.NOT.member_Flag .And. search_counter < points_len + 1)
  
    search_flag = .FALSE.
    coord_counter = 1
    Do While (.NOT.search_flag .AND. coord_counter < f + 1)
      If(test_point(coord_counter) == points(coord_counter, search_counter)) then
        search_flag = .FALSE.
      Else
        search_flag = .TRUE.
      End If
      coord_counter = coord_counter + 1
    End Do
  
    If(.NOT.search_flag) then
      member_flag = .TRUE.
    End If
  
    search_counter = search_counter + 1
  End Do

End Subroutine SearchList
!-----------------------------------------------------------------------------




!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
! Subroutine PrintVals(Emax)
!   - PrintVals subroutine allows the user to output a subset of the
!     eigenvalues to the screen. The user inputs a range, <min, max>.
!
!     Emax
!     - Double precision. Input parameter given by the user.  
!
!-----------------------------------------------------------------------------
  Subroutine WriteFile()
    
   !----Modules-----
    Use Outputs
    Use UserInputs
   !----------------  

	Implicit NONE
    Integer	:: i,j
   

    Print *
    Print *,'*****************************************************************'
    Print *,'       SUMMARY'
    Print *
    Print 1000, (Dmax+6)/3
	Print 1001, EMax
	Print 1002, iMax

 	Open(UNIT = 15, FILE = Out_FN, ACTION='WRITE')

    Write(15, "(2x,i1,3x,f13.5,'d0',3x,i7)") mtype, EMax, iMax
    
	Do i = 1, iMax
      Do j = 1, 2*Dmax
        Write (15,"(f5.1,'d0', 2x)", Advance='NO') points(j,i)
      End Do
      Write (15,*)
    End Do
  
    Close(UNIT = 15)
    
    Print 1003, Out_FN
    Print *

    1000 Format (2x,i2, ' Atom System') 
    1001 Format ('  Energy cut off:                     ', 3x, f13.5)
    1002 Format ('  Total number of phase space points: ', 3x, i7)
    1003 Format ('  Points read to file:                ', a128)

  End Subroutine WriteFile
!-----------------------------------------------------------------------------





     !**********************************************************!
     !********************    Functions   **********************!
     !**********************************************************!


 
!-------------------------------------------------------------------------------
Double Precision Function Maker_general(exes, A, B, m, n)
  !-----------------------------------------------------------------------------!
  ! Arguments:                                                                  !
  !                                                                             !
  ! exes-                                                                       !
  ! Double precision array of length n. The point at which the function         ! 
  ! is to be evaluated                                                          !
  !                                                                             !
  ! A                                                                           !
  ! Integer Array of dimeniosn m x n. The ordered list of powers to be used     !
  ! in the polynomial.                                                          !
  !                                                                             !
  ! B                                                                           !
  ! Double precision array of length m. The coefficient assoicated with the     !  
  ! mth term in the polynomial expression                                       !
  !                                                                             !
  ! m                                                                           !
  ! Integer value. Represents the total number of terms in the polynomial       !
  !                                                                             !
  ! n                                                                           !
  ! Integer value. Represents the total number of independent variables         !
  ! in the polynomial fucntion                                                  !
  !-----------------------------------------------------------------------------! 


  Implicit NONE
  Integer, Intent(IN) :: m, n   
  Double Precision, Intent(IN), Dimension(n)   :: exes
  Double Precision, Intent(IN), Dimension(m)   :: B 
  Integer,          Intent(IN), Dimension(m,n) :: A 
 

  Double precision :: dummy, pot, point
  Integer :: i,j
  integer :: power

  pot = 0.0d0

  Do i = 1, m
    dummy = 1.0d0
    Do j = 1 , n
      power = A(i,j)
      point = exes(j)
      If (point .LT. 0.0d0) then
        dummy = dummy * (point**power)
      Else If (point .GT. 0.0d0) then
        dummy = dummy * (point**power)
      Else
        dummy = 0.0d0
      End If
    End Do
    pot = pot + B(i)*dummy
  End Do

  Maker_general = pot
  
End Function  Maker_general
!-------------------------------------------------------------------------------


!-------------------------------------------------------------------------------
Double Precision Function Maker_specific(exes, A, B, C, m, n, l)
!-----------------------------------------------------------------------------!
  ! Arguments:                                                                  !
  !                                                                             !
  ! exes-                                                                       !
  ! Double precision array of length n. The point at which the function         ! 
  ! is to be evaluated                                                          !
  !                                                                             !
  ! A                                                                           !
  ! Integer Array of dimeniosn m x n. The ordered list of powers to be used     !
  ! in the polynomial.                                                          !
  !                                                                             !
  ! B                                                                           !
  ! Double precision array of length m. The coefficient assoicated with the     !  
  ! mth term in the polynomial expression                                       !
  !                                                                             !
  ! C                                                                           !
  ! Double Precision array of dimension 5 x 2*Size. Contains the specific       ! 
  ! values of the individual powers of x at the grid points                     !
  !                                                                             !
  ! m                                                                           !
  ! Integer value. Represents the total number of terms in the polynomial       !
  !                                                                             !
  ! n                                                                           !
  ! Integer value. Represents the total number of independent variables         !
  ! in the polynomial fucntion                                                  !
  !                                                                             !
  ! l                                                                           !
  !  Integer value. The length of the gid from the mid point to the maximum     !
  !  x grid value.                                                              !
  !-----------------------------------------------------------------------------! 

  Implicit NONE
  Integer, Intent(IN) :: m, n, l   
  Integer,          Intent(IN), Dimension(n)     :: exes
  Double Precision, Intent(IN), Dimension(m)     :: B 
  Integer,          Intent(IN), Dimension(m,n)   :: A
  Double Precision, Intent(IN), Dimension(5,2*l) :: C  

  Double precision :: dummy, pot
  Integer :: i,j
  integer :: power, point

 
  pot = 0.0d0

  Do i = 1, m
    dummy = 1.0d0
    Do j = 1 , n
      power = A(i,j) + 1
      point = exes(j)
      dummy = dummy * C( power, point )
    End Do
    pot = pot + B(i)*dummy
  End Do

  Maker_specific = pot

End Function  Maker_specific
!==============================================================
!==============================================================