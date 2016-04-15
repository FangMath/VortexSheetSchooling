#include <fintrf.h>
C
C LEGEWHTSG.FOR - Gateway function for LEGEWHTS.FOR
C
C This is an example of the FORTRAN code required for interfacing
C a .MEX file to MATLAB.
C
C This subroutine is the main gateway to MATLAB.  When a MEX function
C  is executed MATLAB calls the MEXFUNCTION subroutine in the corresponding
C  MEX file.  
C
C Copyright 1984-2006 The MathWorks, Inc.
C 
C
      SUBROUTINE MEXFUNCTION(NLHS, PLHS, NRHS, PRHS)
C-----------------------------------------------------------------------
C     (pointer) Replace integer by integer*8 on 64-bit platforms
C
      MWPOINTER PLHS(*), PRHS(*)

C-----------------------------------------------------------------------
C

      INTEGER NLHS, NRHS
C
C-----------------------------------------------------------------------
C     (pointer) Replace integer by integer*8 on 64-bit platforms
C
      MWPOINTER MXCREATEDOUBLEMATRIX, MXGETPR

C-----------------------------------------------------------------------
C

      MWSIZE MXGETM, MXGETN
C
C KEEP THE ABOVE SUBROUTINE, ARGUMENT, AND FUNCTION DECLARATIONS FOR USE
C IN ALL YOUR FORTRAN MEX FILES.
C---------------------------------------------------------------------
C
C-----------------------------------------------------------------------
C     (pointer) Replace integer by integer*8 on 64-bit platforms
C
      MWPOINTER ns_ptr, vals_ptr_r, vals_ptr_i
      MWPOINTER nt_ptr, z_ptr_r, z_ptr_i
      MWPOINTER pot_ptr_r, pot_ptr_i

C-----------------------------------------------------------------------
C

      MWSIZE NEL, local_int_ns, local_int_nt
      REAL*8 local_ns, local_vals_r(10000), local_vals_i(10000)
      REAL*8 local_nt, local_z_r(10000), local_z_i(10000) 
      complex*16 local_vals(10000), local_pot(10000), ima, z, pot
      integer*4 iter
      ima = (0,1)

#if defined MSWIND
C For Windows only!
C This resets the floating point exception to allow divide by zero,
C overflow and invalid numbers. 
C
	INTEGER(2) CONTROL
	CALL GETCONTROLFPQQ(CONTROL)
	CONTROL = CONTROL .OR. FPCW$ZERODIVIDE
      CONTROL = CONTROL .OR. FPCW$INVALID
      CONTROL = CONTROL .OR. FPCW$OVERFLOW
	CALL SETCONTROLFPQQ(CONTROL)
#endif

C
C CHECK FOR PROPER NUMBER OF ARGUMENTS
C
      IF (NRHS .NE. 6) THEN
        CALL MEXERRMSGTXT('LEGEWHTS requires six input arguments')
      ELSEIF (NLHS .GT. 3) THEN
        CALL MEXERRMSGTXT('LEGEWHTS requires two output argument')
      ENDIF
C
C CHECK THE DIMENSIONS OF Y.  IT CAN BE 4 X 1 OR 1 X 4.
C
C
      ns_ptr  = MXGETPR(PRHS(1))
      vals_ptr_r = MXGETPR(PRHS(2)) 
      vals_ptr_i = MXGETPR(PRHS(3)) 
      nt_ptr  = MXGETPR(PRHS(4))
      z_ptr_r = MXGETPR(PRHS(5)) 
      z_ptr_i = MXGETPR(PRHS(6)) 
      
	    if(mxGetNumberOfElements(prhs(1)) .ne. 1) then
	      call mexErrMsgIdAndTxt(
     +       'MATLAB:hilbert:invalidSize',
     +       '1st argument ns must be integer scalar.')
	    end if
	    if(mxGetNumberOfElements(prhs(4)) .ne. 1) then
	      call mexErrMsgIdAndTxt(
     +       'MATLAB:hilbert:invalidSize',
     +       '4th argument nt must be integer scalar.')
	    end if
C
C
C COPY RIGHT HAND ARGUMENTS TO LOCAL ARRAYS OR VARIABLES
      NEL = 1
      CALL MXCOPYPTRTOREAL8(ns_ptr, local_ns, NEL)
      CALL MXCOPYPTRTOREAL8(nt_ptr, local_nt, NEL)
C
      local_int_ns = int(local_ns)
      if (local_int_ns .gt. 900) then
	      call mexErrMsgIdAndTxt(
     +       'Bomb!',
     +       '1st argument ns is too large.')
      endif
      local_int_nt = int(local_nt)
      if (local_int_nt .gt. 10000) then
	      call mexErrMsgIdAndTxt(
     +       'Bomb!',
     +       '4th argument nt is too large.')
      endif
C
C
	    if((mxGetNumberOfElements(prhs(2)) .ne. local_int_ns) .or. 
     +    (mxGetNumberOfElements(prhs(3)) .ne. local_int_ns)) then
	      call mexErrMsgIdAndTxt(
     +       'MATLAB:hilbert:invalidSize',
     +       '2nd-3rd arguments Vals must be length ns.')
	    end if
	    if((mxGetNumberOfElements(prhs(5)) .ne. local_int_nt) .or. 
     +    (mxGetNumberOfElements(prhs(6)) .ne. local_int_nt)) then
	      call mexErrMsgIdAndTxt(
     +       'MATLAB:hilbert:invalidSize',
     +       '5th-6th arguments Vals must be length nt.')
	    end if
      NEL = local_int_ns
      CALL MXCOPYPTRTOREAL8(vals_ptr_r, local_vals_r, NEL)
      NEL = local_int_ns
      CALL MXCOPYPTRTOREAL8(vals_ptr_i, local_vals_i, NEL)
      local_vals = local_vals_r + ima*local_vals_i
C
C
      NEL = local_int_nt 
      CALL MXCOPYPTRTOREAL8(z_ptr_r, local_z_r, NEL)
      CALL MXCOPYPTRTOREAL8(z_ptr_i, local_z_i, NEL)
C
C
C CREATE A MATRIX FOR RETURN ARGUMENT
C
      PLHS(1) = MXCREATEDOUBLEMATRIX(1,local_int_nt,0)
      PLHS(2) = MXCREATEDOUBLEMATRIX(1,local_int_nt,0)
C
C ASSIGN POINTERS TO THE VARIOUS PARAMETERS
C
      pot_ptr_r = MXGETPR(PLHS(1))
      pot_ptr_i = MXGETPR(PLHS(2))
C
C DO THE ACTUAL COMPUTATIONS IN A SUBROUTINE
C       CREATED ARRAYS.  
C
        do iter = 1,local_int_nt
        z = local_z_r(iter) + ima*local_z_i(iter)
        call hilbert_legendre(z, 
     +       local_int_ns, local_vals, pot)
        local_pot(iter) = pot
        enddo
C
C COPY OUTPUT WHICH IS STORED IN LOCAL ARRAY TO MATRIX OUTPUT

      NEL = local_int_nt 
      CALL MXCOPYREAL8TOPTR(real(local_pot), pot_ptr_r, NEL)
      NEL = local_int_nt 
      CALL MXCOPYREAL8TOPTR(imag(local_pot), pot_ptr_i, NEL)
C
      RETURN
      END

C        call hilbert_legendre(z, k, vals, pot)
C        call legewhts(k, xs, whts, ifwhts)
C      prhs(1)->ns; prhs(2,3)->vals; prhs(4)->nt; prhs(5,6)->z; plhs(1,2)->pot
C      k -> integer 
C      vals-> input vector of length k
C      z-> input vector of length 1
C      pot-> output vector of length 1 
