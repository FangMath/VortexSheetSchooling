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
      MWPOINTER k_ptr, ifwhts_ptr, xs_ptr, whts_ptr

C-----------------------------------------------------------------------
C

      MWSIZE NEL, local_int_k, local_int_ifwhts 
      REAL*8 local_xs(10000), local_whts(10000)
      REAL*8 local_k, local_ifwhts

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
      IF (NRHS .NE. 2) THEN
        CALL MEXERRMSGTXT('LEGEWHTS requires two input arguments')
      ELSEIF (NLHS .GT. 3) THEN
        CALL MEXERRMSGTXT('LEGEWHTS requires one output argument')
      ENDIF
C
C CHECK THE DIMENSIONS OF Y.  IT CAN BE 4 X 1 OR 1 X 4.
C
C
      k_ptr  = MXGETPR(PRHS(1))
      ifwhts_ptr  = MXGETPR(PRHS(2))
      
	    if(mxGetNumberOfElements(prhs(1)) .ne. 1) then
	      call mexErrMsgIdAndTxt(
     +       'MATLAB:mxsetdimensionsf:invalidSize',
     +       'Size arguments must be integer scalar.')
	    end if
C
	    if(mxGetNumberOfElements(prhs(2)) .ne. 1) then
	      call mexErrMsgIdAndTxt(
     +       'MATLAB:mxsetdimensionsf:invalidSize',
     +       'Ifwhts argument must be integer scalar.')
	    end if
C
C COPY RIGHT HAND ARGUMENTS TO LOCAL ARRAYS OR VARIABLES
      NEL = 1
      CALL MXCOPYPTRTOREAL8(k_ptr, local_k, NEL)
      local_int_k = int(local_k)
      
      NEL = 1
      CALL MXCOPYPTRTOREAL8(ifwhts_ptr, local_ifwhts, NEL)
      local_int_ifwhts = int(local_ifwhts)
C
C
C CREATE A MATRIX FOR RETURN ARGUMENT
C
      PLHS(1) = MXCREATEDOUBLEMATRIX(1,local_int_k,0)
      PLHS(2) = MXCREATEDOUBLEMATRIX(1,local_int_k,0)
C
C ASSIGN POINTERS TO THE VARIOUS PARAMETERS
C
      xs_ptr = MXGETPR(PLHS(1))
      whts_ptr = MXGETPR(PLHS(2))
C
C DO THE ACTUAL COMPUTATIONS IN A SUBROUTINE
C       CREATED ARRAYS.  
C
      CALL legewhts(local_int_k,
     +        local_xs, local_whts, local_int_ifwhts)
C
C COPY OUTPUT WHICH IS STORED IN LOCAL ARRAY TO MATRIX OUTPUT

      NEL = local_int_k 
      CALL MXCOPYREAL8TOPTR(local_xs, xs_ptr, NEL)
      NEL = local_int_k 
      CALL MXCOPYREAL8TOPTR(local_whts, whts_ptr, NEL)
C
      RETURN
      END

C        call legewhts(k, xs, whts, ifwhts)
C      prhs(1)->k; prhs(2)->ifwhts; plhs(1)->xs; plhs(2)->whts
C      k -> integer 
C      ifwhts -> flag
C      xs -> output vector of length k
C      whts -> output vector of length k 
