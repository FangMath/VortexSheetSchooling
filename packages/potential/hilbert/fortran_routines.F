cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c         frome legeexps.f
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine legewhts(n,ts,whts,ifwhts)
        implicit real *8 (a-h,o-z)
        dimension ts(1),whts(1),ws2(1000),rats(1000)
c
c        this subroutine constructs the nodes and the
c        weights of the n-point gaussian quadrature on 
c        the interval [-1,1]
c
c                input parameters:
c
c  n - the number of nodes in the quadrature
c
c                output parameters:
c
c  ts - the nodes of the n-point gaussian quadrature
c  w - the weights of the n-point gaussian quadrature
c
c       . . . construct the array of initial approximations
c             to the roots of the n-th legendre polynomial
c
        eps=1.0d-14
        ZERO=0
        DONE=1
        pi=datan(done)*4
        h=pi/(2*n) 
        do 1200 i=1,n
        t=(2*i-1)*h
        ts(n-i+1)=dcos(t)
1200  CONTINUE
c
c         use newton to find all roots of the legendre polynomial
c
        ts(n/2+1)=0
        do 2000 i=1,n/2
c
        xk=ts(i)
        ifout=0
        deltold=1
        do 1400 k=1,10
        call legepol_sum(xk,n,pol,der,sum)
        delta=-pol/der
        xk=xk+delta
        if(abs(delta) .lt. eps) ifout=ifout+1
c
        if(ifout .eq. 3) goto 1600
 1400 continue
 1600 continue
        ts(i)=xk
        ts(n-i+1)=-xk
 2000 continue
c     
c        construct the weights via the orthogonality relation
c
        if(ifwhts .eq. 0) return
c
        do 2400 i=1,(n+1)/2
        call legepol_sum(ts(i),n,pol,der,sum)
        whts(i)=1/sum
        whts(n-i+1)=whts(i)
 2400 continue
c
        return
        end

c

        subroutine legepol_sum(x,n,pol,der,sum)
        implicit real *8 (a-h,o-z)
c
        done=1
        sum=0 
c
        pkm1=1
        pk=x
        sum=sum+pkm1**2 /2
        sum=sum+pk**2 *(1+done/2)
c
        pk=1
        pkp1=x
c
c        if n=0 or n=1 - exit
c
        if(n .ge. 2) goto 1200

        sum=0 
c
        pol=1
        der=0
        sum=sum+pol**2 /2
        if(n .eq. 0) return
c
        pol=x
        der=1
        sum=sum+pol**2*(1+done/2)
        return
 1200 continue
c
c       n is greater than 1. conduct recursion
c
        do 2000 k=1,n-1
        pkm1=pk
        pk=pkp1
        pkp1=( (2*k+1)*x*pk-k*pkm1 )/(k+1)
        sum=sum+pkp1**2*(k+1+done/2)
 2000 continue
c
c        calculate the derivative
c
        pol=pkp1
        der=n*(x*pkp1-pk)/(x**2-1)
        return
        end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        


cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c        used for hilbert_legendre()
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

        subroutine legeexps(itype,n,x,u,v,whts)
        implicit real *8 (a-h,o-z)
        dimension x(1),whts(1),u(n,n),v(n,n)
c
c         this subroutine constructs the gaussiaqn nodes 
c         on the interval [-1,1], and the weights for the 
c         corresponding order n quadrature. it also constructs
c         the matrix v converting the coefficients
c         of a legendre expansion into its values at the n
c         gaussian nodes, and its inverse u, converting the
c         values of a function at n gaussian nodes into the
c         coefficients of the corresponding legendre series.
c         no attempt has been made to make this code efficient, 
c         but its speed is normally sufficient, and it is 
c         mercifully short.
c
c                 input parameters:
c
c  itype - the type of the calculation to be performed
c          itype=0 means that only the gaussian nodes are 
c                  to be constructed. 
c          itype=1 means that only the nodes and the weights 
c                  are to be constructed
c          itype=2 means that the nodes, the weights, and
c                  the matrices u, v are to be constructed
c  n - the number of gaussian nodes and weights to be generated
c  
c                 output parameters:
c
c  x - the order n gaussian nodes - computed independently
c          of the value of itype.
c  u - the n*n matrix converting the  values at of a polynomial of order
c         n-1 at n legendre nodes into the coefficients of its 
c         legendre expansion - computed only in itype=2
c  v - the n*n matrix converting the coefficients
c         of an n-term legendre expansion into its values at
c         n legendre nodes (note that v is the inverse of u)
c          - computed only in itype=2
c  whts - the corresponding quadrature weights - computed only 
c         if itype .ge. 1
c
c       . . . construct the nodes and the weights of the n-point gaussian 
c             quadrature
c
        ifwhts=0
        if(itype. gt. 0) ifwhts=1
        call legewhts(n,x,whts,ifwhts)
c
c       construct the matrix of values of the legendre polynomials
c       at these nodes        
c
        if(itype .ne. 2) return
        do 1400 i=1,n
c
        call legepols(x(i),n-1,u(1,i) )
 1400 continue
c
        do 1800 i=1,n
        do 1600 j=1,n
        v(i,j)=u(j,i)
 1600 continue
 1800 continue
c
c       now, v converts coefficients of a legendre expansion
c       into its values at the gaussian nodes. construct its 
c       inverse u, converting the values of a function at 
c       gaussian nodes into the coefficients of a legendre 
c       expansion of that function
c
        do 2800 i=1,n
        d=1
        d=d*(2*i-1)/2
        do 2600 j=1,n
        u(i,j)=v(j,i)*whts(j)*d
 2600 continue
 2800 continue
        return
        end


        subroutine legepols(x,n,pols)
        implicit real *8 (a-h,o-z)
        dimension pols(*)
c
        pkm1=1
        pk=x
c
        pk=1
        pkp1=x
c
c
c        if n=0 or n=1 - exit
c
        if(n .ge. 2) goto 1200
        pols(1)=1
        if(n .eq. 0) return
c
        pols(2)=x
        return
 1200 continue
c
        pols(1)=1
        pols(2)=x
c
c       n is greater than 2. conduct recursion
c
        do 2000 k=1,n-1
        pkm1=pk
        pk=pkp1
        pkp1=( (2*k+1)*x*pk-k*pkm1 )/(k+1)
        pols(k+2)=pkp1
 2000 continue
c
        return
        end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc









cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c         frome qlegefuns.f
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine hilbert_legendre(z, k, vals, pot)
        implicit real *8 (a-h,o-z)
        complex *16 ima, z, vals(k), pot
        real *8 time_start, time_end, time_diff
c
        data ima/(0,1)/
        real *8 xs(1000), ys(1000), xnodes(1000)
        real *8 u(1000000), v(1000000), whts(1000)
        real *8 coefs_real(0:1000), coefs_imag(0:1000)
        complex *16 coefs(0:1000), qfuns(0:1000), zterms(0:1000)
c
c       this routine computes the integral:
c
c          pot = \int_{-1}^1  f(x)/(z - x) dx
c
c       where f(x) is a (complex) function sampled at k legendre 
c       nodes on [-1,1]. The target z can be arbitrary, as the
c       Hilbert transform is evaluated using the following P_n, Q_n
c       identity
c
c          Q_n = \frac{1}{2} \int_{-1}^1 P_n / (z - x) dx
c
c       all of the following calculations are done assuming that we
c       define the branch cut in Q_n along the interval [-1,1] so that
c       there is a jump in the Hilbert transform across [-1,1].
c
c       input:
c         z - target point
c         k - number of legendre nodes at which f is evaluated
c         vals - the values of f at the k legendre nodes, complex valued
c
c       output:
c         pot - the Hilbert transform of f
c
c
c       get the expansion coefficients of vals
c
ccc        if (k .gt. 900) then
ccc          call prinf('k = *', k, 1)
ccc          call prinf('bomb!! k is too large*', ima, 0)
ccc          stop
ccc        endif
c
        itype = 2
        call legeexps(itype, k, xnodes, u, v, whts)
c        
        do i = 1,k
          xs(i) = vals(i)
          ys(i) = -ima*vals(i)
        enddo
c
        call matvec(k, u, xs, coefs_real)
        call matvec(k, u, ys, coefs_imag)
c
        do i = 0,k-1
          coefs(i) = coefs_real(i) + ima*coefs_imag(i)
        enddo

c
c       evaluate the first k-1 Q_n's
c
c        call cpu_time(time_start)
c        do f = 1, 100000
        call zqneval(z, k-1, qfuns)
c        enddo
c
cccc        call prin2('coefs = *', coefs, 2*nterms)
cccc        call prin2('qfuns = *', qfuns, 2*nterms)

c
c       and evaluate the expansion
c
        pot = 0
        do i = 0,k-1
          pot = pot + 2*coefs(i)*qfuns(i)
cccc          zterms(i) = 2*coefs(i)*qfuns(i)
        enddo
c        call cpu_time(time_end)
c        call prin2('subroutine time_diff = *', time_end - time_start, 1)
c
        return
        end


        subroutine matvec(n, a, x, y)
        implicit real *8 (a-h,o-z)
        real *8 a(n,n), x(n), y(n)
c
        do i = 1,n
          d = 0
          do j = 1,n
            d = d + a(i,j)*x(j)
          enddo
          y(i) = d
        enddo
c
        return
        end


        subroutine zqneval(z, n, qfuns)
        implicit real *8 (a-h,o-z)
        complex *16 z, qfuns(0:n), d, q0, q1, qnext, qn, qfuns2(0:10010)
        complex *16 ratio
c
c       Evaluate Q_0 .. Q_n at the point z in the complex plane.
c       Taking the branch cut along [-1,1] means that Q_n decays
c       as z \to \infty. We must use an up, down, and normalize
c       recursion.
c
c       input:
c         z - the complex target point
c         n - max order of Q_n to evaluate
c
c       ouput:
c         qfuns - returns Q_0(z), ..., Q_n(z) with the INDEX 
c             STARTING AT ZERO
c
c       first, construct Q_0(z) and Q_1(z)
c
c
        call zqlege01(z, q0, q1)
        qfuns(0) = q0
        if (n .eq. 0) return
c
        qfuns(1) = q1
        if (n .eq. 1) return
c
c       recurse up until things have blown up sufficiently
c 
        qlarge = 1.0d16/abs(q0)
        qlarge = qlarge**2
c
        nmax = 10000
        ntop = 0
c
        do i=1,nmax
          qn = ( (2*i+1)*z*q1-i*q0 ) /(i+1) 
          q0 = q1
          q1 = qn
          if (abs(qn) .gt. qlarge) then
            ntop = i+1
            goto 1300
          endif 
        enddo
c 
 1300 continue

c
c       if it did not blowup, just recurse up
c
        if (ntop .eq. 0) then
          call zqneval_up(z, n, qfuns)
          return
        endif
c
c       otherwise recurse down
c
        ntop = ntop+5
c
        qfuns2(ntop) = 1
        qfuns2(ntop+1) = 0
c
        do i = ntop,1,-1
          qfuns2(i-1) = (2*i+1)*z*qfuns2(i)/i - (i+1)*qfuns2(i+1)/i
        enddo
c
c       and normalize
c
        call zqlege01(z, q0, q1)
        ratio = q0/qfuns2(0)
c
        do i = 0,ntop
          qfuns2(i) = qfuns2(i)*ratio
        enddo

        do i = 0,n
            qfuns(i) = qfuns2(i)
        enddo

        return
        end




        subroutine zqlege01(z, q0, q1)
        implicit real *8 (a-h,o-z)
        complex *16 z, q0, q1, d
c
c       this routine evaluates Q_0 and Q_1 for complex argument
c       *off* of the real line - this assumes the branch cut is
c       on [-1,1]
c
c       input:
c         z - a complex number
c
c       output:
c         q0 - the value of Q_0(z)  
c         q1 - the value of Q_1(z)  
c
        done = 1
        d = log((done+z)/(z - done))
        q0 = d/2
        q1 = z/2*d-1
c
        return
        end









        subroutine zqneval_up(z, n, qfuns)
        implicit real *8 (a-h,o-z)
        complex *16 z, qfuns(0:n), d, q0, q1
c
c       Evaluate Q_0 .. Q_n at the point z in the complex plane
c       ONLY USING AN UPWARD RECURRENCE!!! (see zqneval).
c       Branch cut is assumed to be on [-1,1].
c
c       input:
c         z - the complex target point
c         n - max order of Q_n to evaluate
c
c       ouput:
c         qfuns - returns Q_0(z), ..., Q_n(z) with the INDEX 
c             STARTING AT ZERO
c
c       NOTE: 
c       If you want to evaluate Q_n on the real-line, use the 
c       routine qneval - this routine will NOT return the principal 
c       value of Q_n for z with arbitrarily small imaginary parts.
c
c       NOTE 2: The branch cut is taken to be on [-1,1], consistent with
c       the identity:
c
c               Q_n(z) = 1/2 \int_{-1}^1 P_n(x)/(z - x) dx
c
c       first, construct Q_0(z) and Q_1(z)
c 
        call zqlege01(z, q0, q1)
        qfuns(0) = q0
        if (n .eq. 0) return
c
        qfuns(1) = q1
        if (n .eq. 1) return
c
c       recurse up until 
c 
        do i=1,n-1
          qfuns(i+1)=( (2*i+1)*z*qfuns(i)-i*qfuns(i-1) ) /(i+1)
        enddo
c 
        return
        end
