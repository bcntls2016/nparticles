program nparticles
	implicit none
	character (len = 40)	:: fileden = "density.jogger.dat"	! input file containing the wavefunction
	character (len = 1)		:: cchar = "#"			! character to be skipped on routine titols
	logical (kind = 4)		:: limp = .false.	! treat impurity cm coord quantum/classical
	logical (kind = 4)		:: omp_dynamic_enable = .false.	! enable/disable openmp dynamic resource allocation
	integer (kind = 4)		:: nx,ny,nz			! number of points of the grid in each direction
	integer (kind = 4)		:: isalto			! number of lines to skip at reading
	integer (kind = 4)		:: ninvar			! number of components of lambda
	integer (kind = 4)		:: ix, iy, iz, ir			! counters for grid loops
	integer (kind = 4)		:: iter
	integer (kind=4)		:: readmode = 42	! density reading mode 
	integer (kind=4)		:: nthreads = 1		! number of cpu cores to use 
	real (kind = 8)			:: xmax,ymax,zmax! maximum values of x, y and z
	real (kind = 8) 		:: hx,hy,hz			! x, y and z step for the grid
	real (kind = 8), allocatable	:: den(:,:,:)			! helium density
	real (kind = 8), allocatable 	:: x(:),y(:),z(:)		! values in x, y and z
	real (kind = 8)			:: rxden(3)			! vector from density element to impurity
	real (kind = 8)			:: vimp(3)			! vector storing the velocity of the impurity
	real (kind = 8)			:: rimp(3)			! vector storing the position of the impurity
	real (kind = 8)			:: sum,r	! sum of densities, radius between point and impurity
	real (kind = 8)			:: radius
	real (kind = 8)			:: dxyz            	! grid-element volume. (dxyz=hx*hy*hz)
	complex(kind = 8), allocatable	:: invar(:)			! lambda (internal variables)
	complex(kind = 8), allocatable	:: psi(:,:,:)			! wave function. density = mod(psi)**2

	namelist /input/ fileden,nthreads,readmode
	read(5,nml=input)

	!$ call omp_set_dynamic(omp_dynamic_enable)
	!$ call omp_set_num_threads(nthreads)

	select case (readmode)
		case (1)
			open (unit=1, file=fileden)
			call titols(1,cchar,isalto)
			read(1,*) xmax,ymax,zmax,hx,hy,hz,nx,ny,nz,limp,rimp
			allocate (den(nx,ny,nz))
			allocate (x(nx))
			allocate (y(ny))
			allocate (z(nz))
			read(1,*) den
			close(1)
		case (2)
			open (unit=1, file=fileden)
			call titols(1,cchar,isalto)
			read(1,*) xmax,ymax,zmax,hx,hy,hz,nx,ny,nz,limp,rimp
			allocate (psi(nx,ny,nz))
			allocate (den(nx,ny,nz))
			allocate (x(nx))
			allocate (y(ny))
			allocate (z(nz))
			read(1,*) psi
			close(1)
			den = conjg(psi) * psi
		case (3)
			open (unit=1, file=fileden)
			call titols(1,cchar,isalto)
			read(1,*) xmax,ymax,zmax,hx,hy,hz,nx,ny,nz,rimp,vimp
			allocate (psi(nx,ny,nz))
			allocate (den(nx,ny,nz))
			allocate (x(nx))
			allocate (y(ny))
			allocate (z(nz))
			read(1,*) psi
			close(1)
			den = conjg(psi) * psi
		case (4)
			open (unit=1, file=fileden)
			call titols(1,cchar,isalto)
			read(1,*) xmax,ymax,zmax,hx,hy,hz,nx,ny,nz,rimp,vimp,ninvar
			allocate (invar(ninvar))
			allocate (psi(nx,ny,nz))
			allocate (den(nx,ny,nz))
			allocate (x(nx))
			allocate (y(ny))
			allocate (z(nz))
			read(1,*) invar
			read(1,*) psi
			close(1)
			den = conjg(psi) * psi
		case default
			write(*,*)
			write(*,*) "you have chosen a 'denmode' unequal to {1,2,3,4}. please modify '2dden.settings' and choose one of:"
			write(*,*)
			write(*,*) "denmode = 1:	static helium density with/without impurity."
			write(*,*) "denmode = 2:	static helium wave function with/without a vortex/impurity."
			write(*,*) "denmode = 3:	dynamic helium wave function with impurity in an excited 'simple'"
			write(*,*) "		symmetric electronic state, ionised ground state or ground state."
			write(*,*) "denmode = 4:	dynamic helium wave function with impurity in an excited"
			write(*,*) "		anisotropic electronic state / internal state."
			call exit(10)
	end select

	dxyz = hx * hy * hz

	!$omp parallel
	!$omp do private(ix)    
	do ix = 1,nx  !.................... grid x
		x(ix) = -xmax + hx * (ix - 1)
	end do
	!$omp end do nowait

	!$omp do private(iy)    
	do iy = 1,ny  !.................... grid y
		y(iy) = -ymax + hy * (iy - 1)
	end do
	!$omp end do nowait

	!$omp do private(iz)    
	do iz = 1,nz  !.................... grid  z
		z(iz) = -zmax + hz * (iz - 1)
	end do
	!$omp end do
	!$omp end parallel

	do ir = 0,75 ! 30 AA / 0.4 AA = 75
		radius = 0.4*ir
		sum = 0
		!$omp parallel private(ix,iy,iz,rxden,r)
		!$omp do reduction(+:sum)
		do iz = 1,nz
			do iy = 1,ny
				do ix = 1,nx
					rxden(1) = rimp(1) - x(ix)
					rxden(2) = rimp(2) - y(iy)
					rxden(3) = rimp(3) - z(iz)
					r = sqrt(rxden(1)*rxden(1) + rxden(2)*rxden(2) + rxden(3)*rxden(3))
					if (r <= radius) then
						sum = sum+den(ix,iy,iz)
					end if
				end do
			end do
		end do
		!$omp end do
		!$omp end parallel
		sum = sum * dxyz
		write (6,*) radius, sum
	end do
end program


subroutine titols(ulog,cchar,isalto)
	implicit none
	character (len=1) :: pchar
	character (len=1) :: pcolumn
	character (len=1) :: cchar
	integer :: nl
	integer :: i
	integer :: ulog
	integer :: isalto
	nl    = 0
	pcolumn=cchar
	do while(pcolumn.eq.cchar)
	  read(ulog,5000,end=9000) pchar
	  pcolumn=pchar
	  nl=nl+1
	end do
	rewind(ulog)
	isalto=nl-1
	do i=1,isalto
	  read(ulog,*)
	end do
	return
	9000 print *,' ey mister this file is too short ...'  
		 stop 'stop due to severe errors in routine titols'
	5000 format(a1)
end
