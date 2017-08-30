PROGRAM ROTMOL

	implicit none
	integer i, j, k, natoms, nop
	logical :: dbg = .false.
	real*8 pi
	parameter(pi = 3.14159265359)
	character*3, allocatable :: atom(:)
	integer, allocatable :: nstep(:)
	real*8, allocatable :: R(:,:), XYZ(:,:), NEWXYZ(:,:)
	real*8, allocatable :: phi(:), theta(:), gama(:)
	
	! Ler arquivo "input"
	call readInput()
	
	!abrir "output.xyz"
	open(unit=2, file="output.xyz")	
	if(dbg) then
		open(unit=99, file="log")
		call writeParam()
	end if
	
	
	do i=1,nop
	
		R = RMAT(phi(i),theta(i),gama(i))
		do j=1,nstep(i)
		
			do k=1,natoms
				NEWXYZ(k,:) = MULTVEC(XYZ(k,:),R)
			end do
			
			XYZ = NEWXYZ
		end do
	end do
	
	
	
	
	
	print *, "----------matrix de rotacao---------"
	
	do i=1,3
		print "(3F9.5)", (R(i,j), j=1,3)
	end do
	
	print *, "--------posicoes iniciais--------------"
	
	do i=1,natoms
		print "(3F9.5)", (XYZ(i,j), j=1,3)
	end do
	
	
	print *, "------novas posicoes----------------"
	
	do i=1,natoms
		print "(3F9.5)", (NEWXYZ(i,j), j=1,3)
	end do

	contains


real*8 function RMAT(phi,theta,gama)

! FUNCAO QUE RETORNA A OPERACAO (MATRIZ) DE ROTACAO RESULTANTE A PARTIR DAS TRES
! OPERACOES PRIMARIAS DE ROTACAO EM PHI(X) THETA(Y) E GAMMA(Z)

	implicit none
	real*8, intent(in) :: phi,theta,gama
	real*8 :: RxMat, RyMat, RzMat, c, s
	dimension RxMat(3,3), RyMat(3,3), RzMat(3,3), RMAT(3,3)
	
	
	RxMat(1,:) = (/ 1.0d0, 0.0d0, 0.0d0 /)
	RxMat(2,:) = (/ 0.0d0, cos(phi), -sin(phi) /)
	RxMat(3,:) = (/ 0.0d0, sin(phi), cos(phi) /)
		
	RyMat(1,:) = (/ cos(theta), 0.0d0, sin(theta) /)
	RyMat(2,:) = (/ 0.0d0, 1.0d0, 0.0d0 /)
	RyMat(3,:) = (/ -sin(theta), 0.0d0, cos(theta) /)
	
	RzMat(1,:) = (/ cos(gama), -sin(gama), 0.0d0 /)
	RzMat(2,:) = (/ sin(gama), cos(gama), 0.0d0 /)
	RzMat(3,:) = (/ 0.0d0, 0.0d0, 1.0d0 /)
	
	RMAT = MULTMAT(RzMat,RyMat)
	RMAT = MULTMAT(RMAT,RxMat)
	

end function RMAT

real*8 function MULTMAT(A,B)

!	FUNCAO PARA MULTIPLICACAO DE DUAS MATRIZES(3,3)

	implicit none
	real*8, dimension(3,3), intent(in) :: A, B
	dimension MULTMAT(3,3)
	integer :: i, j, k
	
	MULTMAT = 0.d0
	do i=1,3
		do j=1,3
			do k=1,3
				MULTMAT(i,j) = MULTMAT(i,j) + A(i,k)*B(k,j)
			end do
		end do
	end do
	
end function MULTMAT

real*8 function MULTVEC(V,A)
	
! FUNCAO DE MULTIPLICACAO DE UM VETOR(3) POR UMA MATRIZ(3,3)
	
	implicit none
	real*8, dimension(3,3), intent(in) :: A
	real*8, dimension(3), intent(in) :: V
	dimension MULTVEC(3)
	integer :: i, j

	MULTVEC = 0.d0
	
	do i=1,3
		do j=1,3
			MULTVEC(i) = MULTVEC(i) + V(j)*A(i,j)
		end do
	end do


end function MULTVEC

subroutine readInput()

	open(unit=1, file="input")
	read(1,*) natoms, nop
	if(natoms.lt.0) then
		natoms = -natoms
		dbg = .true.
	end if
	ALLOCATE(XYZ(1:natoms,1:3),ATOM(1:natoms),NEWXYZ(1:natoms,1:3))
	ALLOCATE(nstep(nop),phi(nop),theta(nop),gama(nop))
	
	do i=1,nop
		read(1,*) nstep(i), phi(i), theta(i), gama(i)
		phi(i) = phi(i)*pi/180
		theta(i) = theta(i)*pi/180
		gama(i) = gama(i)*pi/180
	end do
	
	do i=1,natoms
		read(1,*) atom(i), (XYZ(i,j), j=1,3)
	end do
	
	close(1)

end subroutine readInput

subroutine writeParam()
	write(99,*) "-------ROTMOL---DEBUG---------------"
	write(99,*) "----------mcquintao-----------------"
	write(99,*) ""
	write(99,*) "Numero de atomos e operaçoes:"
	write(99,"(2I4)") natoms, nop

	write(99,*) "Operaçoes: (passos, phi theta gama)"
	do i=1,nop
		write(99,"(I3, 3F9.3)") nstep(i), phi(i)*180/pi, theta(i)*180/pi, gama(i)*180/pi
	end do
	write(99,*) "Coordenadas iniciais:"
	do i=1,natoms
		write(99,"(A3, 3F8.4)") atom(i), (XYZ(i,j), j=1,3)
	end do
	write(99,*) "-----------------------------------"
	write(99,*) ""
end subroutine

END PROGRAM ROTMOL
