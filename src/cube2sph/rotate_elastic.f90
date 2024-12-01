
subroutine rotate_aniso_cijkl(theta,phi,d11,d12,d13,d14,d15,d16, &
    d22,d23,d24,d25,d26, &
    d33,d34,d35,d36,d44,d45,d46,d55,d56,d66, &
    c11,c12,c13,c14,c15,c16,c22,c23,c24,c25,c26, &
    c33,c34,c35,c36,c44,c45,c46,c55,c56,c66)
  implicit none

  double precision,intent(out):: d11,d12,d13,d14,d15,d16,d22,d23,d24,d25,d26, &
                                d33,d34,d35,d36,d44,d45,d46,d55,d56,d66
  double precision,intent(in) ::  c11,c12,c13,c14,c15,c16,c22,c23,c24,c25,c26, &
                                  c33,c34,c35,c36,c44,c45,c46,c55,c56,c66
  double precision, intent(in) :: theta,phi 

  ! local
  double precision,dimension(3,3,3,3) :: C,Cprime
  double precision,dimension(3,3) :: R
  double precision :: cost,cosp,sint,sinp,colat,lon,s  
  double precision,parameter :: pi = atan(1.0) * 4.
  integer,dimension(3,3)      :: vid 

  integer :: i,j,k,l,p,q,m,n
  vid = reshape((/0,5,4,5,1,3,4,3,2/), shape(vid))
  vid = transpose(vid) + 1

  ! construct R tensor
  colat = theta * pi / 180. 
  lon = phi * pi / 180.
  cost = cos(colat); sint = sin(colat)
  cosp = cos(lon); sinp = sin(lon)
  R(1,:) = (/cost * cosp, -sinp, sint * cosp/)
  R(2,:) = (/cost * sinp, cosp, sint * sinp/)
  R(3,:) = (/-sint, dble(0.0), cost/)

  ! copy data to C tensor
  call c_ij2cijkl(c11,c12,c13,c14,c15,c16,c22,c23,c24,c25,c26, &
                  c33,c34,c35,c36,c44,c45,c46,c55,c56,c66,C,.false.)
                  
  ! rotate
  do i=1,3; do j=1,3; do k=1,3; do l=1,3
    s = 0.
    do p=1,3; do q=1,3; do m=1,3; do n = 1,3
      s = s + c(p,q,m,n) * R(i,p) * R(j,q) * R(k,m) * R(l,n)
    enddo;enddo;enddo;enddo;

    Cprime(i,j,k,l) = s
  enddo;enddo;enddo;enddo;

  ! COPY BACK
  call c_ij2cijkl(d11,d12,d13,d14,d15,d16,d22,d23,d24,d25,d26, &
                  d33,d34,d35,d36,d44,d45,d46,d55,d56,d66,Cprime,.true.)


end subroutine rotate_aniso_cijkl

subroutine ij2id(i,j,s)
  implicit none
  integer,intent(in) :: i,j 
  integer :: s 

  integer :: i1,j1 
  i1 = i -1 
  j1 = j - 1
  if(i > j) then 
    i1 = j - 1
    j1 = i - 1
  endif

  s = i1 *3 + j1 - (i1 +1) * i1 / 2 +1
end subroutine ij2id

subroutine c_ij2cijkl(c11,c12,c13,c14,c15,c16,c22,c23,c24,c25,c26, &
  c33,c34,c35,c36,c44,c45,c46,c55,c56,c66,C,inv)
  implicit none
  double precision ::  c11,c12,c13,c14,c15,c16,c22,c23,c24,c25,c26, &
                c33,c34,c35,c36,c44,c45,c46,c55,c56,c66
  double precision,dimension(3,3,3,3) :: C 
  logical,intent(in) :: inv 

  ! local
  double precision :: m(21)
  integer          :: i,j,k,l,i1,j1,id

  ! copy 21 tensor to m21
  m(:) = (/c11,c12,c13,c14,c15,c16,c22,c23,c24,c25,c26, &
  c33,c34,c35,c36,c44,c45,c46,c55,c56,c66/)

  if(.not.inv) then 
    do i=1,3
      do j=1,3
        call ij2id(i,j,i1)
        do k=1,3
          do l=1,3
            call ij2id(k,l,j1)
            call ij2id(i1,j1,id)
            c(i,j,k,l) = m(id)
          enddo
        enddo
      enddo
    enddo
  else
    do i=1,3
      do j=1,3
        call ij2id(i,j,i1)
        do k=1,3
          do l=1,3
            call ij2id(k,l,j1)
            call ij2id(i1,j1,id)
            m(id) = c(i,j,k,l)
          enddo
        enddo
      enddo
    enddo
  endif
end subroutine c_ij2cijkl