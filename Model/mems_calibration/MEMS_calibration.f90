Module MEMS_calibration
  implicit none

  contains

subroutine MEMS(inPara,input,soil,LQ,C)
implicit none

double precision, dimension(30), intent(in) :: inPara
double precision, dimension(365,4), intent(in) :: input
double precision, dimension(7), intent(in) :: soil
double precision, dimension(4), intent(in) :: LQ
double precision, dimension(11), intent(out) :: C

integer :: ii, jj, kk
integer :: Nyear, time
double precision :: NPP, tmp, pi
double precision :: mean_temp, tmp_range, hemis
double precision :: fsi, fdoc, flig
double precision :: LCI, LitN, tmp_mod
double precision :: tmod_opt, tmod_lag, tmod_shp, tmod_Q10, tmod_Tref
double precision, dimension(4) :: litter
double precision, dimension(11) :: dC

pi = 3.1415927

Nyear = 700
time = 1

C = 0.
dC = 0.

C(1) = 0.35
C(2) = 0.35
C(3) = 0.15
C(6) = 0.15
C(9) = 3000.

! Average temperature and hemisphere value
hemis = -1.5 ! Nothern latitude

mean_temp = sum((input(:,1)+input(:,2))/2)/365.

! Litter fractions
LitN = LQ(1)
fsi = LQ(2)
flig = LQ(3)
fdoc = LQ(4)

jj = 0
LCI = 0.

! Values from the manual
tmod_opt = 45.
tmod_lag = 4.
tmod_shp = 15.
tmod_Q10 = 2.
tmod_Tref = 13.5

do ii = 1, Nyear*365

   jj = jj + 1

   if(jj .gt. 365) then
      jj = 1
   end if

   C(7) = 0.

   NPP = input(jj,4)
   litter(1) = NPP*fsi - NPP*fsi*fdoc
   litter(2) = NPP - NPP*(fsi+flig)
   litter(3) = NPP*flig
   litter(4) = NPP*fsi*fdoc

   tmp_range = (input(jj,1) - input(jj,2))/2.
   tmp = tmp_range*sin((2*jj/365.*pi)+hemis) + mean_temp

   tmp_mod = exp(-(tmp/(tmod_opt+tmod_lag))**tmod_shp) &
        *(tmod_Q10**((tmp-tmod_Tref)/tmod_Tref))

   call MEMS_func(inPara,C,litter,LitN,tmp_mod,soil,dC)
   C = C + dC*time

   do kk = 1, 11
      if(C(kk) .lt. 0) then
         C(kk) = 0.
      end if
   end do

end do

end subroutine MEMS


subroutine MEMS_func(inPara,C,litter,LitN,tmp_mod,soil,dC)
implicit none

double precision, dimension(30), intent(in) :: inPara ! Parameters
double precision, dimension(11), intent(in) :: C   ! Carbon pool values
double precision, dimension(4), intent(in) :: litter ! Incoming litter
double precision, intent(in) :: LitN ! Nitrogen content of litter
double precision, intent(in) :: tmp_mod ! Temperature decay rate modifier
double precision, dimension(7), intent(in) :: soil ! Soil properties
double precision, dimension(11), intent(out) :: dC ! Carbon pool change rates

integer :: ii
double precision :: LCI, ek, eb, uk, ub, gamma
double precision :: la1, la4
double precision :: Klm, scConc, Qmax, sorption
double precision :: POMsplit
double precision, dimension(25) :: Para

Para(1) = inPara(9)  ! k1
Para(2) = inPara(10) ! k2
Para(3) = 0.         ! k3, calculated during the loops
Para(4) = inPara(11) ! k4
Para(5) = 0.0005 ! inPara(2)*exp(-2.1) ! k5 = k2*exp(-2.1)
Para(6) = inPara(23) ! k8
Para(7) = inPara(24) ! k9
Para(8) = inPara(25) ! k10
Para(9) = inPara(6)  ! Nmax
Para(10) = inPara(7) ! Nmid
Para(11) = inPara(8) ! LCmax
Para(12) = inPara(29) ! LITfrg
Para(13) = inPara(28) ! DOCfrg
Para(14) = inPara(12) ! EHmin
Para(15) = inPara(13) ! EHmax
Para(16) = inPara(14) ! ESmin
Para(17) = inPara(15) ! ESmax
Para(18) = inPara(18) ! B1
Para(19) = inPara(19) ! B2
Para(20) = inPara(20) ! B3
Para(21) = inPara(16) ! la2
Para(22) = inPara(17) ! la3
Para(23) = inPara(1)  ! scIcept
Para(24) = inPara(2)  ! scSlope
Para(25) = inPara(22) ! dlch

POMsplit = 0.3

if(C(3) .gt. 0.) then
   LCI = C(3)/(C(2)+C(3))
else
   LCI = 0.
end if

gamma = 1/(1+exp(-Para(9)*(LitN-Para(10))))

! Decay rate for pool 3 and the Decay rate modifier based on lignin
if(LCI .lt. 0.7) then
   Para(3) = Para(2)*(0.2/(1+(200/exp(8.15*LCI))))
   ek = exp(-3*LCI)
else
   Para(3) = Para(2)*exp(-2.1)
   ek = exp(-3*0.7)
end if

! Microbial uptake modifier based on lignin
if(LCI .lt. 0.7) then
   eb = 1-exp(-0.7*(0.7-LCI)*10)
else
   eb = 0.
end if

uk = min(gamma,ek)
ub = min(gamma,eb)

! DOC generation from leaching
la1 = min((Para(15)-((Para(15)-Para(14))/Para(11))*LCI), &
     (Para(15)-((Para(15)-Para(14))/Para(9))*LitN))

la4 = min((Para(17)-((Para(17)-Para(16))/Para(11))*LCI), &
     (Para(17)-((Para(17)-Para(16))/Para(9))*LitN))

if(la1 .lt. 0) then
   la1 = 0.
end if

if(la4 .lt. 0) then
   la4 = 0.
end if

! Calculate the sorption and binding capacity of the soil

Klm = exp(-0.186*soil(4)-0.216) ! Binding affinity
ScConc = Para(24)*(100-soil(1)) + Para(23) ! Maximum clay+silt concentration for these processes
Qmax = soil(7)*100*100*soil(6)/1000*scConc*(1-soil(3)) ! Maximum sorption
sorption = C(8)*((Klm*Qmax*C(8))/(1+Klm*C(8))-C(9))/Qmax

dC = 0.

dC(1) = litter(1) - Para(1)*uk*tmp_mod*C(1)
dC(2) = litter(2) - Para(2)*uk*tmp_mod*C(2) - Para(12)*C(2)
dC(3) = litter(3) - Para(3)*tmp_mod*C(3) - Para(12)*C(3)
dC(4) = ub*Para(19)*(1-la1)*Para(2)*uk*tmp_mod*C(2) &
     + ub*Para(18)*(1-la4)*Para(1)*uk*tmp_mod*C(1)  &
     - Para(4)*tmp_mod*C(4)
dC(5) = (Para(12)*C(2) + Para(12)*C(3))*POMsplit &
     + Para(20)*(1-Para(21))*Para(4)*tmp_mod*C(4) &
     - Para(5)*tmp_mod*C(5)
dC(6) = litter(4) + la1*Para(2)*uk*tmp_mod*C(2) &
     + Para(22)*Para(3)*tmp_mod*C(3) &
     + Para(21)*Para(4)*tmp_mod*C(4) &
     + la4*Para(1)*uk*tmp_mod*C(1) &
     - Para(13)*C(6)
dC(7) = ((1-ub*Para(18))*(1-la4))*Para(1)*uk*tmp_mod*C(1) &
     + ((1-ub*Para(19))*(1-la1))*Para(2)*uk*tmp_mod*C(2) &
     + (1-Para(22))*Para(3)*tmp_mod*C(3) &
     + ((1-Para(20))*(1-Para(21)))*Para(4)*tmp_mod*C(4) &
     + (1-Para(22))*Para(5)*tmp_mod*C(5) &
     + Para(6)*tmp_mod*C(8) + Para(7)*tmp_mod*C(9) &
     + (1-Para(22))*Para(8)*tmp_mod*C(10)
dC(8) = C(6)*Para(13) &
     + Para(22)*Para(5)*C(5)*tmp_mod &
     + Para(22)*Para(8)*C(10)*tmp_mod &
     - sorption &
     - Para(25)*C(8) &
     - Para(6)*tmp_mod*C(8)
dC(9) = sorption &
     - Para(7)*tmp_mod*C(9)
dC(10) = (1-POMsplit)*Para(12)*C(2) &
     + (1-POMsplit)*Para(12)*C(3) &
     - Para(8)*tmp_mod*C(10)
dC(11) = Para(25)*C(8)

end subroutine MEMS_func

end Module MEMS_calibration
