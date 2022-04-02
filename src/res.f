      subroutine res
      
!!    ~ ~ ~ PURPOSE ~ ~ ~
!!    this subroutine routes water and sediment through reservoirs
!!    computes evaporation and seepage from the reservoir.

!!    ~ ~ ~ INCOMING VARIABLES ~ ~ ~
!!    name         |units         |definition
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ 
!!    br1(:)       |none          |1st shape parameter for reservoir surface
!!                                |area equation
!!    br2(:)       |none          |2nd shape parameter for reservoir surface
!!                                |area equation
!!    curyr        |none          |current year of simulation
!!    evrsv(:)     |none          |lake evaporation coefficient
!!    iflod1r(:)   |none          |beginning month of non-flood season
!!                                |(needed if IRESCO=2)
!!    iflod2r(:)   |none          |ending month of non-flood season
!!                                |(needed if IRESCO=2)
!!    inum1        |none          |reservoir number
!!    iresco(:)    |none          |outflow simulation code:
!!                                |0 compute outflow for uncontrolled reservoir
!!                                |  with average annual release rate
!!                                |1 measured monthly outflow
!!                                |2 simulated controlled outflow-target release
!!                                |3 measured daily outflow
!!                                |4 stage/volume/outflow relationship 
!!    i_mo         |none          |current month of simulation
!!    ndtargr(:)   |days          |number of days to reach target storage from
!!                                |current reservoir storage
!!                                |(needed if IRESCO=2)
!!    oflowmn(:,:) |m^3/day       |minimum daily ouflow for the month
!!    oflowmx(:,:) |m^3/day       |maximum daily ouflow for the month
!!    pet_day      |mm H2O        |potential evapotranspiration on day
!!    res_evol(:)  |m**3          |volume of water needed to fill the reservoir
!!                                |to the emergency spillway
!!    res_k(:)     |mm/hr         |hydraulic conductivity of the reservoir 
!!                                |bottom
!!    res_nsed(:)  |kg/L          |normal amount of sediment in reservoir
!!    res_pvol(:)  |m**3          |volume of water needed to fill the reservoir
!!                                |to the principal spillway 
!!    res_rr(:)    |m**3/day      |average daily principal spillway release
!!                                |volume
!!    res_sed(:)   |kg/L (ton/m^3)|amount of sediment in reservoir
!!    res_sub(:)   |none          |number of subbasin reservoir is in
!!    res_vol(:)   |m^3 H2O       |reservoir volume
!!    resflwi      |m^3 H2O       |water entering reservoir on day
!!    res_out(:,:,:)|m**3/day      |measured average daily outflow from the
!!                                |reservoir for the month
!!    ressedi      |metric tons   |sediment entering reservoir during time step
!!    sub_subp(:)  |mm H2O        |precipitation for day in subbasin
!!    sub_sumfc(:) |mm H2O        |amount of water in subbasin soil at field 
!!                                |capacity
!!    sub_sw(:)    |mm H2O        |amount of water in soil profile in subbasin
!!    starg(:,:)   |m**3          |monthly target reservoir storage
!!    wuresn(:,:)  |m**3          |average amount of water withdrawn from
!!                                |reservoir each month for consumptive water 
!!                                |use
!!    wurtnf(:)    |none          |fraction of water removed from the reservoir
!!                                |via WURESN which is returned and becomes flow
!!                                |from the reservoir outlet
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ 

!!    ~ ~ ~ OUTGOING VARIABLES ~ ~ ~
!!    name        |units         |definition
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ 
!!    res_sed(:)  |kg/L (ton/m^3)|amount of sediment in reservoir
!!    res_vol(:)  |m^3 H2O       |reservoir volume
!!    resev       |m^3 H2O       |evaporation from reservoir on day
!!    resflwo     |m^3 H2O       |water leaving reservoir on day
!!    respcp      |m^3 H2O       |precipitation on reservoir for day
!!    ressa       |ha            |surface area of reservoir on day
!!    ressep      |m^3 H2O       |seepage from reservoir on day
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ 

!!    ~ ~ ~ LOCAL DEFINITIONS ~ ~ ~
!!    name        |units         |definition
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ 
!!    flw         |m^3/s         |reservoir outflow for day
!!    jres        |none          |reservoir number
!!    sed         |kg/L          |concentration of sediment in reservoir at
!!                               |beginning of day
!!    targ        |m^3 H2O       |target reservoir volume for day
!!    vol         |m^3 H2O       |volume of water stored in reservoir at 
!!                               |beginning of day
!!    vvr         |m^3 H2O       |maximum controlled water release for day
!!    xx          |none          |variable to hold intermediate calculation 
!!                               |result
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ 

!!    ~ ~ ~ SUBROUTINES/FUNCTIONS CALLED ~ ~ ~
!!    Intrinsic: Min


!!    ~ ~ ~ ~ ~ ~ END SPECIFICATIONS ~ ~ ~ ~ ~ ~

      use parm

      integer :: jres
      real*8 :: vol, sed, vvr, targ, xx, flw
      real*8 :: san,sil,cla,sag,lag,gra,ndespill
      real*8 :: inised, finsed, setsed, remsetsed


      !! SJ added all variable definitions below until line 138

      !! vo11, vol2, vol3 are all the volumes in reservoirs 1, 2, and 3 respectively
      !! real_doy is the doy converted to a real number
      !! doy1 and doy2 are the sin and cosine transformed DOYs
      real*8 :: vol1, vol2, vol3, real_doy, doy1, doy2

      !! must define pi with a trig function
      !! evol1, evol2, evol3 are emergency volumes (max volumes) of reservoirs 1, 2, and 3 respectively
      !! pvol1, pvol2, pvol3 are principal spillway volumes (min volumes for a release) of reservoirs 1, 2, and 3 respectively
      real*8 :: pi, evol1, evol2, evol3, pvol1, pvol2, pvol3

      !! IA is the array of reservoir volumes and doy variables (input_array)
      !! RA is the array of optimal releases (release_array)
      real*8, dimension(5) :: IA
      real*8, dimension(3):: RA, f_RA

      !! nvol1, nvol2, and novl3 are normalized reservoir volumes over (0,1) for 1, 2, and 3 respectively
      !! W is the weight array, C is the centers array, R is the radii array (may have to define these elsewhere, depending on where read in)
      !! BF is the exponent value

      real*8:: nvol1, nvol2, nvol3
      real*8:: minrel1, minrel2, minrel3, maxrel1, maxrel2, maxrel3
      integer :: M,N,K,k1,n1,m1
      real*8:: BF

      
      !! dummy variables for interpolation
      !! x and y for calculating reservoir head (hxin, hyin)
      !! x and y for calculating release (rxin, ryin)
      real*8, dimension(:), allocatable :: hxin, hyin, rxin, ryin 

      !! maximum allowable release, minimum allowable release, head 
      real*8 :: mxrel,mnrel,h  

      !! include title of file 
      !! character (len=5) :: c_array, r_array, w_array
      !! character (len=7) :: var_array
      !! character (len=6) :: ra_array
      character (len=7) :: ra_array 
     

      !! define array file names 
      ra_array = "ra.txt"
      N=9
      M=5
      K=3

      !! set jres -- reservoir number
      jres = 0
      jres = inum1


!! store initial values
      vol = 0.
      sed = 0.
      inised = 0.
      finsed = 0.
      setsed = 0.
      remsetsed = 0.
      trapres = 0.

 !! SJ added real_doy, vol1, vol2, vol3, pi, N, M, K, doy1, doy2 initializations
      !! a lot of this should probably be moved into the specific case for DPS policy operation
      !! res_vols skip 2 becasue Gibe II has no storage 
      real_doy = real(iida)
      !! DPS rules don't accomodate leap years 
      if (real_doy == 366.) then
         real_doy = 365.
      end if 

      vol = res_vol(jres)
      sed = res_sed(jres)

      !! get volumes in all 3 reservoirs; 1 = GI, 3 =GIII, 4=Koysha
      vol1 = res_vol(1)
      vol2 = res_vol(3)
      vol3 = res_vol(4)
      pi =4.*ATAN(1.)
      
      !! cyclical representation of time
      doy1 = (sin(2.*pi*(real_doy)/365.) + 1.)/2.
      doy2 = (cos(2.*pi*(real_doy)/365.) + 1.)/2.

      !! Emergency spillway volume & principal spillway volume for GI, GIII, K
      evol1 = res_evol(1)
      evol2 = res_evol(3)
      evol3 = res_evol(4)
      pvol1 = res_pvol(1)
      pvol2 = res_pvol(3)
      pvol3 = res_pvol(4)

      !! SJ
      !! normalize storages on (0,1)
      nvol1 = (vol1 - pvol1) / (evol1 - pvol1)
      nvol2 = (vol2 - pvol2) / (evol2 - pvol2)
      nvol3 = (vol3 - pvol3) / (evol3 - pvol3)

      !! min/max releases for GI, GIII, K
      !! oflowmn/oflowmx give a release volume, so convert to rate 
      minrel1 = oflowmn(i_mo,1) / 86400.
      minrel2 = oflowmn(i_mo,3) / 86400.
      minrel3 = oflowmn(i_mo,4) / 86400.
      maxrel1 = oflowmx(i_mo,1) / 86400.
      maxrel2 = oflowmx(i_mo,3) / 86400.
      maxrel3 = oflowmx(i_mo,4) / 86400.

      !! SJ
      !! create input array
      IA = [nvol1, nvol2, nvol3, doy1, doy2]



    
!!    Storage by particle sizes
      san = res_san(jres)
      sil = res_sil(jres)
      cla = res_cla(jres)
      sag = res_sag(jres)
      lag = res_lag(jres)
      gra = res_gra(jres)

!! calculate surface area for day
      ressa = br1(jres) * res_vol(jres) ** br2(jres)

!! calculate water balance for day
      resev = 10. * evrsv(jres) * pet_day * ressa
      ressep = res_k(jres) * ressa * 240.
      respcp = sub_subp(res_sub(jres)) * ressa * 10.

!! new water volume for day
      if (iresco(jres) /= 5) then 
       res_vol(jres) = res_vol(jres) + respcp + resflwi - resev - ressep
      endif
      
!! subtract consumptive water use from reservoir storage
        xx = 0.
        xx = wuresn(i_mo,jres)
        res_vol(jres) = res_vol(jres) - xx
        if (res_vol(jres) < 0.) then
          xx = xx + res_vol(jres)
          res_vol(jres) = 0.
        end if


!! if reservoir volume is greater than zero

        !! determine reservoir outflow
        select case (iresco(jres))
          case (0)                    !! uncontrolled reservoir
            vvr = 0.
            if (res_vol(jres) > res_pvol(jres)) then
              vvr = res_vol(jres) - res_pvol(jres)
              if (res_rr(jres) > vvr) then
                resflwo = resflwo + vvr
              else
                resflwo = resflwo + res_rr(jres)
              endif
            endif

          case (1)                   !! use measured monthly outflow
            resflwo = res_out(jres,i_mo,curyr)
          !! This will override the measured outflow! This is just a check 
          !! should really calibrate inflow or check res volumes

          case (2)                   !! controlled outflow-target release
            targ = 0.
            if (starg(i_mo,jres) > 0.) then
              targ = starg(i_mo,jres)
            else
              !! target storage based on flood season and soil water
              if (iflod2r(jres) > iflod1r(jres)) then
                if (i_mo > iflod1r(jres) .and. i_mo < iflod2r(jres))    
     &                                                              then
                  targ = res_evol(jres)
                else
                xx = Min(sub_sw(res_sub(jres))/sub_sumfc(res_sub(jres)),
     &                                                               1.)
                targ = res_pvol(jres) + .5 * (1. - xx) *                
     &                                 (res_evol(jres) - res_pvol(jres))
                end if
              else
                if (i_mo > iflod1r(jres) .or. i_mo < iflod2r(jres)) then
                  targ = res_evol(jres)
                else
                xx = Min(sub_sw(res_sub(jres))/sub_sumfc(res_sub(jres)),
     &                                                               1.)
                targ = res_pvol(jres) + .5 * (1. - xx) *                
     &                                 (res_evol(jres) - res_pvol(jres))
                end if
              end if
            endif
            if (res_vol(jres) > targ) then
              resflwo = (res_vol(jres) - targ) / ndtargr(jres)
            else
              resflwo = 0.
            end if

          case (3)                   !! use measured daily outflow
            flw = 0.
            read (350+jres,5000) flw
            resflwo = 86400. * flw
            
          case (4)
            targ = res_pvol(jres) * starg_fps(jres)
            if (res_vol(jres) > targ) then
              resflwo = (res_vol(jres) - targ) / ndtargr(jres)
            else
              resflwo = 0.
            end if
            if (resflwo < oflowmn_fps(jres)) resflwo = oflowmn_fps(jres)
            
          case (5)
            resflwo = 0.
            do jj = 1, nostep
              !! solve quadratic to find new depth
              !testing relationship res_vol(jres) = float(jj) * .1 * res_pvol(jres)
              x1 = bcoef(jres) ** 2 + 4. * ccoef(jres) * (1. - 
     &                                  res_vol(jres) / res_pvol(jres))
              if (x1 < 1.e-6) then
                res_h = 0.
              else
                res_h1 = (-bcoef(jres) - sqrt(x1)) / (2. * ccoef(jres))
                res_h = res_h1 + bcoef(jres)
              end if

              !! calculate water balance for timestep with new surface area
              ressa = res_psa(jres) * (1. + acoef(jres) * res_h)
              resev = 10. * evrsv(jres) * pet_day * ressa
              ressep = res_k(jres) * ressa * 240.
              respcp = sub_subp(res_sub(jres)) * ressa * 10.

              if(res_h <= 1.e-6) then
                res_qi = 0.
                res_h = 0.
              else
                res_qi = weirc(jres) * weirk(jres) * weirw(jres) * 
     &                                                    (res_h ** 1.5)
              end  if
              resflwo = resflwo + res_qi
              res_vol(jres) = res_vol(jres) + (respcp + resflwi - resev 
     &                                                - ressep) / nostep
              res_vol(jres) = res_vol(jres) - res_qi
            enddo

            !! SJ defines case 6
            !! DPS optimized release
            case(6)
            !! only calculate release matrix for a given day for reservoir 1
            !! then pull the rest of the values from the matrix 

                !! base release
                do k1 = 1,K
                  RA(k1) = var(k1)
                end do

                !! calculate release matrix 
                do k1 = 1,K
                  do n1 = 1,N
                    BF = 0.
                    do m1 = 1,M
                      if (R_arr(m1,n1) > ((10.)**(-6.))) then
                        BF = BF + ((IA(m1)-C_arr(m1,n1))/R_arr(m1,n1))**2.
                        !! if (jres == 1) then
                        !! print*, IA(m1), C(m1,n1), R(m1,n1)
                        !!  print*, BF
                        !! end if
                      else
                        BF = BF + ((IA(m1)-C_arr(m1,n1))/((10.)**(-6.)) )**2.
                      endif
                    end do
                    RA(k1) = RA(k1) + W_arr(k1,n1)*exp(BF * -1.)
                   !! if (jres == 1) then
                      !! print*, RA(k1), W(k1,n1), BF, exp(BF * -1.)
                   !! end if

                    if (RA(k1) < 0.) then
                      RA(k1) = 0.
                    elseif (RA(k1) > 1.) then
                      RA(k1) = 1.
                    end if
                  end do
                end do

                if (jres == 1) then
                  open(405, file=ra_array)
                  rewind(405)
                  write(405,*) (RA, col=1,K)
                  rewind(405)
                end if

              !! read in RA values
              open(405, file=ra_array)
              read(405, *) (f_RA, col=1,K) 
              rewind(405)

              !! define resflwo -- if jres > 1, then index is jres minus 1 because reservoir 2 (Gibe II) not included 
              !! de-normalize on (0,1)

              if (jres==1) THEN
                resflwo = minrel1 + f_RA(jres)*(maxrel1 - minrel1)
              elseif (jres==3) then
                resflwo = minrel2 +f_RA(jres-1)*(maxrel2 - minrel2)
              elseif (jres==4) then
                resflwo = minrel3 +f_RA(jres-1)*(maxrel3 - minrel3)
              endif

              !! get full volume from rate
              resflwo = resflwo * 86400.



        end select
          
            !! ndespill = ndtargr(jres)
            !! if (ndespill <= 0.) ndespill = 10.
            !! only release donw to the principal spillway volume 
            !! if ((res_vol(jres) - resflwo) < res_pvol(jres)) then
            !!  resflwo = res_vol(jres) - res_pvol(jres)
            !! endif  
            !! if volume is too great, release all of it 
            !! if (res_vol(jres) > res_evol(jres)) then
            !!  resflwo = resflwo+(res_vol(jres)-res_evol(jres))/ndespill
            !! endif


            !! vvol = res_vol(jres)!! dummy vble 
            h = 0. !! initialize h 
            mxrel = 0.
            mnrel = 0.
            !! SJ ADDED: check if release is valid 
            !! Gibe I 
            if (jres == 1) THEN
              !! calculate head
              !! row 1 is head row 3 is storage 
              hyin = gi_rc(1,:)
              hxin = gi_rc(3,:)              
              h = interp_lin(hxin, hyin, vol) !! vvol)

              !! calculate max allowable release
              rxin = gi_max_rel(1,:)
              ryin = gi_max_rel(2,:)
              mxrel = interp_lin(rxin, ryin, h)

              !! calculate min allowable release 
              !! if (h > 37.) then
              !!   mnrel = mxrel
              !! elseif (h >= 36.) then 
              !!   h = 37. 
              !!   mnrel = interp_lin(rxin, ryin, h)
              !!   mnrel = (mnrel - 1.3) * (h - 36) + 1.3
              if (h < 18.) then 
                mnrel = 0. 
              else
                mnrel = 1.3  !! MEF
              end if

            !! GIBE III
            elseif (jres == 3) THEN
              !! calculate head 
              hxin = giii_rc(3,:)
              hyin = giii_rc(1,:)
              h = interp_lin(hxin, hyin, vol) !! vvol)

              !! max allowable release
              rxin = giii_max_rel(1,:)
              ryin = giii_max_rel(2,:)
              mxrel = interp_lin(rxin, ryin, h)

              !! min allowable release
              if (h < 140.) then
                mnrel = 0.
              !! elseif (h >= 232.) then 
                !! mnrel = mxrel
                !! mnrel = (vol - 14700000000.) / 86400.
              else
                mnrel = 70. !! MEF
              end if

            !! Koysha 
            elseif (jres == 4) THEN
              !! calculate head 
              hxin = k_rc(3,:)
              hyin = k_rc(1,:)
              h = interp_lin(hxin,hyin,vol)!! vvol)

              !! calculate max allowable release
              rxin = k_max_rel(1,:)
              ryin = k_max_rel(2,:)
              mxrel = interp_lin(rxin, ryin, h)

              !! calc min allowable relase
              if (h < 106.) then
                mnrel = 0.
              !! elseif (h >= 167.) then 
                !! mnrel = mxrel
                !! mnrel = (vol - 5700000000.) / 86400. 
              else
                mnrel = 25.  !! MEF 
              end if
            end if

            !! SJ -- convert from rate to volume
            !! mxrel = mxrel * 84600. 
            !! mnrel = mnrel * 84600.

            !! ensure releases are valid 
            if (resflwo / 86400. > mxrel) then 
              resflwo = mxrel * 86400.
            elseif (resflwo / 86400. < mnrel) then
              resflwo = mnrel * 86400. 
            end if

            ndespill = ndtargr(jres)
            if (ndespill <= 0.) ndespill = 10.
            !! only release donw to the principal spillway volume 
            if ((res_vol(jres) - resflwo) < res_pvol(jres)) then
              resflwo = res_vol(jres) - res_pvol(jres)
            endif  

            !! ASK ABOUT THIS 
            !! if volume is too great, release all of it 
            !! if (res_vol(jres) > res_evol(jres)) then
              !! resflwo = resflwo + (res_vol(jres)-res_evol(jres))/ndespill --> FIRST
            !!  resflwo = res_vol(jres) - res_evol(jres) -- > SECOND
            !! THIRD: if not releaseing enough, release down to full
            if ((res_vol(jres) - resflwo) > res_evol(jres)) then
              resflwo = res_vol(jres) - res_evol(jres)
            endif


!! if reservoir volume is zero
      if (res_vol(jres) < 0.001) then

        !! if volume deficit in reservoir exists, reduce seepage so
        !! that reservoir volume is zero
        ressep = ressep + res_vol(jres)
        res_vol(jres) = 0.

        !! if seepage is less than volume deficit, take remainder
        !! from evaporation
        if (ressep < 0.) then
          resev = resev + ressep
          ressep = 0.
        end if
        res_sed(jres) = 0.

      else

        !! check calculated outflow against specified max and min values
        if (resflwo < oflowmn(i_mo,jres)) resflwo = oflowmn(i_mo,jres)
        if (resflwo > oflowmx(i_mo,jres) .and. oflowmx(i_mo,jres) > 0.) 
     &                                                              then
          resflwo = oflowmx(i_mo,jres)
        endif
           
        !! subtract outflow from reservoir storage
        if(iresco(jres) /= 5) then
          res_vol(jres) = res_vol(jres) - resflwo
          if (res_vol(jres) < 0.) then
             resflwo = resflwo + res_vol(jres)
             res_vol(jres) = 0.
          end if
        end if  

        !! add spillage from consumptive water use to reservoir outflow
        resflwo = resflwo + xx * wurtnf(jres)

        !! compute new sediment concentration in reservoir
        if (ressedi < 1.e-6) ressedi = 0.0      !!nbs 02/05/07
        if (ressa == 0.) ressa = 1.e-6     !! MJW added 040711
        velofl = (resflwo / ressa) / 10000.  !!m3/d / ha * 10000. = m/d
        !!    velsetl = 1.35      !! for clay particle m/d
        if (velofl > 1.e-6) then
          trapres = velsetlr(jres) / velofl
        if (trapres > 1.) trapres = 1.  !! set to nres
          susp = 1. - trapres
        else
          susp = 1.
        end if

        if (res_vol(jres) > 0.) then                         !!MJW added 040811
        res_sed(jres) = (ressedi * susp + sed * vol) / res_vol(jres)
        res_san(jres) = (ressani + san * vol) / res_vol(jres)
        res_sil(jres) = (ressili + sil * vol) / res_vol(jres)
        res_cla(jres) = (resclai + cla * vol) / res_vol(jres)
        res_sag(jres) = (ressagi + sag * vol) / res_vol(jres)
        res_lag(jres) = (reslagi + lag * vol) / res_vol(jres)
        res_gra(jres) = (resgrai + gra * vol) / res_vol(jres)

        res_sed(jres) = Max(1.e-6,res_sed(jres))
        res_san(jres) = Max(1.e-6,res_san(jres))
        res_sil(jres) = Max(1.e-6,res_sil(jres))
        res_cla(jres) = Max(1.e-6,res_cla(jres))
        res_sag(jres) = Max(1.e-6,res_sag(jres))
        res_lag(jres) = Max(1.e-6,res_lag(jres))
        res_gra(jres) = Max(1.e-6,res_gra(jres))
        else
        res_sed(jres) = 1.e-6             !!MJW added 040711
        res_san(jres) = 1.e-6
        res_sil(jres) = 1.e-6
        res_cla(jres) = 1.e-6
        res_sag(jres) = 1.e-6
        res_lag(jres) = 1.e-6
        res_gra(jres) = 1.e-6
        endif
        
        !! compute change in sediment concentration due to settling
        if (res_sed(jres) < 1.e-6) res_sed(jres) = 0.0    !!nbs 02/05/07
        if (res_sed(jres) > res_nsed(jres)) then
      inised = res_sed(jres)
          res_sed(jres) = (res_sed(jres) - res_nsed(jres)) *            
     &                                   sed_stlr(jres) + res_nsed(jres)
      finsed = res_sed(jres)
      setsed = inised - finsed

        if (res_gra(jres) >= setsed) then
      res_gra(jres) = res_gra(jres) - setsed
      remsetsed = 0.
      else
      remsetsed = setsed - res_gra(jres)
          res_gra(jres) = 0.
      if (res_lag(jres) >= remsetsed) then
        res_lag(jres) = res_lag(jres) - remsetsed
        remsetsed = 0.
      else
        remsetsed = remsetsed - res_lag(jres)
        res_lag(jres) = 0.
        if (res_san(jres) >= remsetsed) then
          res_san(jres) = res_san(jres) - remsetsed
          remsetsed = 0.
        else
          remsetsed = remsetsed - res_san(jres)
          res_san(jres) = 0.
              if (res_sag(jres) >= remsetsed) then
            res_sag(jres) = res_sag(jres) - remsetsed
            remsetsed = 0.
          else
            remsetsed = remsetsed - res_sag(jres)
            res_sag(jres) = 0.
                if (res_sil(jres) >= remsetsed) then
                res_sil(jres) = res_sil(jres) - remsetsed
              remsetsed = 0.
            else
              remsetsed = remsetsed - res_sil(jres)
              res_sil(jres) = 0.
                  if (res_cla(jres) >= remsetsed) then
                res_cla(jres) = res_cla(jres) - remsetsed
                remsetsed = 0.
              else
                remsetsed = remsetsed - res_cla(jres)
                res_cla(jres) = 0.
              end if
                end if
          end if
        end if
        end if
        endif

        end if

        !! compute sediment leaving reservoir
        ressedo = res_sed(jres) * resflwo
        ressano = res_san(jres) * resflwo
        ressilo = res_sil(jres) * resflwo
        resclao = res_cla(jres) * resflwo
        ressago = res_sag(jres) * resflwo
        reslago = res_lag(jres) * resflwo
        resgrao = res_gra(jres) * resflwo

        !! net change in amount of sediment in reservoir for day
        ressedc = vol * sed + ressedi - ressedo - res_sed(jres) *       
     &                                                     res_vol(jres)
!      write (130,5999) i, jres, res_sed(jres), sed_stlr(jres),          
!     & res_nsed(jres), ressedi, ressedo, resflwi, resflwo
!5999  format (2i4,7e12.4)
      
      end if

!!    update surface area for day
      ressa = br1(jres) * res_vol(jres) ** br2(jres)

     !!  print*, iida, curyr, jres, resflwo
      return
 5000 format (f8.2)
      end