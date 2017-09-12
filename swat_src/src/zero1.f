      subroutine zero1

!!    ~ ~ ~ PURPOSE ~ ~ ~
!!    this subroutine initializes the values for some of the arrays 

      use parm
      real :: sub_dgra !! TODO, this should be array like sub_dlag, by lj.
      real :: vp !! seems not used in SWAT, by lj.
      real :: bio_fecal, grwat_veg,plq_rt,rsp_rt,slg_rt
      real :: sub_petco !! seems not used in SWAT, by lj.
      real :: silt  !! do not know what's meaning, by lj.

!!    added for manure Armen Jan 2009
      sol_mc = 0.
      sol_mn = 0.
      sol_mp = 0.
      sol_n = 0.

!!    added by zhang for CSWAT == 2

      sol_BMC = 0.
      sol_BMN = 0.
      sol_HSC = 0.
      sol_HSN = 0.
      sol_HPC = 0.
      sol_HPN = 0.
      sol_LM = 0.
      sol_LMC = 0.
      sol_LMN = 0.
      sol_LS = 0.
      sol_LSL = 0.
      sol_LSC = 0.
      sol_LSN = 0.
      sol_RNMN = 0.
      sol_LSLC = 0.
      sol_LSLNC = 0.
      sol_RSPC = 0.
      sol_WOC = 0.
      sol_WON = 0.
      sol_HP = 0.
      sol_HS = 0.
      sol_BM = 0.
!!    added by zhang for CSWAT == 2
      


!!  septic changes 6/07/10  jaehak
      bio_amn = 0.
      bio_bod = 0.
      fcoli = 0.  
      bio_ntr = 0.
      bio_fecal = 0.  
      bio_ntr = 0.
      biom = 0.
      rbiom = 0.
      bz_perc = 0.
      plqm = 0.
      failyr = 0
      qstemm = 0
      i_sep = 0
      percp = 0
      sep_cap = 0
      bz_area = 0
!!    isep_typ = 1
      isep_typ = 0
      isep_opt= 1 ! fixed by lj, the old code is sep_opt=1
      sep_tsincefail = 0
      isep_tfail = 0
      coeff_bod_dc = 0
      coeff_bod_conv = 0
      coeff_fc1 = 0
      coeff_fc2 = 0
      coeff_fecal = 0
      coeff_plq = 0
      coeff_mrt = 0
      coeff_rsp = 0
      coeff_slg1 = 0
      coeff_slg2 = 0
      coeff_nitr = 0
      coeff_denitr = 0
      coeff_pdistrb = 0
      coeff_psorpmax = 0
      coeff_solpslp = 0
      coeff_solpintc = 0
      sptqs = 0
      sptbodconcs = 0
      spttssconcs = 0
      spttnconcs = 0
      sptnh4concs = 0
      sptno3concs = 0
      sptno2concs = 0
      sptorgnconcs = 0
      spttpconcs = 0
      sptminps = 0
      sptorgps = 0
      sptfcolis = 0
      bz_z = 0
      bz_thk = 0
      bio_bd = 0
      isep_iyr = 0
      isep_opt = 0
      sep_strm_dist = 0
      sep_den = 0
!!  septic changes 6/07/10  jaehak
     
      cont_cn = 0.
      cont_p = 0.
      dep_chan = 0.
      drain_d = 0.
      drain_t = 0.
      drain_g = 0.
      dr_sub = 0.
      filt_w = 0.
      fname = '             '       ! CB 8/24/09
      grwat_n = 0.
      grwat_i = 0.
      grwat_l = 0.
      grwat_w = 0.
      grwat_d = 0.
      grwat_veg = 0.
      gwati = 0.
      gwatn = 0.
      gwatl = 0.
      gwatw = 0. 
      gwatd = 0.
      gwatveg = 0.
      harg_petco = 0.       ! CB 8/24/09
      hru_slp = 0.
      hruaao = 0.
      hruyro = 0.
      icont = 0
      ida_lup = 0           !CB 8/24/09
      ifilt = 0
      ngrwat = 0
      nop = 1
      orig_phu = 0.
      orig_phuacc = 0.
      orig_pltpst = 0.
      orig_pndsed = 0.
      orig_pndvol = 0.
      orig_potno3 = 0.
      orig_potsed = 0.
      orig_potvol = 0.
      orig_sedpstconc = 0.
      orig_shallst = 0.
      orig_snoeb = 0.
      orig_snohru = 0.
      orig_solactp = 0.
      orig_solaorgn = 0.
      orig_solcov = 0.
      orig_solfon = 0.
      orig_solfop = 0.
      orig_solno3 = 0.
      orig_solorgn = 0.
      orig_solorgp = 0.
      orig_solpst = 0.
      orig_solrsd = 0.
      orig_solsolp = 0.
      orig_solstap = 0.
      orig_sumix = 0.
!!    orig_tnylda = 0.
      orig_wetsed = 0.
      orig_wetvol = 0.
      opcp_stat = 0.
      opr_w = 0.
      otmpmx = 0.
      otmpmn = 0.
      otmpstdmn = 0.
      otmpstdmx = 0.
      pcp_stat = 0.
      pperco_sub = 0.
      phu_plt = 0.
      phuacc = 0.
      phug = 300.
      phu_op = 0.
      phutot = 0.
      plaps = 0.
!!  septic changes 1/29/09 
      plqm = 0.
      plq_rt = 0.
!!  septic changes 1/29/09
      plt_pst = 0.
      pname = ""
      pnd_esa = 0.
      pnd_evol = 0.
      pnd_fr = 0.
      pnd_k = 0.
      pnd_nsed = 0.
      pnd_psa = 0.
      pnd_pvol = 0.
      pnd_sed = 0.

      pnd_san = 0.
      pnd_sil = 0.
      pnd_cla = 0.
      pnd_sag = 0.
      pnd_lag = 0.

      pnd_vol = 0.
      pot_fr = 0.
      pot_no3 = 0.
      pot_no3l = 0.
      pot_solpl = 0.
      pot_k = -1.
      pot_nsed = 0.
      pot_sed = 0.

      pot_san = 0.
      pot_sil = 0.
      pot_cla = 0.
      pot_sag = 0.
      pot_lag = 0.

      pot_tile = 0.
      pot_vol = 0.
      pr_w = 0.
      pst_enr = 0.
      pst_kg = 0.
      pst_dep = 0.
      pst_lag = 0.
      pst_wof = 0.
      pst_wsol = 0.
      radinc = 0.
      rch_dakm = 0.
      rchrg = 0.
      rchrg_dp = 0.
      rdmx = 0.
      resouty = 0.
      revapmn = 0.
      rfile = ""
      rfinc = 0.
      rip_fr = 0.
      rk1 = 0.
      rk2 = 0.
      rk3 = 0.
      rk4 = 0.
      rk5 = 0.
      rk6 = 0.
      rnum1s = 0.
!!  septic changes 1/29/09
      rsp_rt = 0.
      slg_rt = 0.
!!  septic changes 1/29/09 
      sol_rock = 0.
      rs1 = 0.
      rs2 = 0.
      rs3 = 0.
      rs4 = 0.
      rs5 = 0.
      rs6 = 0.
      rs7 = 0.
      rsdin = 0.
      rvar_orig = 0.
      sedpst_act = 0.
      sedpst_bry = 0.
      sedpst_conc = 0.
      sedpst_rea = 0.
      sedst = 0.
      shallst = 0.
      shallst_n = 0.
      silt = 0.
      skoc = 0.
      slsoil = 0.
      slsubbsn = 0.
      snam = ""
      sno_hru = 0.
      snoeb = 0.
      sol_actp = 0.
      sol_alb = 0.
      sol_aorgn = 0.
      sol_bd = 0.
      sol_cbn = 0.
      sol_cov = 0.
      sol_crk = 0.
      sol_fon = 0.
      sol_fop = 0.
      sol_hum = 0.
      sol_k = 0.
      sol_kp = 0.
      sol_nh3 = 0.
      sol_nly = 0
      sol_no3 = 0.
      sol_orgn = 0.
      sol_orgp = 0.
      sol_por = 0.
      sol_pst = 0.
      sol_rsd = 0.
      sol_solp = 0.
      sol_stap = 0.
      sol_sumwp = 0.
      sol_up = 0.
      sol_wp = 0.
      sol_z = 0.
      sol_zmx = 0.
      solarav = 0.
      solpstcnst = 0.
      solpstmon = 0.
      solpstyr = 0.
      srbpstcnst = 0.
      srbpstmon = 0.
      srbpstyr = 0.
      strip_n = 0.
      strip_cn = 0.
      strip_c = 0.
      strip_p = 0.
      rchaao = 0.
      rchyro = 0.
      rchmono = 0.
      rch_san = 0.
      rch_sil = 0.
      rch_cla = 0.
      rch_sag = 0.
      rch_lag = 0.
      rch_gra = 0.
      resouta = 0.
      sub_fr = 0.
      sub_km = 0.
      sub_dsan = 0.
      sub_dsil = 0.
      sub_dcla = 0.
      sub_dsag = 0.
      sub_dlag = 0.
      sub_dgra = 0.
      sub_petco = 0.       ! CB 8/24/09
      sub_smfmx = 0.       ! CB 8/24/09
      sub_smfmn = 0.       ! CB 8/24/09
      sub_sftmp = 0.       ! CB 8/24/09
      sub_smtmp = 0.       ! CB 8/24/09
      sub_timp = 0.       ! CB 8/24/09

      subaao = 0.
      subyro = 0.
      sumix = 0.
      sweepeff = 0.
      swtrg = 0
      t_base = 0.
      t_opt = 0.
      t_ov = 0.
      tconc = 0.
      tdrain = 0.
      terr_p = 0.
      terr_cn = 0.
      terr_sl = 0.
      tfile = ""
      thalf = 0.
      tillnm = ""
      title = ""
      tlaps = 0.
      tmp_an = 0.
      tmpinc = 0.
      tmpmn = 0.
      tmpmx = 0.
      tmpstdmn = 0.
      tmpstdmx = 0.
      tnconc = 0.
      tno3conc = 0.
      tnyld = 0.
      tnylda = 0.
      tpconc = 0.
      trapeff = 0.
      urbcoef = 0.
      urbcn2 = 0.
      usle_cfac = 0.
      usle_eifac = 0.
      usle_k = 0.
      usle_ls = 0.
      usle_p = 0.
      vel_chan = 0.
!!  septic changes 1/29/09
      vp = 0. 
!!  septic changes 1/29/09
      vfsi = 0.    !CB 8/24/09
      vpd2 = 0.
      wac21 = 0.
      wac22 = 0.
      wavp = 0.
      wet_fr = 0.
      iwetgw = 0
      iwetile = 0
      wet_k = 0.
      wet_mxsa = 0.
      wet_mxvol = 0.
      wet_nsa = 0.
      wet_nsed = 0.
      wet_nvol = 0.
      wet_sed = 0.

      wet_san = 0.
      wet_sil = 0.
      wet_cla = 0.
      wet_sag = 0.
      wet_lag = 0.

      wet_vol = 0.
      wfsh = 0.
      manure_kg = 0.
      wndav = 0.
      wshd_pstap = 0.

      wshd_sepno3 = 0.
      wshd_sepnh3 = 0.
      wshd_seporgn = 0.
      wshd_sepfon = 0.
      wshd_seporgp = 0.
      wshd_sepfop = 0.
      wshd_sepsolp = 0.
      wshd_sepbod = 0.
      wshd_sepmm = 0.

      wsyf = 0.
      wudeep = 0.
      wurch = 0.
      wushal = 0.
      yldkg = 0.

      depch = 0.
      depfp = 0.
      depprch = 0.
      depprfp = 0.
      depsanfp = 0.
      depsilfp = 0.
      depclafp = 0.
      depsagfp = 0.
      deplagfp = 0.
      depgrafp = 0.
      depsanch = 0.
      depsilch = 0.
      depclach = 0.
      depsagch = 0.
      deplagch = 0.
      depgrach = 0.

      sanst = 0.
      silst = 0.
      clast = 0.
      sagst = 0.
      lagst = 0.
      grast = 0.

      return
      end