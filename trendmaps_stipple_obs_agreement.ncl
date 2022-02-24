;Commonly loaded

load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/contrib/calendar_decode2.ncl"

begin

modelversion = (/"UKESM1","ukesm1p1"/)
modelname_long = (/"UKESM1","UKESM1.1"/)
plotname = "overlap_w_obs_temp_trends_"+modelversion(1)
wks=gsn_open_wks("png",plotname)
plot = new(9,graphic)
rmse = new(9,float)
nyear=165
syear = (/1901,1941,1981/)
eyear = (/1940,1980,2014/)

do i = 0,2
  
  colormap = read_colormap_file("MPL_RdBu")
  seas = "ann"

  runid = (/"1","2","3","4","5","6","7","8","9","10","11","12","16","17","18","19"/)


; =================== control run ukesm1 ===================                                                                                              

  branch_off2 = (/2121,2165,2210,2250,2285,2340/)
  inctrl2 = addfile("tas_aw310.nc","r")
  ctrl2 = inctrl2->tas
  lonc = inctrl2->lon
  if (lonc(0) .ge. 0) then
    ctrl2 = lonFlip(ctrl2)
  end if
  ctrl2!1 = "latitude"
  ctrl2!2 = "longitude"
  ann_ctrl2 = month_to_annual(ctrl2(time|0:,latitude|:,longitude|:),1)
  ann_ctrl2&year = ispan(1960,3839,1)
  
  drift_trend2 = new((/6,165,dimsizes(ann_ctrl2&latitude),dimsizes(ann_ctrl2&longitude)/),float)
  drift_trend2!0 = "ncl_join"
  drift_trend2!1 = "year"
  drift_trend2!2 = "latitude"
  drift_trend2!3 = "longitude"
;  slope = new((/6,dimsizes(ann_ctrl2&latitude),dimsizes(ann_ctrl2&longitude)/),float)
  test = ispan(branch_off2(0),branch_off2(0)+164,1)
;  printVarSummary(test)

  do c = 0,5
    slope := regCoef_n(ispan(branch_off2(c),branch_off2(c)+164,1),ann_ctrl2({year|branch_off2(c):(branch_off2(c)+164)},latitude|:,longitude|:),0,0)
    intercept := slope@yintercept
;    printVarSummary(intercept)                                                                                                                           
    drift_trend2(c,:,:,:) = conform(drift_trend2(c,:,:,:),slope,(/1,2/))*conform(drift_trend2(c,:,:,:),ispan(branch_off2(c),branch_off2(c)+164,1),0)+conform(drift_trend2(c,:,:,:),onedtond(slope@yintercept,(/dimsizes(ann_ctrl2&latitude),dimsizes(ann_ctrl2&longitude)/)),(/1,2/))
  end do
;  copy_VarCoords(ann_ctrl2(0,:,:),slope(0,:,:))
; =================== read in model data ============

  files := systemfunc("ls -1 UKESM1/historical/tas*nc")
  in = addfiles(files,"r")
  ListSetType(in,"join")

  paper_member_ind = (/4,3,1,2,5,0/) ; indices for member selection in SO2-dep paper
  ens := in[paper_member_ind]->tas(:,:,:,:)
  print(files(paper_member_ind))
  print(branch_off2)
  units = ens@units
  printVarSummary(ens)
  ens!2="latitude"
  ens!3="longitude"
  lonm = ens&longitude
 if (lonm(0) .ge. 0) then
   ens = lonFlip(ens)
 end if

  ann1p0 := month_to_annual(ens(time|:,ncl_join|:,latitude|:,longitude|:),1)
  DJF1p0 := month_to_season(ens(time|:,ncl_join|:,latitude|:,longitude|:),"DJF")
  JJA1p0 := month_to_season(ens(time|:,ncl_join|:,latitude|:,longitude|:),"JJA")

  ann1p0!0 = "year"
  printVarSummary(ann1p0)
  dims = dimsizes(ann1p0)
  nmembers = dims(1)
  ann1p0&year = ispan(1850,2014,1)
  DJF1p0&time = ispan(1850,2014,1)
  JJA1p0&time = ispan(1850,2014,1)
  
;  ann1p0_anom = ann1p0 - conform(ann1p0,dim_avg_n_Wrap(ctrl2,0),(/2,3/)

  ddrifted := ann1p0(year|:,ncl_join|:,latitude|:,longitude|:)-drift_trend2(year|:,ncl_join|:,latitude|:,longitude|:)
  ddrifted!0 = "year"
  ddrifted!1 = "ncl_join"
  ddrifted!2 = "latitude"
  ddrifted!3 = "longitude"
  copy_VarCoords(ann1p0(year|:,ncl_join|:,latitude|:,longitude|:),ddrifted(year|:,ncl_join|:,latitude|:,longitude|:))
  print(ann1p0(year|164,ncl_join|2,latitude|25,longitude|25))
  ann1p0 := ddrifted
  print(drift_trend2(year|164,ncl_join|2,latitude|25,longitude|25))
  print(ddrifted(ncl_join|2,year|164,latitude|25,longitude|25))

  tval1 = new(dimsizes(ann1p0(0,0,:,:)),double)
  nptxy1 = new(dimsizes(ann1p0(0,0,:,:)),integer)
  trends_1p0 = regcoef(ispan(syear(i),eyear(i),1),dim_avg_n_Wrap(ann1p0(ncl_join|:,latitude|:,longitude|:,{year|syear(i):eyear(i)}),0),tval1,nptxy1)*10
  copy_VarCoords(ann1p0(0,0,:,:),trends_1p0)

  trends_1p0_all = regCoef(ispan(syear(i),eyear(i),1),ann1p0(ncl_join|:,latitude|:,longitude|:,{year|syear(i):eyear(i)}))*10
  copy_VarCoords(ann1p0(ncl_join|:,year|0,latitude|:,longitude|:),trends_1p0_all)
  out = addfile("trends_1p0.nc","c")
  out->trends_1p0 = trends_1p0
  
; =================== control run ukesm1p1 ===================                                                                                              

  branch_off = (/2811,2851,2891,2931,2971,3011/)
  inctrl = addfile("tas_by230.nc","r")
  ctrl = inctrl->UM_m01s03i236_vn1107
  lonc := inctrl->longitude
  if (lonc(0) .ge. 0) then
    ctrl = lonFlip(ctrl)
  end if

  ann_ctrl = month_to_annual(ctrl(time|0:5555,latitude|:,longitude|:),1)
  ann_ctrl&year = ispan(2742,3204,1)
  drift_trend = new((/6,165,dimsizes(ann_ctrl&latitude),dimsizes(ann_ctrl&longitude)/),float)
  drift_trend!0 = "ncl_join"
  drift_trend!1 = "year"
  drift_trend!2 = "latitude"
  drift_trend!3 = "longitude"
   do c = 0,5
     slope := regCoef_n(ispan(branch_off(c),branch_off(c)+164,1),ann_ctrl({year|branch_off(c):(branch_off(c)+164)},latitude|:,longitude|:),0,0)
     copy_VarCoords(ann_ctrl(0,:,:),slope)
     intercept := slope@yintercept
     drift_trend(c,:,:,:) = conform(drift_trend(c,:,:,:),slope,(/1,2/))*conform(drift_trend(c,:,:,:),ispan(branch_off(c),branch_off(c)+164,1),0)+conform(drift_trend(c,:,:,:),onedtond(slope@yintercept,(/dimsizes(ann_ctrl&latitude),dimsizes(ann_ctrl&longitude)/)),(/1,2/))
   end do
  
; ----- new ukesm aka ukesm1.1 -------- 

  files := systemfunc("ls -1v ukesm1p1/historical/tas*nc")
  print(files)
  in := addfiles(files,"r")
  ListSetType(in,"join")

  ens := in[:]->UM_m01s03i236_vn1107(:,:,:,:)
  units := ens@units

  lonm := ens&longitude
  if (lonm(0) .ge. 0) then
    ens := lonFlip(ens)
  end if

  ann1p1 := month_to_annual(ens(time|:,ncl_join|:,latitude|:,longitude|:),1)
  DJF1p1 := month_to_season(ens(time|:,ncl_join|:,latitude|:,longitude|:),"DJF")
  JJA1p1 := month_to_season(ens(time|:,ncl_join|:,latitude|:,longitude|:),"JJA")
  ann1p1&year = ispan(1850,2014,1)
  DJF1p1&time = ispan(1850,2014,1)
  JJA1p1&time = ispan(1850,2014,1)

  copy_VarCoords(ann1p1(year|:,ncl_join|:,latitude|:,longitude|:),drift_trend(year|:,ncl_join|:,latitude|:,longitude|:))
  ddrifted := ann1p1(year|:,ncl_join|:,latitude|:,longitude|:)-drift_trend(year|:,ncl_join|:,latitude|:,longitude|:)
  ddrifted!0 = "year"
  ddrifted!1 = "ncl_join"
  ddrifted!2 = "latitude"
  ddrifted!3 = "longitude"
  copy_VarCoords(ann1p1(year|:,ncl_join|:,latitude|:,longitude|:),ddrifted(year|:,ncl_join|:,latitude|:,longitude|:))
  print(dim_avg_n_Wrap(ddrifted(ncl_join|:,year|:,latitude|:,longitude|:),(/1,2,3/)))

  ann1p1 := ddrifted
  
;  trends_1p1 = new((/6,144,192/),float)
;  do n = 0,5

  tval := new(dimsizes(ann1p1(0,0,:,:)),double)
  nptxy := new(dimsizes(ann1p1(0,0,:,:)),integer)
  
  trends_1p1 = regcoef(ispan(syear(i),eyear(i),1),dim_avg_n_Wrap(ann1p1(ncl_join|:,latitude|:,longitude|:,{year|syear(i):eyear(i)}),0),tval,nptxy)*10
  trends_1p1_all = regCoef(ispan(syear(i),eyear(i),1),ann1p1(ncl_join|:,latitude|:,longitude|:,{year|syear(i):eyear(i)}))*10
  printVarSummary(trends_1p1_all)
  copy_VarCoords(ann1p1(ncl_join|:,year|0,latitude|:,longitude|:),trends_1p1_all)
;  end do
;  trends_1p1 = regCoef_n(ispan(syear(i),eyear(i),1),dim_avg_n_Wrap((ann1p1({year|syear(i):eyear(i)},ncl_join|:,latitude|:,longitude|:)-drift_trend({year|syear(i):eyear(i)},ncl_join|:,latitude|:,longitude|:)),1),0,0)*10
  copy_VarCoords(ann1p1(0,0,:,:),trends_1p1)
  
; =================== read in observations ==================

  inobs = addfile("/gws/nopw/j04/smurphs/E/adittus/analysis/obs/HadCRUT5/HadCRUT.5.0.1.0.analysis.anomalies.ensemble_mean.nc", "r") ;obs/HadCR
  obs = inobs->tas_mean(0:1979,:,:)
  printVarSummary(obs)
  obs!0 = "time"
  lono = obs&longitude
  if (lono(0) .ge. 0) then
    obs = lonFlip(obs)
  end if
  obs_ann = month_to_annual(obs,1)
  obs_DJF = month_to_season(obs,"DJF")
  obs_JJA = month_to_season(obs,"JJA")
  obs_ann&year = ispan(1850,2014,1)  
  obs_DJF&time = ispan(1850,2014,1)
  obs_JJA&time = ispan(1850,2014,1)

;  ============ interpolate & smooth model data to obs resolution =====================

  ; DJF_LR = linint2_Wrap(DJF&longitude,DJF&latitude,DJF,True,obs_ann&longitude,obs_ann&latitude,0)
  ; JJA_LR = linint2_Wrap(JJA&longitude,JJA&latitude,JJA,True,obs_ann&longitude,obs_ann&latitude,0)

  ; DJF_bias = dim_avg_n_Wrap(DJF_LR({time|syear(i):eyear(i)},ncl_join|:,latitude|:,longitude|:),(/0,1/)) - dim_avg_n_Wrap(obs_DJF({time|syear(i):eyear(i)},latitude|:,longitude|:),0)
  ; JJA_bias = dim_avg_n_Wrap(JJA_LR({time|syear(i):eyear(i)},ncl_join|:,latitude|:,longitude|:),(/0,1/)) - dim_avg_n_Wrap(obs_JJA({time|syear(i):eyear(i)},latitude|:,longitude|:),0)
  ; copy_VarCoords(DJF_LR(0,0,:,:),JJA_bias)
  ; copy_VarCoords(JJA_LR(0,0,:,:),DJF_bias)
  otval := new(dimsizes(obs_ann(0,:,:)),double)
  onptxy := new(dimsizes(obs_ann(0,:,:)),integer)
  otrends = regcoef(ispan(syear(i),eyear(i),1),obs_ann(latitude|:,longitude|:,{year|syear(i):eyear(i)}),otval,onptxy)*10
  copy_VarCoords(obs_ann(0,:,:),otrends)
;  onptxy = where(onptxy .ge. 0.9*dimsizes(obs_ann(latitude|0,longitude|0,{year|syear(i):eyear(i)})),onptxy,onptxy@_FillValue)
  onptxy = where(onptxy .ge. 3,onptxy,onptxy@_FillValue)
  printVarSummary(otrends)
  printVarSummary(onptxy)
  otrends = where(onptxy .eq. onptxy@_FillValue,otrends@_FillValue,otrends)
  copy_VarCoords(obs_ann(0,:,:),otrends)
  rad     = 4.*atan(1.)/180.    ;
  lat = obs_ann&latitude
  clat = cos(lat*rad)

;  rmse(0+2*i) =  wgt_arearmse(dim_avg_n_Wrap(DJF_LR({time|syear(i):eyear(i)},ncl_join|:,latitude|:,longitude|:),(/0,1/)),dim_avg_n_Wrap(obs_DJF({time|syear(i):eyear(i)},latitude|:,longitude|:),0),clat,1.0,0)
;;  rmse(1) =  wgt_arearmse(JJA_oLR({time|syear(i):eyear(i)},latitude|:,longitude|:),JJA({time|syear(i):eyear(i)},ncl_join|0,latitude|:,longitude|:),clat,1.0,0)
;  rmse(1+2*i) =  wgt_arearmse(dim_avg_n_Wrap(JJA_LR({time|syear(i):eyear(i)},ncl_join|:,latitude|:,longitude|:),0),dim_avg_n_Wrap(obs_JJA({time|syear(i):eyear(i)},latitude|:,longitude|:),(/0,1/)),clat,1.0,0)

  diff = trends_1p1 - trends_1p0
  printMinMax(trends_1p1,True)

  probadiff = ttest(dim_avg_n_Wrap(trends_1p0_all,0),dim_stddev_n_Wrap(trends_1p0_all,0),6,dim_avg_n_Wrap(trends_1p1_all,0),dim_stddev_n_Wrap(trends_1p1_all,0),6,False,False)
  printVarSummary(probadiff)
  copy_VarCoords(trends_1p1,diff)
  copy_VarCoords(trends_1p1,probadiff)

  maps_all_LR = linint2_Wrap(trends_1p1_all&longitude,trends_1p1_all&latitude,trends_1p1_all,True,otrends&longitude,otrends&latitude,0)
  maps_all_LR1p0 = linint2_Wrap(trends_1p0_all&longitude,trends_1p0_all&latitude,trends_1p0_all,True,otrends&longitude,otrends&latitude,0)
  consistent_w_obs1p1 = new((/dimsizes(otrends&latitude),dimsizes(otrends&longitude)/),float) ; at the moment missing values in obs are categorised as inconsistent with obs. need to add 3rd category                                                                             
  consistent_w_obs1p0 = new((/dimsizes(otrends&latitude),dimsizes(otrends&longitude)/),float)

  do lati = 0,dimsizes(otrends&latitude)-1
    do lonj = 0,dimsizes(otrends&longitude)-1
      if (.not. ismissing(otrends(lati,lonj)) .and. (min(maps_all_LR(:,lati,lonj)) .le. otrends(lati,lonj)) .and. (max(maps_all_LR(:,lati,lonj)) .ge. otrends(lati,lonj))) then
        consistent_w_obs1p1(lati,lonj) = 1
      else
        consistent_w_obs1p1(lati,lonj) = 0
      end if
    end do
  end do

  do lati = 0,dimsizes(otrends&latitude)-1
    do lonj = 0,dimsizes(otrends&longitude)-1
      if (.not. ismissing(otrends(lati,lonj)) .and. (min(maps_all_LR1p0(:,lati,lonj)) .le. otrends(lati,lonj)) .and. (max(maps_all_LR1p0(:,lati,lonj)) .ge. otrends(lati,lonj))) then
        consistent_w_obs1p0(lati,lonj) = 1
      else
        consistent_w_obs1p0(lati,lonj) = 0
      end if
    end do
  end do

  copy_VarCoords(otrends,consistent_w_obs1p1)
  copy_VarCoords(otrends,consistent_w_obs1p0)

; ======================== plotting ==========================

  res = True
  res@gsnDraw = False
  res@gsnFrame = False
  res@mpPerimOn = False
  res@vpKeepAspect = True
  res@vpWidthF= 0.8
  res@vpHeightF= 0.6
  res@cnFillOn              = True     ; turn on color fill                                                                                                                                                                                                              
  res@cnFillMode            = "AreaFill"
  res@cnLinesOn             = False    ; turn off contour lines                                                                                                                                                                                                                          
  res@cnLineLabelsOn        = False

  res@pmTickMarkDisplayMode = "Always" 
  res@mpMinLatF             = -90
  res@mpLandFillColor           = "gray"
  res@gsnAddCyclic          = True
  res@lbLabelBarOn = False
  res@cnFillPalette = colormap(::-1,:)
  res@cnLevelSelectionMode = "ManualLevels"
  res@tmXBOn = False
  res@tmYLOn = False
  res@tmXTOn = False
  res@tmYROn = False
;  res@mpMinLatF = lat1
;  res@mpMaxLatF = lat2
  res@mpCenterLonF = 0.

  res@cnLevelSelectionMode = "ManualLevels"
  res@cnMinLevelValF = -0.6
  res@cnMaxLevelValF = 0.6
  res@cnLevelSpacingF = 0.05
  res@tiMainJust = "TopCenter"
  res@tiMainFontHeightF = 0.035
  res@tiMainFont = 21
  if (i .eq. 0) then
    res@tiMainString = "UKESM1.0"
  else
    res@tiMainString = ""
    delete(res@tiMainString)
  end if
  res@tiYAxisString = syear(i)+"-"+eyear(i)
  res@tiYAxisFontHeightF = 0.035
  printVarSummary(diff)
  plot(0+3*i) = gsn_csm_contour_map(wks,trends_1p0,res)
  delete(res@tiYAxisString)
  if (i .eq. 0) then
    res@tiMainString = "UKESM1.1"
  else
    res@tiMainString = ""
    delete(res@tiMainString)
  end if
  printVarSummary(trends_1p1)
  plot(1+3*i) = gsn_csm_contour_map(wks,trends_1p1,res)
  if (i .eq. 0) then
    res@tiMainString = "HadCRUT5"
  else
    res@tiMainString = ""
    delete(res@tiMainString)
  end if
  printVarSummary(otrends)
  plot(2+3*i) = gsn_csm_contour_map(wks,otrends,res)

  plots = new(6,graphic)
  do n = 0,5
    res@tiMainString = ""
;    plots(n) = gsn_csm_contour_map(wks,trends_1p1_all(n,:,:),res)
    plots(n) = gsn_csm_contour_map(wks,maps_all_LR(n,:,:),res)
  end do

  res2 = True
  res2@gsnDraw             = False
  res2@gsnFrame            = False
  res2@cnLevelSelectionMode = "ManualLevels"
  res2@cnMinLevelValF = 0
  res2@cnMaxLevelValF = 3.
  res2@cnLevelSpacingF = 0.05
  res2@cnInfoLabelOn       = False
  res2@cnLinesOn           = False
  res2@cnLineLabelsOn      = False
  res2@cnFillOn            = True
  res2@cnFillMode = "AreaFill"
  res2@cnFillScaleF        = 0.8
  res2@cnMonoFillColor = True
  res2@cnFillColor = "black"
;  res2@mpCenterLonF        = 180
  res2@tmXBOn = False
  res2@tmYLOn = False
  res2@tmXTOn = False
  res2@tmYROn = False
  res2@lbLabelBarOn = False
;  res2@gsnAddCyclic = True                                                                                       \
  
  res2@cnFillMode = "AreaFill"

  plot2   := gsn_csm_contour(wks,gsn_add_cyclic_point(consistent_w_obs1p0(:,:)), res2)
  opt     = True
  opt@gsnShadeFillType = "pattern"
  opt@gsnShadeLow = 17
;      opt@gsnShadeHigh = 6                                                                                        
  plot2   = gsn_contour_shade(plot2, 0.05, 0.05, opt)
  overlay(plot(0+3*i),plot2)
;  draw(plot2)
;  frame(wks)
; ====== observations =====
  
  plot2   := gsn_csm_contour(wks,gsn_add_cyclic_point(consistent_w_obs1p1(:,:)), res2)
  opt     = True
  opt@gsnShadeFillType = "pattern"
  opt@gsnShadeLow = 17
;      opt@gsnShadeHigh = 6                                                                                     \
                                                                                                                 
  plot2   = gsn_contour_shade(plot2, 0.05, 0.05, opt)
  overlay(plot(1+3*i),plot2)

;                                 ; ========= probadiff =====


;   plot2   := gsn_csm_contour(wks,gsn_add_cyclic_point(probadiff(:,:)), res2)
;   printMinMax(probadiff,True)
;   opt     = True
;   opt@gsnShadeFillType = "pattern"
;   opt@gsnShadeHigh = 17
; ;      opt@gsnShadeHigh = 6

;   plot2   = gsn_contour_shade(plot2, 0.05, 0.05, opt)
; ;  overlay(plot(1+3*i),plot2)

  panelres1 = True
  panelres1@txString = "UKESM1.1 ensemble members "+syear(i)+" "+eyear(i)
  panelres1@gsnPanelLabelBar = True
  gsn_panel(wks,plots,(/2,3/),panelres1)
  
end do

                                ; getvalues plot(2)@contour    ; retrieve colors/levels from plot
                                ;   "cnFillColors" : colors
                                ;   "cnLevels"     : levels
                                ; end getvalues

                                ; ncolors = dimsizes(colors)
                                ; colors(ncolors/2-1) = 0    ; set two middle colors to white
                                ; colors(ncolors/2) = 0       ; if you have odd number of colors this

                                ; print(colors)

                                ; setvalues plot(2)@contour
                                ;   "cnFillColors" : colors       ; reset the colors
                                ; end setvalues


pres = True
gsnPanelDebug=True
pres@gsnPanelLeft = 0.05
pres@gsnPanelScalePlotIndex = 4
pres@gsnPanelLabelBar = True
pres@lbTitleString = "K/decade"
pres@lbBoxEndCapStyle           = "TriangleBothEnds"
;pres@lbTitleJust = "TopCenter"
pres@lbTitleOffsetF = 0.12
;pres@gsnPanelFigureStrings    := conform(2,"RMSE = ",-1)+decimalPlaces(rmse(:),2,True)
pres@gsnPanelYWhiteSpacePercent = 2.
pres@gsnPanelXWhiteSpacePercent = 1.
pres@lbTitleFontHeightF = 0.015
pres@pmLabelBarOrthogonalPosF = -0.025
gsn_panel(wks,plot(:),(/3,3/),pres)
end
