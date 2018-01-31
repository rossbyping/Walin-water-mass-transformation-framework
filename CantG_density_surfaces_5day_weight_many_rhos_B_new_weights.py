import numpy as np
import numpy.ma as ma
import matplotlib.pyplot as plt
import scipy.io as sio
import scipy.interpolate as ip
import os
import time
from datetime import datetime
from mpl_toolkits.basemap import Basemap, addcyclic
from scipy.ndimage.filters import minimum_filter, maximum_filter
from netCDF4 import Dataset
ticc=time.time()

#################costant parameters
Cp=3992.1
const1=1000.0
rho0=1025.0
deltarho=0.1
deltn='01'

#################
filename_area_t='/work/Ping.Zhai/MOM5_CORE_B_new/data/19800101.ocean_grid.nc'
f_area_t=Dataset(filename_area_t,'r')

lons1=f_area_t.variables['xt_ocean'][:]
lons1[(lons1<0)]=lons1[(lons1<0)]+360
lats1=f_area_t.variables['yt_ocean'][:]

indexlon=np.where((lons1>100) & (lons1<290))
indexlat=np.where((lats1>-20) & (lats1<20))
lons=lons1[indexlon[0]]
lats=lats1[indexlat[0]]
rhos=np.arange(20.,28.7,deltarho)

nrho=rhos.size
nlon=lons.size
nlat=lats.size
nz=50

area_t=f_area_t.variables['area_t'][indexlat[0],indexlon[0]]
area=np.zeros((nz,nlat,nlon))
for k in range(0,nz):
    area[k,:,:]=area_t
f_area_t.close()

filename_ocean='/work/Ping.Zhai/MOM5_CORE_B_new/data/19800101.ocean_new.nc'
f_ocean=Dataset(filename_ocean,'r')
z=f_ocean.variables['st_ocean'][:]
zw=f_ocean.variables['sw_ocean'][:]
zw=np.insert(zw,0,0)
f_ocean.close()

#####prepare the sig0ab
sig0ab=np.zeros((nrho*2))
for n in range(0,nrho):
    sig0ab[n*2]=rhos[n]-deltarho/2
    sig0ab[n*2+1]=rhos[n]+deltarho/2

yearexp=1980
for yearn in range(yearexp,yearexp+1):

    filename_budgets='/work/Ping.Zhai/MOM5_CORE_B_new/data/'+str(yearn)+'0101.ocean_budgets_new.nc'
    filename_bdy='/work/Ping.Zhai/MOM5_CORE_B_new/data/'+str(yearn)+'0101.ocean_bdy_flux.nc'
    filename_dsigma0='/work/Ping.Zhai/MOM5_CORE_B_new/data/'+str(yearn)+'0101.ocean_dsigma0dTS.nc'
    filename_sigma0='/work/Ping.Zhai/MOM5_CORE_B_new/data/'+str(yearn)+'0101.ocean_sigma0.nc'
    filename_ocean='/work/Ping.Zhai/MOM5_CORE_B_new/data/'+str(yearn)+'0101.ocean_new.nc'
    filename_dic_B='/work/Ping.Zhai/MOM5_CORE_B_new/data/'+str(yearn)+'0101.ocean_topaz_term_dic.nc'
    filename_dic_A='/work/Ping.Zhai/MOM5_CORE_A_new/data/'+str(yearn)+'0101.ocean_topaz_term_dic.nc'

    f_budgets =Dataset(filename_budgets,'r')
    f_bdy =Dataset(filename_bdy,'r')
    f_dsigma0 =Dataset(filename_dsigma0,'r')
    f_sigma0 =Dataset(filename_sigma0,'r')
    f_ocean=Dataset(filename_ocean,'r')
    f_dic_B=Dataset(filename_dic_B,'r')
    f_dic_A=Dataset(filename_dic_A,'r')

    datetime=f_budgets.variables['time'][:]
    nt=datetime.size

######write the result into netcdf file
    CantG_isop = Dataset('/work/Ping.Zhai/EqPacific_work/B_new/CantG_density_surface_'+str(yearn)+'_5day_weights_many_rhos_B_new_weights_add_new_interps_'+deltn+'.nc', 'w', format='NETCDF4')
    CantG_isop.description = 'Cant diapycnal transport due to water mass transformation'
    CantG_isop.createDimension('time', nt)
    CantG_isop.createDimension('rhos', nrho)
    CantG_isop.createDimension('lat', nlat)
    CantG_isop.createDimension('lon', nlon)
    timefile = CantG_isop.createVariable('time', 'f4', ('time',))
    rhosfile = CantG_isop.createVariable('rhos', 'f4', ('rhos',))
    latfile = CantG_isop.createVariable('lat', 'f4', ('lat',))
    lonfile = CantG_isop.createVariable('lon', 'f4', ('lon',))
    timefile.units="days since 1860-01-01 00:00:00, no leap year"
    timefile[:]=datetime
    rhosfile[:]=rhos
    latfile[:]=lats
    lonfile[:]=lons

    CantG_S_tendencyfile=CantG_isop.createVariable('CantG_S_tendency', 'f4', ('time','rhos', 'lat', 'lon',),fill_value=-1.0e30)
    CantG_S_advectionfile=CantG_isop.createVariable('CantG_S_advection', 'f4', ('time','rhos', 'lat', 'lon',),fill_value=-1.0e30)
    CantG_S_submesofile=CantG_isop.createVariable('CantG_S_submeso', 'f4', ('time','rhos', 'lat', 'lon',),fill_value=-1.0e30)
    CantG_S_gmfile=CantG_isop.createVariable('CantG_S_gm', 'f4', ('time','rhos', 'lat', 'lon',),fill_value=-1.0e30)
    CantG_S_vdifffile=CantG_isop.createVariable('CantG_S_vdiff', 'f4', ('time','rhos', 'lat', 'lon',),fill_value=-1.0e30)
    CantG_S_vdiff_K33file=CantG_isop.createVariable('CantG_S_vdiff_K33', 'f4', ('time','rhos', 'lat', 'lon',),fill_value=-1.0e30)
    CantG_S_vdiff_cbtfile=CantG_isop.createVariable('CantG_S_vdiff_cbt', 'f4', ('time','rhos', 'lat', 'lon',),fill_value=-1.0e30) 
    CantG_S_hdifffile=CantG_isop.createVariable('CantG_S_hdiff', 'f4', ('time','rhos', 'lat', 'lon',),fill_value=-1.0e30)
    CantG_S_KPPfile=CantG_isop.createVariable('CantG_S_KPP', 'f4', ('time','rhos', 'lat', 'lon',),fill_value=-1.0e30)
    CantG_S_surf_pme_file=CantG_isop.createVariable('CantG_S_surf_pme', 'f4', ('time','rhos', 'lat', 'lon',),fill_value=-1.0e30)
    CantG_S_resifile=CantG_isop.createVariable('CantG_S_resi', 'f4', ('time','rhos', 'lat', 'lon',),fill_value=-1.0e30)
    CantG_S_surf_flux_file=CantG_isop.createVariable('CantG_S_surf_flux', 'f4', ('time','rhos', 'lat', 'lon',),fill_value=-1.0e30)

    CantG_T_tendencyfile=CantG_isop.createVariable('CantG_T_tendency', 'f4', ('time','rhos', 'lat', 'lon',),fill_value=-1.0e30)
    CantG_T_advectionfile=CantG_isop.createVariable('CantG_T_advection', 'f4', ('time','rhos', 'lat', 'lon',),fill_value=-1.0e30)
    CantG_T_submesofile=CantG_isop.createVariable('CantG_T_submeso', 'f4', ('time','rhos', 'lat', 'lon',),fill_value=-1.0e30)
    CantG_T_gmfile=CantG_isop.createVariable('CantG_T_gm', 'f4', ('time','rhos', 'lat', 'lon',),fill_value=-1.0e30)
    CantG_T_vdifffile=CantG_isop.createVariable('CantG_T_vdiff', 'f4', ('time','rhos', 'lat', 'lon',),fill_value=-1.0e30)
    CantG_T_vdiff_K33file=CantG_isop.createVariable('CantG_T_vdiff_K33', 'f4', ('time','rhos', 'lat', 'lon',),fill_value=-1.0e30)
    CantG_T_vdiff_cbtfile=CantG_isop.createVariable('CantG_T_vdiff_cbt', 'f4', ('time','rhos', 'lat', 'lon',),fill_value=-1.0e30)
    CantG_T_hdifffile=CantG_isop.createVariable('CantG_T_hdiff', 'f4', ('time','rhos', 'lat', 'lon',),fill_value=-1.0e30)
    CantG_T_KPPfile=CantG_isop.createVariable('CantG_T_KPP', 'f4', ('time','rhos', 'lat', 'lon',),fill_value=-1.0e30)
    CantG_T_surffile=CantG_isop.createVariable('CantG_T_surf', 'f4', ('time','rhos', 'lat', 'lon',),fill_value=-1.0e30)
    CantG_T_surf_sensible_file=CantG_isop.createVariable('CantG_T_surf_sensible', 'f4', ('time','rhos', 'lat', 'lon',),fill_value=-1.0e30)
    CantG_T_surf_latent_file=CantG_isop.createVariable('CantG_T_surf_latent', 'f4', ('time','rhos', 'lat', 'lon',),fill_value=-1.0e30)
    CantG_T_surf_shortwave_file=CantG_isop.createVariable('CantG_T_surf_shortwave', 'f4', ('time','rhos', 'lat', 'lon',),fill_value=-1.0e30)
    CantG_T_swfile=CantG_isop.createVariable('CantG_T_sw', 'f4', ('time','rhos', 'lat', 'lon',),fill_value=-1.0e30)
    CantG_T_resifile=CantG_isop.createVariable('CantG_T_resi', 'f4', ('time','rhos', 'lat', 'lon',),fill_value=-1.0e30)


########load the salinity data
    for t in range(0,nt):
#        dzt=f_budgets.variables['dzt'][t,:,indexlat[0],indexlon[0]]
        drhodS=f_dsigma0.variables['dsigma0dsalinity'][t,:,indexlat[0],indexlon[0]]/1.0
        drhodT=f_dsigma0.variables['dsigma0dtheta'][t,:,indexlat[0],indexlon[0]]
        Cant=f_dic_A.variables['dic'][t,:,indexlat[0],indexlon[0]]/1.0-f_dic_B.variables['dic'][t,:,indexlat[0],indexlon[0]]/1.0

        salt_tendency=f_budgets.variables['salt_tendency'][t,:,indexlat[0],indexlon[0]]*const1*area*drhodS*Cant        
        salt_vdiffuse_impl=f_budgets.variables['salt_vdiffuse_impl'][t,:,indexlat[0],indexlon[0]]*const1*area*drhodS*Cant
        salt_vdiffuse_K33=f_budgets.variables['salt_vdiffuse_k33'][t,:,indexlat[0],indexlon[0]]*const1*area*drhodS*Cant
        salt_vdiffuse_cbt=f_budgets.variables['salt_vdiffuse_diff_cbt'][t,:,indexlat[0],indexlon[0]]*const1*area*drhodS*Cant
        salt_advection=f_budgets.variables['salt_advection'][t,:,indexlat[0],indexlon[0]]*const1*area*drhodS*Cant
        salt_nonlocal_KPP=f_budgets.variables['salt_nonlocal_KPP'][t,:,indexlat[0],indexlon[0]]*const1*area*drhodS*Cant
        salt_submeso=f_budgets.variables['salt_submeso'][t,:,indexlat[0],indexlon[0]]*const1*area*drhodS*Cant
        neutral_diffusion_salt=f_budgets.variables['neutral_diffusion_salt'][t,:,indexlat[0],indexlon[0]]*const1*area*drhodS*Cant
        neutral_gm_salt=f_budgets.variables['neutral_gm_salt'][t,:,indexlat[0],indexlon[0]]*const1*area*drhodS*Cant
        S_resi=f_budgets.variables['salt_rivermix'][t,:,indexlat[0],indexlon[0]]*const1*area*drhodS*Cant
#        S_resi=S_resi+f_budgets.variables['salt_runoffmix'][t,:,indexlat[0],indexlon[0]]*const1*area*drhodS*Cant
#        S_resi=S_resi+f_budgets.variables['salt_calvingmix'][t,:,indexlat[0],indexlon[0]]*const1*area*drhodS*Cant
#        S_resi=S_resi+f_budgets.variables['salt_xland'][t,:,indexlat[0],indexlon[0]]*const1*area*drhodS*Cant
#        S_resi=S_resi+f_budgets.variables['salt_xlandinsert'][t,:,indexlat[0],indexlon[0]]*const1*area*drhodS*Cant
#        S_resi=S_resi+f_budgets.variables['salt_sigma_diff'][t,:,indexlat[0],indexlon[0]]*const1*area*drhodS*Cant
#        S_resi=S_resi+f_budgets.variables['mixdownslope_salt'][t,:,indexlat[0],indexlon[0]]*const1*area*drhodS*Cant
#        salt_eta_smooth=f_budgets.variables['salt_eta_smooth'][t,indexlat[0],indexlon[0]]*const1*area[0,:,:]*drhodS[0,:,:]*Cant[0,:,:]
#        salt_eta_smooth=[]
        salt_surface_flux=np.nansum(salt_vdiffuse_impl,axis=0)
        salt=f_ocean.variables['salt'][t,0,indexlat[0],indexlon[0]]/1.0
        pme=f_ocean.variables['pme_river'][t,indexlat[0],indexlon[0]]/1.0-f_ocean.variables['river'][t,indexlat[0],indexlon[0]]/1.0
        salt=np.float64(salt)
        pme=np.float64(pme) 
        pme_salt=pme*salt*area[0,:,:]*drhodS[0,:,:]*Cant[0,:,:]       
        S_resi[0,:,:]=S_resi[0,:,:]+pme_salt
        salt_vdiffuse_impl[0,:,:]=salt_vdiffuse_impl[0,:,:]-salt_surface_flux


    ########load the temperature data
        temp_tendency=f_budgets.variables['temp_tendency'][t,:,indexlat[0],indexlon[0]]/Cp*area*drhodT*Cant
        temp_vdiffuse_impl=f_budgets.variables['temp_vdiffuse_impl'][t,:,indexlat[0],indexlon[0]]/Cp*area*drhodT*Cant
        temp_vdiffuse_K33=f_budgets.variables['temp_vdiffuse_k33'][t,:,indexlat[0],indexlon[0]]/Cp*area*drhodT*Cant
        temp_vdiffuse_cbt=f_budgets.variables['temp_vdiffuse_diff_cbt'][t,:,indexlat[0],indexlon[0]]/Cp*area*drhodT*Cant
        temp_advection=f_budgets.variables['temp_advection'][t,:,indexlat[0],indexlon[0]]/Cp*area*drhodT*Cant
        temp_nonlocal_KPP=f_budgets.variables['temp_nonlocal_KPP'][t,:,indexlat[0],indexlon[0]]/Cp*area*drhodT*Cant
        temp_submeso=f_budgets.variables['temp_submeso'][t,:,indexlat[0],indexlon[0]]/Cp*area*drhodT*Cant
        neutral_diffusion_temp=f_budgets.variables['neutral_diffusion_temp'][t,:,indexlat[0],indexlon[0]]/Cp*area*drhodT*Cant
        neutral_gm_temp=f_budgets.variables['neutral_gm_temp'][t,:,indexlat[0],indexlon[0]]/Cp*area*drhodT*Cant
        sw_heat=f_budgets.variables['sw_heat'][t,:,indexlat[0],indexlon[0]]/Cp*area*drhodT*Cant
    
        T_resi=f_budgets.variables['temp_rivermix'][t,:,indexlat[0],indexlon[0]]/Cp*area*drhodT*Cant
#        T_resi=T_resi+f_budgets.variables['temp_runoffmix'][t,:,indexlat[0],indexlon[0]]/Cp*area*drhodT*Cant
#        T_resi=T_resi+f_budgets.variables['temp_calvingmix'][t,:,indexlat[0],indexlon[0]]/Cp*area*drhodT*Cant
#         T_resi=T_resi+f_budgets.variables['temp_xland'][t,:,indexlat[0],indexlon[0]]/Cp*area*drhodT*Cant
#        T_resi=T_resi+f_budgets.variables['temp_xlandinsert'][t,:,indexlat[0],indexlon[0]]/Cp*area*drhodT*Cant
#        T_resi=T_resi+f_budgets.variables['temp_sigma_diff'][t,:,indexlat[0],indexlon[0]]/Cp*area*drhodT*Cant
#        T_resi=T_resi+f_budgets.variables['mixdownslope_temp'][t,:,indexlat[0],indexlon[0]]/Cp*area*drhodT*Cant
#        temp_eta_smooth=f_budgets.variables['temp_eta_smooth'][t,indexlat[0],indexlon[0]]/Cp*area[0,:,:]*drhodT[0,:,:]*Cant[0,:,:]
#         frazil_2d=f_budgets.variables['frazil_2d'][t,indexlat[0],indexlon[0]]/Cp*area[0,:,:]*drhodT[0,:,:]*Cant[0,:,:]
        sfc_hflux_pme=f_budgets.variables['sfc_hflux_pme'][t,indexlat[0],indexlon[0]]/Cp*area[0,:,:]*drhodT[0,:,:]*Cant[0,:,:]
        T_resi[0,:,:]=T_resi[0,:,:]+sfc_hflux_pme
        temp_eta_smooth=[]
        frazil_2d=[]
        sfc_hflux_pme=[]
    
        sfc_heat_coupler=f_budgets.variables['sfc_hflux_coupler'][t,indexlat[0],indexlon[0]]/Cp*area[0,:,:]*drhodT[0,:,:]*Cant[0,:,:]
        sensible=f_bdy.variables['sens_heat'][t,indexlat[0],indexlon[0]]/Cp*area[0,:,:]*drhodT[0,:,:]*Cant[0,:,:]
        latent=f_bdy.variables['evap_heat'][t,indexlat[0],indexlon[0]]/Cp*area[0,:,:]*drhodT[0,:,:]*Cant[0,:,:]
        shortwave=f_bdy.variables['swflx'][t,indexlat[0],indexlon[0]]/Cp*area[0,:,:]*drhodT[0,:,:]*Cant[0,:,:]

        temp_vdiffuse_impl[0,:,:]=temp_vdiffuse_impl[0,:,:]-sfc_heat_coupler 
        sig0=np.around(f_sigma0.variables['sigma0'][t,:,indexlat[0],indexlon[0]], decimals=4)
        sig0[(ma.getmask(sig0))]=np.nan
       
        depthab=np.zeros((nrho*2,nlat,nlon))+np.nan
        for j in range(0,nlat):
            for i in range(0,nlon):

                indextmp=np.where((sig0[:,j,i]>3) & (sig0[:,j,i]<44) & np.logical_not(ma.getmask(sig0[:,j,i])))
                sig0tmp=sig0[indextmp[0],j,i]
                ztmp=z[indextmp[0]]

                indextmp=np.argsort(sig0tmp)
                sig0tmp=sig0tmp[indextmp]
                ztmp=ztmp[indextmp]

                difftmp=np.diff(sig0tmp)
                indextmpi=np.where((np.abs(difftmp)<8.0e-5))
                ztmp=np.delete(ztmp,indextmpi[0])
                sig0tmp=np.delete(sig0tmp,indextmpi[0])

                if sig0tmp.size>2:
                   depthab[:,j,i]=ip.interp1d(sig0tmp,ztmp,kind='linear',bounds_error=False,assume_sorted=True)(sig0ab)

        #################loop density

        for n in range(0,nrho):
            zweights=np.zeros((nz,nlat,nlon))

            for j in range(0,nlat):
                for i in range(0,nlon):

                    indexz=np.where((sig0[:,j,i]>sig0ab[n*2]) & (sig0[:,j,i]<sig0ab[n*2+1]) & np.logical_not(ma.getmask(sig0[:,j,i])))
                    indexz=indexz[0]

                    if indexz.size>0:

                       if indexz[-1]+1<nz:
                          if sig0[indexz[-1]+1,j,i]<sig0ab[n*2+1]:
                             indexz=np.append(indexz,indexz[-1]+1)

                       if indexz[0]-1>-1:
                          if sig0[indexz[0]-1,j,i]>sig0ab[n*2]:
                             indexz=np.append(indexz[0]-1,indexz)
                       zweights[indexz,j,i]=1

                       if indexz[0]>0 and sig0ab[n*2]>0 :

                          if depthab[n*2,j,i]>=zw[indexz[0]]:
                             zweights[indexz[0],j,i]=(zw[indexz[0]+1]-depthab[n*2,j,i])/(zw[indexz[0]+1]-zw[indexz[0]])
                          if depthab[n*2,j,i]<zw[indexz[0]]:
                             zweights[indexz[0]-1,j,i]=(zw[indexz[0]]-depthab[n*2,j,i])/(zw[indexz[0]]-zw[indexz[0]-1])

                       if indexz[-1]<(nz-1) and np.logical_not(ma.getmask(sig0[indexz[-1]+1,j,i])) and sig0ab[n*2+1]<30:

                          if depthab[n*2+1,j,i]<=zw[indexz[-1]+1]:
                             zweights[indexz[-1],j,i]=(depthab[n*2+1,j,i]-zw[indexz[-1]])/(zw[indexz[-1]+1]-zw[indexz[-1]])
                          if depthab[n*2+1,j,i]>zw[indexz[-1]+1]:
                             zweights[indexz[-1]+1,j,i]=(depthab[n*2+1,j,i]-zw[indexz[-1]+1])/(zw[indexz[-1]+2]-zw[indexz[-1]+1])

                    else:
                        indexupper=np.where((sig0[:,j,i]<sig0ab[n*2]))
                        indexlower=np.where((sig0[:,j,i]>sig0ab[n*2+1]))
                        if indexupper[0].size>0 and indexlower[0].size>0 and indexlower[0][0]>indexupper[0][-1] and np.nansum(zweights[:,j,i])<0.1:

                           if depthab[n*2,j,i]>zw[indexupper[0][-1]+1]:
                              zweights[indexlower[0][0],j,i]=(depthab[n*2+1,j,i]-depthab[n*2,j,i])/(zw[indexlower[0][0]+1]-zw[indexlower[0][0]])
                           if depthab[n*2,j,i]<zw[indexupper[0][-1]+1] and depthab[n*2+1,j,i]>zw[indexupper[0][-1]+1]:
                              zweights[indexupper[0][-1],j,i]=(zw[indexupper[0][-1]+1]-depthab[n*2,j,i])/(zw[indexupper[0][-1]+1]-zw[indexupper[0][-1]])
                              zweights[indexlower[0][0],j,i]=(depthab[n*2+1,j,i]-zw[indexlower[0][0]])/(zw[indexlower[0][0]+1]-zw[indexlower[0][0]])
                           if depthab[n*2+1,j,i]<zw[indexupper[0][-1]+1]:
                              zweights[indexupper[0][-1],j,i]=(depthab[n*2+1,j,i]-depthab[n*2,j,i])/(zw[indexupper[0][-1]+1]-zw[indexupper[0][-1]])

            CantG_S_tendencyfile[t,n,:,:]=np.nansum(salt_tendency*zweights,axis=0)/deltarho            
            CantG_S_advectionfile[t,n,:,:]=-np.nansum(salt_advection*zweights,axis=0)/deltarho
            CantG_S_submesofile[t,n,:,:]=-np.nansum(salt_submeso*zweights,axis=0)/deltarho
            CantG_S_gmfile[t,n,:,:]=-np.nansum(neutral_gm_salt*zweights,axis=0)/deltarho
            CantG_S_vdifffile[t,n,:,:]=np.nansum(salt_vdiffuse_impl*zweights,axis=0)/deltarho
            CantG_S_vdiff_K33file[t,n,:,:]=np.nansum(salt_vdiffuse_K33*zweights,axis=0)/deltarho
            CantG_S_vdiff_cbtfile[t,n,:,:]=np.nansum(salt_vdiffuse_cbt*zweights,axis=0)/deltarho
            CantG_S_hdifffile[t,n,:,:]=np.nansum(neutral_diffusion_salt*zweights,axis=0)/deltarho
            CantG_S_KPPfile[t,n,:,:]=np.nansum(salt_nonlocal_KPP*zweights,axis=0)/deltarho
            CantG_S_resifile[t,n,:,:]=-np.nansum(S_resi*zweights,axis=0)/deltarho

            CantG_T_tendencyfile[t,n,:,:]=np.nansum(temp_tendency*zweights,axis=0)/deltarho
            CantG_T_advectionfile[t,n,:,:]=-np.nansum(temp_advection*zweights,axis=0)/deltarho
            CantG_T_submesofile[t,n,:,:]=-np.nansum(temp_submeso*zweights,axis=0)/deltarho
            CantG_T_gmfile[t,n,:,:]=-np.nansum(neutral_gm_temp*zweights,axis=0)/deltarho
            CantG_T_vdifffile[t,n,:,:]=np.nansum(temp_vdiffuse_impl*zweights,axis=0)/deltarho
            CantG_T_vdiff_K33file[t,n,:,:]=np.nansum(temp_vdiffuse_K33*zweights,axis=0)/deltarho
            CantG_T_vdiff_cbtfile[t,n,:,:]=np.nansum(temp_vdiffuse_cbt*zweights,axis=0)/deltarho
            CantG_T_hdifffile[t,n,:,:]=np.nansum(neutral_diffusion_temp*zweights,axis=0)/deltarho
            CantG_T_KPPfile[t,n,:,:]=np.nansum(temp_nonlocal_KPP*zweights,axis=0)/deltarho
            CantG_T_swfile[t,n,:,:]=np.nansum(sw_heat*zweights,axis=0)/deltarho
            CantG_T_resifile[t,n,:,:]=-np.nansum(T_resi*zweights,axis=0)/deltarho
        
            CantG_S_surf_pme_file[t,n,:,:]=-pme_salt*zweights[0,:,:]/deltarho
            CantG_S_surf_flux_file[t,n,:,:]=salt_surface_flux*zweights[0,:,:]/deltarho
            CantG_T_surffile[t,n,:,:]=sfc_heat_coupler*zweights[0,:,:]/deltarho
            CantG_T_surf_sensible_file[t,n,:,:]=sensible*zweights[0,:,:]/deltarho
            CantG_T_surf_latent_file[t,n,:,:]=latent*zweights[0,:,:]/deltarho
            CantG_T_surf_shortwave_file[t,n,:,:]=shortwave*zweights[0,:,:]/deltarho

        salt_tendency=[]
        salt_advection=[]
        salt_submeso=[]
        neutral_gm_salt=[]
        salt_vdiffuse_impl=[]
        salt_vdiffuse_K33=[]
        salt_vdiffuse_cbt=[]
        neutral_diffusion_salt=[]
        salt_nonlocal_KPP=[]
        S_resi=[]
        temp_tendency=[]
        temp_advection=[]
        temp_submeso=[]
        neutral_gm_temp=[]
        temp_vdiffuse_impl=[]
        temp_vdiffuse_K33=[]
        temp_vdiffuse_cbt=[]
        neutral_diffusion_temp=[]
        temp_nonlocal_KPP=[]
        sw_heat=[]
        T_resi=[]


    f_budgets.close()
    f_dsigma0.close()
    f_sigma0.close()
    f_ocean.close() 
    f_dic_A.close()
    f_dic_B.close()
    CantG_isop.close()

#file.close()
elipsetime= time.time()-ticc
sio.savemat('elptimeCant.mat',{'elipsetime':elipsetime})

