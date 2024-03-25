# -*- coding: utf-8 -*-
"""
Created on Thu Oct 27 15:32:17 2022

@author: stefano.pauletta
"""
import xlwt
import xlrd
import pvlib
import pandas as pd
import pyarrow.feather as feather
import numpy as np
import geopandas
from shapely.geometry import Point
import os
import scipy
from scipy.cluster import vq
from xlutils.copy import copy
from optihood.calpinage_dev2 import Calpinage_light as cp
import math

class weather:
    """
    Class to download Tipical Meteorogical Year data from PVGIS for a given location.
    Data are saved in the target file while location is assigned from the source excel file ("scenario.xls" by default)
    Soil probe temperature is derived from surface soil temperature at the nearest MeteoSchweiz station, the 0°C limit and a sinusoidal profile along the year.
    A file target is written in a csv.
    """
    
    def __init__(self, 
                 source="scenario.xls",               
                 cluster=True,
                 n_clusters=20,
                 soil_param=True,
                 T_soil_min=0,
                 Day_soil_min=155,
                 T_soil_max=10,
                 clustering_vars=[],
                 save_file=False,
                 load_file=False,
                 ):
        """
        Class constructor. Geographic info are taken from the source file
        Methods are applied at the end in sequence to create weather file for optimease
        """
        self.get_scenario(source)
        self.tilt=0
        self.az=0
        self.coordinates = [self.lat,self.long,self.tilt,self.az]
        self.w=1
        self.l=2
        self.time_y="time.yy"
        self.time_m="time.mm"
        self.time_d="time.dd"
        self.time_h="time.hh"
        self.lb_air="tre200h0"          # tre200s0
        self.lb_ghi="gls"               # gre000z0  
        self.lb_dhi="str.diffus"        # ods000z0  
        self.lb_gT="ground_temp"
        self.lb_RH="Rel_Hum"            # ure200s0  
        self.lb_wnd="wind_speed"        # fkl010z0  
        self.lb_p="pressure"            # prestas0 
        self.lb_wnd_d="wind_direction"  # dkl010z0  
        self.lb_IR="IR"                 # oli000z0    
        self.n_cluster=n_clusters
        
        if soil_param==True:
            self.T_soil_min=T_soil_min
            self.Day_soil_min=Day_soil_min
            self.T_soil_max=T_soil_max
        else:
            self.T_soil_min=0
            self.Day_soil_min=58
            self.T_soil_max=999
        self.load_file=load_file
        self.get_TMY_formatted()
        if self.load_file==False:
            self.get_soil()
        self.set_solar_scenarios()
        if cluster==True:
            self.do_clustering(clustering_vars)
        if save_file!=False:
            self.save_DB(self.target)
        
        return None
    
    def get_scenario(self,source):
        """
        Method to access the source file and locate data concerning the solar simulations
        
        """
        workbook = xlrd.open_workbook(source)

        worksheet = workbook.sheet_by_name('profiles')
        self.solar_f_adr = worksheet.cell(3, worksheet.row_values(0, start_colx=0, end_colx=None).index('path')).value
        self.demand_f_adr = worksheet.cell(2, worksheet.row_values(0, start_colx=0, end_colx=None).index('path')).value
        if self.solar_f_adr!='' or self.solar_f_adr!='PVGIS':
            self.solar_file=True
        else:
            self.solar_file=False
        self.target=worksheet.cell(1, worksheet.row_values(0, start_colx=0, end_colx=None).index('path')).value
        
        
        """
        Load new sheet, compute centroid of building portfolio
        and use it for solar radiation data calculations
        """
        worksheet = workbook.sheet_by_name('hood')
        df_columns=worksheet.row_values(0,start_colx=0,end_colx=None)
        self.df_hood=pd.DataFrame(columns=df_columns)
        i=0
        for i in range(len(df_columns)):
            self.df_hood[df_columns[i]]=worksheet.col_values(i,start_rowx=1,end_rowx=None)
        
        geometry=geopandas.GeoSeries(geopandas.points_from_xy(self.df_hood.loc[:,'longitude'], 
                                                              self.df_hood.loc[:,'latitude']))
        self.lat=geometry.unary_union.centroid.y
        self.long=geometry.unary_union.centroid.x
        self.W_S_D=-23.45
        self.W_S_A=np.abs(self.lat+self.W_S_D)
        # self.set_solar_scenarios()
        
        # worksheet = workbook.sheet_by_name('solar')
        # self.tilt = float(worksheet.cell(1, worksheet.row_values(0, start_colx=0, end_colx=None).index('tilt')).value)
        # self.az = float(worksheet.cell(1, worksheet.row_values(0, start_colx=0, end_colx=None).index('azimuth')).value)
        return None
    
    def load_bld_demand(self,bld_name):
        source=self.demand_f_adr + r"\Building2"+str(bld_name)+".csv"
        demand=pd.read_csv(source,sep=";")
        
        return demand
        
    def minimum_distance_tilt(self):
        def func(tilt,w,dist,W_S_A):
            value=self.w*math.sin(math.radians(self.tilt))*1/np.tan(math.radians(self.W_S_A))
            return value
        return None
        
        
    
    def set_solar_scenarios(self):
            """
            Method to set the solar scenarios with optimized azimuth and tilt 
            """
            # df_columns=worksheet.row_values(0,start_colx=0,end_colx=None)
            # df=pd.DataFrame(columns=df_columns)
            # i=0
            # for i in range(len(df_columns)):
            #     df[df_columns[i]]=worksheet.col_values(i,start_rowx=1,end_rowx=None)
            
            # geometry=geopandas.GeoSeries(geopandas.points_from_xy(df.loc[:,'longitude'], 
            #                                                       df.loc[:,'latitude']))
            # self.lat=geometry.unary_union.centroid.y
            # self.long=geometry.unary_union.centroid.x
            
            '''
            we subtract the busy area as a section of the long side
            '''
            self.df_hood.loc[:,'long_side']=self.df_hood.loc[:,'long_side']-self.df_hood.loc[:,'busy_area']/self.df_hood.loc[:,'short_side']
            self.df_hood.loc[:,'roof_area']=self.df_hood.loc[:,'long_side']*self.df_hood.loc[:,'short_side']
            """
            once the building geometry is defined, the class Calpinage can 
            be called to compute the roof coverage ratio, the number of solar panels
            and the total receiving surface for each optimization scenario for 
            orientation: same of building on short side
                        same of buidlgin on long side
                        full south (180° according to used convention)
            and structure:  portait, with optimal tilt for each case as computed by PVlib
                            East West (EW) with fixed tilt of 20°
            """
            # df['P_Coverage_90']=0
            # df['P_Coverage_0']=0
            # df['P_Coverage_S']=0
            # df['EW_Coverage_90']=0
            # df['EW_Coverage_0']=0
            # df['EW_Coverage_S']=0
            # df['P_PanelN_90']=0
            # df['P_PanelN_0']=0
            # df['P_PanelN_S']=0
            # df['EW_PanelN_90']=0
            # df['EW_PanelN_0']=0
            # df['EW_PanelN_S']=0
            # df['P_Tilt_90']=0
            # df['P_Tilt_0']=0
            # df['P_Tilt_S']=0            
            # df['EW_Tilt_90']=0
            # df['EW_Tilt_0']=0
            # df['EW_Tilt_S']=0
            # df['P_Az_90']=0
            # df['P_Az_0']=0
            # df['P_Az_S']=0
            # df['EW_Az_90']=0
            # df['EW_Az_0']=0
            # df['EW_Az_S']=0
            #for each single building do this:
            self.solar_cases=pd.DataFrame(
                columns=['bld_name','latitude','longitude',
                        'bld_azimut','cll_azimut','tilt',
                        'row_dist','N_panel','ratio','cll_layout',
                        'cll_alignment','demand_coverage']
                )
            for bld in self.df_hood.index:    
                ### method to acquire buiding hourly thermal and electricity demand
                ### gets a DF with columns=['elec_demand','thermal_demad']
                ### and values normalized on annual sum
                demand=self.load_bld_demand(bld)
                #### method

                
                
                # case 1: portait parallel to building on short side
                #         PVGIS optimal tilt
                #         minimum interdistance is 0.6m
                roof_short_opt=cp(orientation=float(self.df_hood.loc[bld,'bld_orientation']),
                   lat=float(self.df_hood.loc[bld,'latitude']),
                   long=float(self.df_hood.loc[bld,'longitude']),
                   w=1,l=2,
                   W=float(self.df_hood.loc[bld,'short_side']),
                   L=float(self.df_hood.loc[bld,'long_side']),
                   tilt_EW=20,f_EW=False,
                   f_plot=False,
                   d_rows=0.6,
                   parallel="short",optimal=True,
                   demand=demand.iloc[:,1],irradiance=[self.irr_TMY.ghi,self.irr_TMY.dhi,self.irr_TMY.dni])
                
                self.solar_cases.loc[
                    self.solar_cases.index.size]=[bld,
                                                  float(self.df_hood.loc[bld,'latitude']),
                                                  float(self.df_hood.loc[bld,'longitude']),
                                                  float(self.df_hood.loc[bld,'bld_orientation']),
                                                  roof_short_opt.roof[0,'cll_azimut'],
                                                  roof_short_opt.roof[0,'tilt'],
                                                  roof_short_opt.roof[0,'row_dist'],
                                                  roof_short_opt.roof[0,'N_panel'],
                                                  roof_short_opt.roof[0,'ratio'],
                                                  roof_short_opt.roof[0,'f_EW'],
                                                  "short","demand_cover"
                                                  ]
                
                # case 2: portait parallel to building on long side
                #         PVGIS optimal tilt
                #         minimum interdistance is 0.6m
                roof_long_opt=cp(orientation=float(self.df_hood.loc[bld,'bld_orientation']),
                   lat=float(self.df_hood.loc[bld,'latitude']),
                   long=float(self.df_hood.loc[bld,'longitude']),
                   w=1,l=2,
                   W=float(self.df_hood.loc[bld,'short_side']),
                   L=float(self.df_hood.loc[bld,'long_side']),
                   tilt_EW=20,f_EW=False,
                   f_plot=False,
                   d_rows=0.6,
                   parallel="long",optimal=True,
                   demand=demand.iloc[:,1],irradiance=[self.irr_TMY.ghi,self.irr_TMY.dhi])
                
                self.solar_cases.loc[self.solar_cases.index.size,]
                
                
                
                df.loc[bld,'P_Coverage_90']=roof.roof.loc[0,'ratio']
                df.loc[bld,'P_Coverage_0']=roof.roof.loc[1,'ratio']
                df.loc[bld,'P_Coverage_S']=roof.roof.loc[2,'ratio']
                df.loc[bld,'P_Tilt_90']=roof.roof.loc[0,'tilt']
                df.loc[bld,'P_Tilt_0']=roof.roof.loc[1,'tilt']
                df.loc[bld,'P_Tilt_S']=roof.roof.loc[2,'tilt']
                df.loc[bld,'P_Az_90']=roof.roof.loc[0,'azimut']
                df.loc[bld,'P_Az_0']=roof.roof.loc[1,'azimut']
                df.loc[bld,'P_Az_S']=roof.roof.loc[2,'azimut']
                # case 2: portait parallel to building on long side, short side and
                #         facing South
                #           optimal tilt based on demand
                
                self.optimize_tilt(orient=df.loc[bld,'P_Az_90'],
                                   demand=demand['electricityDemand'],
                                   lat=float(df.loc[bld,'latitude']),
                                   long=float(df.loc[bld,'longitude']),
                                   tilt=df.loc[bld,'P_Tilt_90'])
                #           (minimum interdistance is 0.6m)
                # case 2: EW parallel to building on long side, short side and
                #         facing South
                #            PVGIS optimal tilt
                #           (minimum interdistance is 0.6m)
                
                roof=cp(orientation=float(df.loc[bld,'bld_orientation']),
                   lat=float(df.loc[bld,'latitude']),
                   long=float(df.loc[bld,'longitude']),
                   w=1,l=2,
                   W=float(df.loc[bld,'short_side']),
                   L=float(df.loc[bld,'long_side']),
                   tilt_EW=20,f_EW=True,
                   f_plot=True,
                   d_rows=0.6)
                df.loc[bld,'EW_Coverage_90']=roof.roof.loc[0,'ratio']
                df.loc[bld,'EW_Coverage_0']=roof.roof.loc[1,'ratio']
                df.loc[bld,'EW_Coverage_S']=roof.roof.loc[2,'ratio']
                df.loc[bld,'EW_Tilt_90']=roof.roof.loc[0,'tilt']
                df.loc[bld,'EW_Tilt_0']=roof.roof.loc[1,'tilt']
                df.loc[bld,'EW_Tilt_S']=roof.roof.loc[2,'tilt']
                df.loc[bld,'EW_Az_90']=roof.roof.loc[0,'azimut']
                df.loc[bld,'EW_Az_0']=roof.roof.loc[1,'azimut']
                df.loc[bld,'EW_Az_S']=roof.roof.loc[2,'azimut']
            self.solar_scen=df

            return None
    def optimize_tilt(self,orient=0,
                          demand=None,
                          lat=42,
                          long=6,
                          tilt=0):
        def loss(x,lat,long,orient):
            PVGIS_data = pvlib.iotools.get_pvgis_hourly(lat, long,components=False,
                                                    surface_azimuth=orient,
                                                    start=2016,
                                                    end=2016,
                                                    surface_tilt=x,
                                                    optimal_surface_tilt=False,
                                                    optimalangles=False,
                                                    map_variables=True)
            normal_irr_profile=PVGIS_data[0].poa_global/PVGIS_data[0].poa_global.sum()
            normal_irr_profile=normal_irr_profile.iloc[:demand.index.size].reset_index(drop=True)
            normal_demand_profile=demand/demand.sum()
            normal_coverage_profile=normal_demand_profile
            for cov in normal_coverage_profile.index:
                if normal_irr_profile.loc[cov]>normal_demand_profile.loc[cov]:
                    normal_coverage_profile.loc[cov]=normal_demand_profile.loc[cov]
                else:
                    normal_coverage_profile.loc[cov]=normal_irr_profile.loc[cov]
            irr_cov_sum=normal_coverage_profile.sum()
            return irr_cov_sum
        # args=
        x_bounds = (tilt, 90)
        results_df=pd.DataFrame(columns=['tilt','diff'])
        for x in range(tilt,80,10):
            results_df.loc[results_df.index.size,'diff']=loss(x,lat,long,orient)
            results_df.loc[results_df.index.size-1,'tilt']=x
        
        optimal_x = results_df.tilt[0]
        
                
    def get_TMY(self):
        """
        Method to request TMY data from PVGIS or load data from file 
        """
        if self.load_file==True:
            self.irr_TMY=pd.read_csv(self.target,sep=";")
            dummy=[]
            for tag in self.irr_TMY.columns:
                if 'time' in tag.lower():
                    dummy.append('time.'+tag[-2:])
                else:
                    dummy.append(tag)
            self.irr_TMY.columns=dummy
            self.irr_TMY.index = pd.to_datetime(self.irr_TMY['time.yy'].astype(str) + "-" 
                    + self.irr_TMY['time.mm'].astype(str) + "-" + 
                    self.irr_TMY['time.dd'].astype(str) + " " + self.irr_TMY["time.hh"].astype(str) + ":00:00")
            self.irr_TMY = self.irr_TMY.drop(columns=["time.yy", 
                                                  "time.mm", 
                                                  "time.dd", 
                                                  "time.hh"])
            self.MaxTMYear=self.irr_TMY.index.year.max()
            self.irr_TMY.index = self.irr_TMY.index.map(lambda t: t.replace(year=self.MaxTMYear))
            self.months=[{'month': 1, 'year': self.MaxTMYear}, 
                         {'month': 2, 'year': self.MaxTMYear}, 
                         {'month': 3, 'year': self.MaxTMYear}, 
                         {'month': 4, 'year': self.MaxTMYear}, 
                         {'month': 5, 'year': self.MaxTMYear}, 
                         {'month': 6, 'year': self.MaxTMYear}, 
                         {'month': 7, 'year': self.MaxTMYear}, 
                         {'month': 8, 'year': self.MaxTMYear}, 
                         {'month': 9, 'year': self.MaxTMYear}, 
                         {'month': 10, 'year': self.MaxTMYear}, 
                         {'month': 11, 'year': self.MaxTMYear}, 
                         {'month': 12, 'year': self.MaxTMYear}]
        
        elif self.solar_file==True:
            path=self.solar_f_adr            
            meteo_data=pd.read_csv(path, sep='\t')
            meteo_data.index = pd.to_datetime(meteo_data["time.yy"].astype(str) + "-" + meteo_data["time.mm"].astype(str) + "-" + 
                    meteo_data["time.dd"].astype(str) + " " + meteo_data["time.hh"].astype(str) + ":00:00")
            meteo_data = meteo_data.drop(columns=["stn",
                                                  "time.yy", 
                                                  "time.mm", 
                                                  "time.dd", 
                                                  "time.hh",
                                                  "rre150h0",		
                                                  "fkl010h1",		
                                                  "tso100hs",	
                                                  "nto000sw",				
                                                  "str.vert.E"	,
                                                  "str.vert.S"	,
                                                  "str.vert.W"	,
                                                  "str.vert.N"	,
                                                  "bodenalbedo"	,	
                                                  "ir.vertikal.S"	,
                                                  "bodenemissivitaet",	
                                                  "dewpt"	,
                                                  "enthalpy"	,
                                                  "mixratio"	,
                                                  "wetbulb"])
            meteo_data.columns=(["temp_air"	,
                                 "pressure"	,
                                 "relative_humidity"	,
                                 "wind_speed"	,
                                 "wind_direction"	,
                                 "ghi"	,
                                 "dhi"	,
                                 "dni"	,
                                 "IR(h)"])

            self.irr_TMY=meteo_data
            self.MaxTMYear=self.irr_TMY.index.year.max()
            self.irr_TMY.index = self.irr_TMY.index.map(lambda t: t.replace(year=self.MaxTMYear))
            self.months=[{'month': 1, 'year': self.MaxTMYear}, 
                         {'month': 2, 'year': self.MaxTMYear}, 
                         {'month': 3, 'year': self.MaxTMYear}, 
                         {'month': 4, 'year': self.MaxTMYear}, 
                         {'month': 5, 'year': self.MaxTMYear}, 
                         {'month': 6, 'year': self.MaxTMYear}, 
                         {'month': 7, 'year': self.MaxTMYear}, 
                         {'month': 8, 'year': self.MaxTMYear}, 
                         {'month': 9, 'year': self.MaxTMYear}, 
                         {'month': 10, 'year': self.MaxTMYear}, 
                         {'month': 11, 'year': self.MaxTMYear}, 
                         {'month': 12, 'year': self.MaxTMYear}]
        else:
            PVGIS_output = pvlib.iotools.get_pvgis_tmy(self.lat, self.long,map_variables=True)
            self.irr_TMY=PVGIS_output[0]
            self.MaxTMYear=self.irr_TMY.index.year.max()
            self.irr_TMY.index = self.irr_TMY.index.map(lambda t: t.replace(year=self.MaxTMYear))
            self.months=PVGIS_output[1]
        if hasattr(self, 'DB_soil'):
            self.align_DB()
        return None
    
    def get_hourly(self,st_date="2005",end_date="2005",opt=True):
        """
        Method to request historical data between 2005 and 2016 from PVGIS
        """
        self.opt=opt
        self.start_h=pd.to_datetime(st_date,format="%Y-%m-%d",exact=True)
        self.end_h=pd.to_datetime(end_date,format="%Y-%m-%d",exact=True)
        output = pvlib.iotools.get_pvgis_hourly(self.lat, self.long,start=self.start_h,end=self.end_h,map_variables=True,surface_tilt=self.tilt,optimal_surface_tilt=self.opt,optimalangles=False)
        self.irr_h= output[0] 
        self.months_bis=output[1]
        if hasattr(self, 'DB_soil'):
            self.align_DBh()
        return None
    
    def get_TMY_formatted(self):
        """
        Method to request TMY data from PVGIS and format the output dataframe as wished
        """
        if ~hasattr(self, 'irr_TMY'):          
            self.get_TMY()
        self.irr_TMYf=self.irr_TMY.copy()
        self.irr_TMYf.index.name = "utc_time"
        self.irr_TMYf[self.time_y] = self.MaxTMYear   #irr_TMY.index.year
        self.irr_TMYf[self.time_m] = self.irr_TMYf.index.month
        self.irr_TMYf[self.time_d] = self.irr_TMYf.index.day
        self.irr_TMYf[self.time_h] = self.irr_TMYf.index.hour
        
        if self.load_file==False:
            self.irr_TMYf.rename(columns={"temp_air": self.lb_air, "ghi": self.lb_ghi,"dhi":self.lb_dhi,
                                          "relative_humidity":self.lb_RH,"wind_speed":self.lb_wnd,"pressure":self.lb_p},inplace=True)
            self.good_list=[self.time_y,self.time_m,self.time_d,self.time_h,
                            self.lb_air,self.lb_ghi,self.lb_dhi,self.lb_RH,self.lb_wnd,self.lb_p]
            self.bad_list=['IR(h)','wind_direction']
            self.irr_TMYf.drop(self.bad_list,axis='columns',inplace=True)
            self.irr_TMYf=self.irr_TMYf.loc[:,self.good_list]
        
        if hasattr(self, 'DB_soil'):
            self.align_DB()
        return None
        
    def get_soil(self ,start_date="2021-01-01", end_date="2021-12-31"):
        """
        Method to compute soil probe temperature hourlyprofile

        Parameters
        ----------
        start_date : TYPE, optional
            DESCRIPTION. The default is "2021-01-01".
        end_date : TYPE, optional
            DESCRIPTION. The default is "2021-12-31".

        Returns
        -------
        None.

        """
        self.get_stn()
        head_tail = os.path.split(os.path.abspath(__file__))

        # dati=feather.read_feather(os.path.join(head_tail[0],'DB_soil.fea'))
        dati=feather.read_feather('DB_soil.fea')
        self.Src_soil=dati.loc[dati['stn']==self.stn,['time','Soil_temperature']].copy() 
        self.Src_soil.set_index('time',inplace=True)
        self.soil_year=self.Src_soil.index.year.max()
        self.soil_mean=pd.to_numeric(self.Src_soil.iloc[:,0]).mean()
        self.T_soil_max=self.soil_mean
        self.compute_soil()
        if hasattr(self, 'irr_TMY'):          
            self.align_DB()
        if hasattr(self, 'irr_TMYf'):          
            self.align_DB()
        if hasattr(self, 'irr_h'):
            self.align_DBh()
        return None
    
    def align_DB(self):
        """
        Method to align time index of TMY dataframe and soil temperature
        """
        dummy_index = pd.date_range('2021-01-01 00:00:00', periods=8760, tz='UTC', freq='H')
        self.irr_TMY.index = dummy_index
        self.DB_soil.index=self.irr_TMY.index
        self.irr_TMY[self.lb_gT] =self.DB_soil['Soil_temperature'].copy()
        if hasattr(self, 'irr_TMYf'):  
            self.irr_TMYf.index = dummy_index
            self.irr_TMYf[self.lb_gT] =self.DB_soil['Soil_temperature'].copy()
        return None
    
    def align_DBh(self):
        """
        Method to align time index of historical dataframe and soil temperature
        """
        self.DB_soilh=self.DB_soil.copy()
        self.DB_soilh.index=self.irr_h.index
        self.irr_h[self.lb_gT] =self.DB_soilh['Soil_temperature'].copy()
        return None
    
    def save_DB(self,target):
        """
        Method to save to target
        """
        self.irr_TMYf.to_csv(target,';',index=False)
        return None
    
    def compute_soil(self):
        """
        Method to compute the geothermal probe temerpature during the year as a function of the annual
        soil average temperature in the nearest Meteoshweiz station
        """
        
        #start=-np.pi*(1-self.Day_soil_min*24/8760) #1400 is index for 3 March
        #end=2*np.pi+start
        start=1
        end=8760
        ind=np.linspace(start,end,8760)
        mean=(self.T_soil_min+self.T_soil_max)/2 #minimum probe temperature is 0°C
        amp=(self.T_soil_max-self.T_soil_min)/2 # maximum probe temperature is the yearly average of surface soil temperature
        ground_temp=mean+amp*np.sin(3/2*np.pi+(ind-self.Day_soil_min*24)/8760*2*np.pi)
        self.DB_soil=self.Src_soil.copy()
        self.DB_soil['Soil_temperature']=ground_temp.tolist()
        return None
    
    def get_stn(self):
        """
        Method to locate the nearest soil temperature station
        (it can be improved...centroid of triangle of 3 nearest stations?)
        """
        #print(os.path.abspath(__file__))
        head_tail = os.path.split(os.path.abspath(__file__))
        # DB_geo= pd.read_csv(os.path.join(head_tail[0],"Soil_stn.csv"),encoding='cp1252',delimiter=';',
        #                                 skip_blank_lines=True, 
        #                                 )
        DB_geo= pd.read_csv("Soil_stn.csv",encoding='cp1252',delimiter=';',
                                        skip_blank_lines=True, 
                                        )
        DB_geo=DB_geo.loc[DB_geo['Data source']=='MeteoSchweiz']
        target=geopandas.GeoSeries([Point(self.long, self.lat)])
        target_arr=geopandas.GeoSeries()
        gdf = geopandas.GeoDataFrame(
            DB_geo, geometry=geopandas.points_from_xy(DB_geo.Longitude, DB_geo.Latitude))
        for i in range(0,gdf.index.size):
            target_arr=pd.concat([target_arr,target],axis="index")
        target_arr.index=gdf.index    
        dist=gdf.distance(target_arr)
        station=dist.loc[dist==dist.min()].index[0]
        self.stn=DB_geo.loc[station,'stn']
        return None
    
    def do_clustering(self,clustering_vars):
        self.meteo_daily = self.irr_TMYf.resample('D').mean()
        self.meteo_daily['week_end'] = [1000 if d.weekday() >= 5 else 0 for d in self.meteo_daily.index]
        if clustering_vars==[]:
            clustering_vars= ['tre200h0', 'gls', 'str.diffus', 'ground_temp', 'week_end'] # columns to use for clustering
        clustering_input = self.meteo_daily.loc[:,clustering_vars]
        
        """Normalize input data and perform clustering
        """
        clustering_input_norm=vq.whiten(clustering_input)
        self.meteo_cluster,self.code_BK=vq.kmeans2(clustering_input_norm,self.n_cluster,iter=100,thresh=1e-5,minit="++")
        
        
        """locate nearest days to clusters and compute bin
        """
        labels=[]
        lab_indx=[]
        lab_d=[]
        for i in range(self.n_cluster):
            cl,d=vq.vq(clustering_input_norm,[self.meteo_cluster[i]])
            labels.append(d.argmin())
            lab_indx.append(i)
            lab_d.append(d.min())
            
        """create clustering result table
        """
        self.results = pd.DataFrame(index=self.meteo_daily.index[labels])
        self.results['labels'] = lab_indx
        self.results['count'] = pd.Series(self.code_BK).value_counts().loc[lab_indx].values
        self.results['distances'] = lab_d
        
        return None
    
    
    
if __name__ == "__main__":    
    workbook = xlrd.open_workbook(r"..\data\excels\basic_example\scenario.xls")
    worksheet = workbook.sheet_by_name('solar')
    worksheet.cell(1, worksheet.row_values(0, start_colx=0, end_colx=None).index('latitude')).value=47.55
    worksheet.cell(2, worksheet.row_values(0, start_colx=0, end_colx=None).index('latitude')).value=47.55
    worksheet.cell(3, worksheet.row_values(0, start_colx=0, end_colx=None).index('latitude')).value=47.55
    worksheet.cell(4, worksheet.row_values(0, start_colx=0, end_colx=None).index('latitude')).value=47.55
    worksheet.cell(1, worksheet.row_values(0, start_colx=0, end_colx=None).index('longitude')).value=7.55
    worksheet.cell(2, worksheet.row_values(0, start_colx=0, end_colx=None).index('longitude')).value=7.55
    worksheet.cell(3, worksheet.row_values(0, start_colx=0, end_colx=None).index('longitude')).value=7.55
    worksheet.cell(4, worksheet.row_values(0, start_colx=0, end_colx=None).index('longitude')).value=7.55
    wb=copy(workbook)
    wb.save("scenario.xls")
    meteo=weather(cluster=True,
                  n_clusters=20,
                  source=r"..\data\excels\basic_example\scenario.xls",
                  save_file=False)
    