import oemof.solph as solph
import numpy as np
import pandas as pd
import pvlib
import optihood.combined_prod as cp
from oemof.thermal.solar_thermal_collector import flat_plate_precalc

class SolarCollector(solph.Transformer):
    def __init__(self, label, buildingLabel, inputs, outputs, connector, electrical_consumption, peripheral_losses, latitude,
                 longitude,
                 collector_tilt, roof_area, zenith_angle, collector_azimuth, eta_0, a_1, a_2, temp_collector_inlet,
                 delta_temp_n, irradiance_global,
                 irradiance_diffuse, temp_amb_col, capacityMin, capacityMax, epc, base, env_capa, env_flow, varc, dispatchMode):

        flatPlateCollectorData = flat_plate_precalc(
            latitude, longitude, collector_tilt, collector_azimuth, eta_0, a_1, a_2, temp_collector_inlet, delta_temp_n,
            irradiance_global, irradiance_diffuse, temp_amb_col
        )

        if not (np.isnan(roof_area) or np.isnan(zenith_angle)):
            self.surface_used = self._calculateArea(zenith_angle, collector_tilt, collector_azimuth)
        else:
            self.surface_used = np.nan

        self.collectors_eta_c = flatPlateCollectorData['eta_c']

        self.collectors_heat = flatPlateCollectorData['collectors_heat']/1000 #flow in kWh per m² of solar thermal panel

        if dispatchMode:
            investArgs = {'ep_costs':epc,
                        'minimum':capacityMin,
                        'maximum':capacityMax,
                        'space':self.surface_used,
                        'roof_area':roof_area,
                        'env_per_capa':env_capa}
        else:
            investArgs = {'ep_costs':epc,
                        'minimum':capacityMin,
                        'maximum':capacityMax,
                        'nonconvex':True,
                        'space':self.surface_used,
                        'roof_area':roof_area,
                        'offset':base,
                        'env_per_capa':env_capa}
        self.__collector_source = solph.Source(
            label='heat_'+label + "__" + buildingLabel,
            outputs={
                connector: solph.Flow(
                    fix=self.collectors_heat,
                    investment=solph.Investment(**investArgs),
                     variable_costs=varc,
                     env_per_flow=env_flow,
                )
            },
        )

        self.__collector_excess_heat = solph.Sink(
            label='excess_solarheat' + "__" + buildingLabel, inputs={connector: solph.Flow()}
        )

        self.__collector_transformer = solph.Transformer(
            label=label + '__' + buildingLabel,
            inputs={connector: solph.Flow(), inputs: solph.Flow()},
            outputs={outputs: solph.Flow()},
            conversion_factors={
                connector: 1,
                inputs: electrical_consumption * (1 - peripheral_losses),
                outputs: 1 - peripheral_losses
            },
        )

    def getSolar(self, type):
        if type == 'source':
            return self.__collector_source
        elif type == 'transformer':
            return self.__collector_transformer
        elif type == 'sink':
            return self.__collector_excess_heat
        else:
            print("Transformer label not identified...")
            return []

    def _calculateArea(self, zenith_angle, collector_tilt, collector_azimuth):
        # Coefficient representing the area used by one unit of capacity of a solar collector
        coeff = -np.sin((zenith_angle+collector_tilt)*np.pi/180)*np.cos(collector_azimuth*np.pi/180)/np.sin(zenith_angle*np.pi/180)
        return coeff

class PVT(solph.Transformer):
    def __init__(self, label, buildingLabel, inputs, outputs, connector, electrical_consumption, peripheral_losses, latitude,
                 longitude,
                 collector_tilt, roof_area, zenith_angle, collector_azimuth, eta_0, a_1, a_2, temp_collector_inlet,temp_amb,
                 delta_temp_n, irradiance_global,
                 irradiance_diffuse, capacityMin, capacityMax, epc, base, env_capa, env_flow, varc, dispatchMode,
                 pv_efficiency=0.15, taualpha=0.85):

        pvdata = self.computePvSolarPosition(irradiance_diffuse, irradiance_global, latitude, longitude, collector_azimuth,
                                           collector_tilt, temp_amb)
        pv_electricity = np.minimum(self.pvtElPrecalc(temp_amb, temp_collector_inlet, pvdata['pv_ira'] / 1000, delta_temp_n), capacityMax + base)
        pvtCollectorData = self.pvtThPrecalc(latitude, longitude, collector_tilt, collector_azimuth, eta_0, a_1, a_2, temp_collector_inlet, delta_temp_n,
                                       irradiance_global, irradiance_diffuse, temp_amb,  pv_efficiency, taualpha)
        self.taualpha = taualpha

        self.collectors_heat = pvtCollectorData['collectors_heat']/1000
        self.collectors_eta_c = pvtCollectorData['eta_c']
        self.surface_used = self.calculateArea(zenith_angle, collector_tilt, collector_azimuth)
        surface_used_el = self.calculateArea(zenith_angle, collector_tilt, collector_azimuth, pv_efficiency)
        if dispatchMode:
            investArgs = {'ep_costs':epc,
                         'minimum':capacityMin,
                         'maximum':capacityMax,
                         'space':self.surface_used,
                         'space_el': surface_used_el,
                         'roof_area':roof_area,
                         'env_per_capa':env_capa}
        else:
            investArgs={'ep_costs':epc,
                        'minimum':capacityMin,
                         'maximum':capacityMax,
                         'nonconvex':True,
                         'space':self.surface_used,
                         'space_el': surface_used_el,
                         'roof_area':roof_area,
                         'offset':base,
                         'env_per_capa':env_capa}
        self.__PVTel_source = solph.Source(label='elSource_' + label + '__' + buildingLabel,
                                 outputs={outputs[1]: solph.Flow(
                                     investment=solph.Investment(**investArgs),
                                     variable_costs=varc,
                                     env_per_flow=env_flow,
                                     max=pv_electricity
                                 )}
                                 )
        self.__PVTheat_source = solph.Source(
            label='heatSource_' + label + "__" + buildingLabel,
            outputs={
                connector: solph.Flow(
                    fix=self.collectors_heat,
                    investment=solph.Investment(**investArgs),
                    variable_costs=varc,
                    env_per_flow=env_flow,
                )
            },
        )
        self.__PVT_excessheat = solph.Sink(
            label='excessheat_' + label + "__" + buildingLabel, inputs={connector: solph.Flow()}
        )
        self.__PVTheat_transformer = solph.Transformer(
            label=label + '__' + buildingLabel,
            inputs={connector: solph.Flow(), inputs: solph.Flow()},
            outputs={outputs[0]: solph.Flow()},
            conversion_factors={
                connector: 1,
                inputs: electrical_consumption * (1 - peripheral_losses),
                outputs[0]: 1 - peripheral_losses
            },
        )

    def pvtThPrecalc(self, latitude, longitude, collector_tilt, collector_azimuth, eta_0, a_1, a_2, temp_collector_inlet, delta_temp_n,
            irradiance_global, irradiance_diffuse, temp_amb,  pv_efficiency, taualpha):
        # function inspired from oemof thermal solar_thermal_collector
        # The calculation of eta_c is overriden
        data = pd.DataFrame(
            {
                'ghi': irradiance_global,
                'dhi': irradiance_diffuse,
                'temp_amb': temp_amb
            }
        )
        solposition = pvlib.solarposition.get_solarposition(
            time=data.index, latitude=latitude, longitude=longitude
        )
        dni = pvlib.irradiance.dni(
            ghi=data['ghi'], dhi=data['dhi'], zenith=solposition['apparent_zenith']
        )
        total_irradiation = pvlib.irradiance.get_total_irradiance(
            surface_tilt=collector_tilt,
            surface_azimuth=collector_azimuth,
            solar_zenith=solposition['apparent_zenith'],
            solar_azimuth=solposition['azimuth'],
            dni=dni.fillna(0),  # fill NaN values with '0'
            ghi=data['ghi'],
            dhi=data['dhi'],
        )
        data['col_ira'] = total_irradiation['poa_global']
        data["eta_c"] = self.calcEtaC(eta_0, a_1, a_2, temp_collector_inlet, delta_temp_n, data['temp_amb'], total_irradiation['poa_global'], pv_efficiency, taualpha)
        collectors_heat = data["eta_c"] * total_irradiation['poa_global']
        data["collectors_heat"] = collectors_heat
        return data

    def calcEtaC(self, eta_0, a_1, a_2, temp_collector_inlet, delta_temp_n, temp_amb, irradiation, pv_efficiency, taualpha):
        delta_t = temp_collector_inlet + delta_temp_n - temp_amb
        eta_c = pd.Series()
        for index, value in irradiation.items():
            if value > 0:
                eta = (eta_0*(1-pv_efficiency/(taualpha))
                        - a_1 * delta_t[index] / value
                        - a_2 * delta_t[index] ** 2 / value)
                eta_c[index] = eta*(eta > 0) + 0
            else:
                eta_c[index] = 0
        return eta_c

    def computePvSolarPosition(self, irradiance_diffuse, irradiance_global, latitude, longitude, pv_azimuth, pv_tilt,
                               temp_amb_pv):
        data = pd.DataFrame(
            {
                'ghi': irradiance_global,
                'dhi': irradiance_diffuse,
                'temp_amb': temp_amb_pv
            }
        )
        solposition = pvlib.solarposition.get_solarposition(
            time=data.index, latitude=latitude, longitude=longitude
        )
        dni = pvlib.irradiance.dni(
            ghi=data['ghi'], dhi=data['dhi'], zenith=solposition['apparent_zenith']
        )
        total_irradiation = pvlib.irradiance.get_total_irradiance(
            surface_tilt=pv_tilt,
            surface_azimuth=pv_azimuth,
            solar_zenith=solposition['apparent_zenith'],
            solar_azimuth=solposition['azimuth'],
            dni=dni.fillna(0),  # fill NaN values with '0'
            ghi=data['ghi'],
            dhi=data['dhi'],
        )
        data['pv_ira'] = total_irradiation['poa_global']
        return data

    # model according to TriHP Energy management algorithms description. D6.2 v1.0
    # nominal power is 1 KW due to the normed optimizer
    def pvtElPrecalc(self, temp_amb, temp_collector_inlet, i_H_t,delta_temp_n, a4=-0.062059,
                   a5=0.04277774, a6=9.692792, a7=-1.885868, a8=6.6):
        temp_cell = temp_collector_inlet + delta_temp_n - temp_amb
        pvPower = np.maximum(0, (a4*i_H_t + a5)*temp_cell + a6*i_H_t + a7) / a8
        return pvPower


    def calculateArea(self, zenith_angle, pv_tilt, pv_azimuth, pv_efficiency=0.0):
        # Coefficient representing the area used by one unit of capacity (kW thermal) of a PVT panel
        coeff = -np.sin((zenith_angle+pv_tilt)*np.pi/180)*np.cos(pv_azimuth*np.pi/180)/np.sin(zenith_angle*np.pi/180)
        if pv_efficiency:
            coeff /= pv_efficiency
        return coeff

    def getPVT(self, type):
        if type == 'heat_source':
            return self.__PVTheat_source
        elif type == 'heat_transformer':
            return self.__PVTheat_transformer
        elif type == 'el_source':
            return self.__PVTel_source
        elif type == 'excess_heat_sink':
            return self.__PVT_excessheat
        else:
            print("Label not identified...")
            return []

class HeatPumpLinear:
    "Information about the model can be found in combined_pro.py CombinedTransformer"
    def __init__(self, buildingLabel, temperatureDHW, temperatureSH, temperatureLow, input, outputSH, outputDHW,
                 capacityMin, capacityMax, nomEff,
                 epc, base, varc, env_flow, env_capa, dispatchMode):
        self.__copDHW = self._calculateCop(temperatureDHW, temperatureLow)
        self.__copSH = self._calculateCop(temperatureSH, temperatureLow)
        self.avgCopSh = (sum(self.__copSH)/len(self.__copSH))
        self.nominalEff = nomEff
        if dispatchMode:
            investArgs = {'ep_costs' : epc * nomEff,
            'minimum' : capacityMin / nomEff,
            'maximum' : capacityMax / nomEff,
            'env_per_capa' : env_capa * nomEff}
        else:
            investArgs = {'ep_costs' : epc * nomEff,
            'minimum' : capacityMin / nomEff,
            'maximum' : capacityMax / nomEff,
            'nonconvex' : True,
            'offset' : base,
            'env_per_capa' : env_capa * nomEff}
        self.__heatpump = cp.CombinedTransformer(label='HP' + '__' + buildingLabel,
                                            inputs={input: solph.Flow(
                                                investment=solph.Investment(**investArgs),
                                            )},
                                            outputs={outputSH: solph.Flow(
                                                          variable_costs=varc,
                                                          env_per_flow=env_flow,
                                                      ),
                                                      outputDHW: solph.Flow(
                                                          variable_costs=varc,
                                                          env_per_flow=env_flow,
                                                      )
                                                  },
                                            efficiencies={outputSH: self.__copSH,
                                                        outputDHW: self.__copDHW})

    def _calculateCop(self, tHigh, tLow):
        coefW = [0.66610, -2.2365, 15.541, 25.705, -17.407, 3.8145]
        coefQ = [11.833, 96.504, 14.496, -50.064, 161.02, -133.60]
        QCondenser = coefQ[0] + (coefQ[1] * tLow / 273.15) + (coefQ[2] * tHigh / 273.15) + (
                coefQ[3] * tLow / 273.15 * tHigh / 273.15) + (
                             coefQ[4] * (tLow / 273.15) ** 2) + (
                             coefQ[5] * (tHigh / 273.15) ** 2)
        WCompressor = coefW[0] + (coefW[1] * tLow / 273.15) + (coefW[2] * tHigh / 273.15) + (
                coefW[3] * tLow / 273.15 * tHigh / 273.15) + (
                             coefW[4] * (tLow / 273.15) ** 2) + (
                              coefW[5] * (tHigh / 273.15) ** 2)
        cop = np.divide(QCondenser, WCompressor)
        return cop

    def getHP(self, type):
        if type == 'sh':
            return self.__heatpump
        else:
            print("Transformer label not identified...")
            return []


class GeothermalHeatPumpLinear:
    "Information about the model can be found in combined_pro.py CombinedTransformer"
    def __init__(self, buildingLabel, temperatureDHW, temperatureSH, temperatureLow, input, outputSH, outputDHW,
                 capacityMin, capacityMax, nomEff,
                 epc, base, varc, env_flow, env_capa, dispatchMode):
        self.__copDHW = self._calculateCop(temperatureDHW, temperatureLow)
        self.__copSH = self._calculateCop(temperatureSH, temperatureLow)
        self.avgCopSh = (sum(self.__copSH)/len(self.__copSH))
        self.nominalEff = nomEff
        if dispatchMode:
            investArgs= {'ep_costs':epc*nomEff,
                        'minimum':capacityMin/nomEff,
                        'maximum':capacityMax/nomEff,
                        'env_per_capa':env_capa*nomEff,
                    }
        else:
            investArgs= {'ep_costs':epc*nomEff,
                        'minimum':capacityMin/nomEff,
                        'maximum':capacityMax/nomEff,
                        'nonconvex':True,
                        'offset':base,
                        'env_per_capa':env_capa*nomEff,
                    }
        self.__geothermalheatpump = cp.CombinedTransformer(label='GWHP' + '__' + buildingLabel,
                                            inputs={input: solph.Flow(
                                                investment=solph.Investment(**investArgs),
                                            )},
                                            outputs={outputSH: solph.Flow(
                                                          variable_costs=varc,
                                                          env_per_flow=env_flow,
                                                      ),
                                                      outputDHW: solph.Flow(
                                                          variable_costs=varc,
                                                          env_per_flow=env_flow,
                                                      )
                                                  },
                                            efficiencies={outputSH: self.__copSH,
                                                        outputDHW: self.__copDHW})

    def _calculateCop(self, tHigh, tLow):
        coefW = [0.1600, -1.2369, 19.9391, 19.3448, 7.1057, -1.4048]
        coefQ = [13.8978, 114.8358, -9.3634, -179.4227, 342.3363, -12.4969]
        QCondenser = coefQ[0] + (coefQ[1] * tLow / 273.15) + (coefQ[2] * tHigh / 273.15) + (
                coefQ[3] * tLow / 273.15 * tHigh / 273.15) + (
                             coefQ[4] * (tLow / 273.15) ** 2) + (
                             coefQ[5] * (tHigh / 273.15) ** 2)
        WCompressor = coefW[0] + (coefW[1] * tLow / 273.15) + (coefW[2] * tHigh / 273.15) + (
                coefW[3] * tLow / 273.15 * tHigh / 273.15) + (
                             coefW[4] * (tLow / 273.15) ** 2) + (
                              coefW[5] * (tHigh / 273.15) ** 2)
        cop = np.divide(QCondenser, WCompressor)
        return cop

    def getHP(self, type):
        if type == 'sh':
            return self.__geothermalheatpump
        else:
            print("Transformer label not identified...")
            return []


class CHP:
    "Information about the model can be found in combined_pro.py CombinedCHP"
    def __init__(self, buildingLabel, input, outputEl, outputSH, outputDHW, efficiencyEl, efficiencySH, efficiencyDHW,
                 capacityMin, capacitySH, epc, base, varc1, varc2, env_flow1, env_flow2, env_capa, timesteps, dispatchMode):
        self._efficiencyEl = [efficiencyEl] * timesteps
        self._efficiencySH = [efficiencySH] * timesteps
        self._efficiencyDHW = [efficiencyDHW] * timesteps
        self.avgEff = efficiencySH
        if dispatchMode:
            investArgs = {'ep_costs': epc * self.avgEff,
                          'minimum': capacityMin / self.avgEff,
                          'maximum': capacitySH / self.avgEff,
                          'env_per_capa': env_capa * self.avgEff}
        else:
            investArgs={'ep_costs':epc*self.avgEff,
                        'minimum':capacityMin/self.avgEff,
                        'maximum':capacitySH/self.avgEff,
                        'nonconvex':True,
                        'offset':base,
                        'env_per_capa':env_capa*self.avgEff}
        self.__CHP = cp.CombinedCHP(
                        label='CHP'+'__'+buildingLabel,
                        inputs={
                            input: solph.Flow(
                                investment=solph.Investment(**investArgs),
                            )
                        },
                        outputs={
                            outputSH: solph.Flow(
                                variable_costs=varc2,
                                env_per_flow=env_flow2,
                            ),
                            outputDHW: solph.Flow(
                                variable_costs=varc2,
                                env_per_flow=env_flow2,
                            ),
                            outputEl: solph.Flow(
                                variable_costs=varc1,
                                env_per_flow=env_flow1,
                            ),
                        },
                        efficiencies={outputSH: self._efficiencySH,
                                      outputDHW: self._efficiencyDHW,
                                      outputEl: self._efficiencyEl,
                                      }
                    )

    def getCHP(self, type):
        if type == 'sh':
            return self.__CHP
        else:
            print("Transformer label not identified...")
            return []

class GasBoiler(cp.CombinedTransformer):
    "Information about the model can be found in combined_pro.py CombinedTransformer"
    def __init__(self, buildingLabel, input, outputSH, outputDHW, efficiencySH, efficiencyDHW,
                 capacityMin, capacityMax, epc, base, varc, env_flow, env_capa, dispatchMode):
        self.__efficiency = efficiencySH
        if dispatchMode:
            investArgs= {'ep_costs':epc*efficiencySH,
                        'minimum':capacityMin/efficiencySH,
                        'maximum':capacityMax/efficiencySH,
                        'env_per_capa':env_capa*efficiencySH}
        else:
            investArgs= {'ep_costs':epc*efficiencySH,
                        'minimum':capacityMin/efficiencySH,
                        'maximum':capacityMax/efficiencySH,
                        'nonconvex':True,
                        'offset':base,
                        'env_per_capa':env_capa*efficiencySH}
        super(GasBoiler, self).__init__(
            label='GasBoiler'+'__'+buildingLabel,
            inputs={
                input: solph.Flow(
                    investment=solph.Investment(**investArgs),
                )
            },
            outputs={
                outputSH: solph.Flow(
                    variable_costs=varc,
                    env_per_flow=env_flow
                ),
                outputDHW: solph.Flow(
                    variable_costs=varc,
                    env_per_flow=env_flow
                ),
            },
            efficiencies={outputSH: efficiencySH,
                          outputDHW: efficiencyDHW}
                 )

class ElectricRod(cp.CombinedTransformer):
    def __init__(self, buildingLabel, input, outputSH, outputDHW, efficiency,
                 capacityMin, capacityMax, epc, base, varc, env_flow, env_capa, dispatchMode):

        self.__efficiency = efficiency
        if dispatchMode:
            investArgs = {'ep_costs': epc * efficiency,
                          'minimum': capacityMin / efficiency,
                          'maximum': capacityMax / efficiency,
                          'env_per_capa': env_capa * efficiency}
        else:
            investArgs={'ep_costs':epc*efficiency,
                        'minimum':capacityMin/efficiency,
                        'maximum':capacityMax/efficiency,
                        'nonconvex':True,
                        'offset':base,
                        'env_per_capa':env_capa*efficiency}
        super(ElectricRod, self).__init__(
            label='ElectricRod'+'__'+buildingLabel,
            inputs={
                input: solph.Flow(investment=solph.Investment(**investArgs) )
            },
            outputs={
                outputSH: solph.Flow(
                    variable_costs=varc,
                    env_per_flow=env_flow,
                ),
                outputDHW: solph.Flow(
                    variable_costs=varc,
                    env_per_flow=env_flow,
                ),
            },
            efficiencies={outputSH: efficiency,
                          outputDHW: efficiency}
                 )

class GeothermalHeatPumpLinearSingleUse:
    "Class implementing a linear model for geothermal heat pump for single use only (either DHW or SH but not both)"
    def __init__(self, buildingLabel, temperatureH, temperatureLow, input, outputH,
                 capacityMin, capacityMax, epc, base, varc, env_flow, env_capa, dispatchMode):
        self.__copH = self._calculateCop(temperatureH, temperatureLow)
        #self.avgCopSh = (sum(self.__copH)/len(self.__copH))
        #self.nominalEff = nomEff
        if dispatchMode:
            investArgs = {'ep_costs':epc,
                            'minimum':0,
                            'maximum':capacityMax,
                            'env_per_capa':env_capa}
        else:
            investArgs = {'ep_costs':epc,
                            'minimum':0,
                            'maximum':capacityMax,
                            'nonconvex':True,
                            'offset':base,
                            'env_per_capa':env_capa}
        self.__geothermalheatpump = solph.Transformer(label=f'GWHP{str(temperatureH)}' + '__' + buildingLabel,
                                            inputs={input: solph.Flow()},
                                            outputs={outputH: solph.Flow(
                                                          variable_costs=varc,
                                                          env_per_flow=env_flow,
                                                          investment=solph.Investment(**investArgs),
                                                      )
                                                  },
                                            conversion_factors={outputH: self.__copH})

    def _calculateCop(self, tHigh, tLow):
        coefW = [0.1600, -1.2369, 19.9391, 19.3448, 7.1057, -1.4048]
        coefQ = [13.8978, 114.8358, -9.3634, -179.4227, 342.3363, -12.4969]
        QCondenser = coefQ[0] + (coefQ[1] * tLow / 273.15) + (coefQ[2] * tHigh / 273.15) + (
                coefQ[3] * tLow / 273.15 * tHigh / 273.15) + (
                             coefQ[4] * (tLow / 273.15) ** 2) + (
                             coefQ[5] * (tHigh / 273.15) ** 2)
        WCompressor = coefW[0] + (coefW[1] * tLow / 273.15) + (coefW[2] * tHigh / 273.15) + (
                coefW[3] * tLow / 273.15 * tHigh / 273.15) + (
                             coefW[4] * (tLow / 273.15) ** 2) + (
                              coefW[5] * (tHigh / 273.15) ** 2)
        cop = np.divide(QCondenser, WCompressor)
        return cop

    def getHP(self, type):
        if type == 'sh':
            return self.__geothermalheatpump
        else:
            print("Transformer label not identified...")
            return []