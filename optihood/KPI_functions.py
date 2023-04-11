
import os
import pandas as pd
from optihood.plot_functions import getData
from optihood.labelDict import labelDict
from optihood.labelDict import labelPositionDict
from optihood.labelDict import fullPositionDict
import numpy as np
import matplotlib.pyplot as plt
import math

# plt.rcParams['figure.figsize'] = [8.0, 8.0]
plt.rcParams['figure.dpi'] = 600

def read_results(resultFileName):
    dataDict = getData(resultFileName)
    keys=dataDict.keys()
    return dataDict


def app_labeldict(labelDict, test_str):
    res = test_str.replace("(", "").replace(")", "").replace("'", "").replace('"', '')
    res = list(map(str, res.split(', ')))
    res = [labelDict.get(item,item)  for item in res]
    return res


def heat_gen(dataDict, buildings, timeStep):
    """
    Function to calculate total heat generation per technology
    :param :
            dataDict: full result file
            buildings: number of buildings
            timeStep: type str hour, month or year

    :return: heatTec: dataframe with total amount of heat generated per technology and selected time-step
    """

    if timeStep == 'year':
        shHeatTec = pd.DataFrame()
        dhwHeatTec = pd.DataFrame()
    elif timeStep == 'month':
        heatTec = pd.DataFrame()
    elif timeStep == 'hour':
        heatTec = pd.DataFrame()


    for b in range(1, buildings+1):
        if "(('ElectricRod__Building" + str(b) + "', 'shSourceBus__Building" + str(b) + "'), 'flow')" in list(dataDict["shSourceBus__Building" + str(b)]):
            shGenList = ["(('CHP__Building" + str(b) + "', 'shSourceBus__Building" + str(b) + "'), 'flow')",
                         "(('GWHP__Building" + str(b) + "', 'shSourceBus__Building" + str(b) + "'), 'flow')",
                         "(('GasBoiler__Building" + str(b) + "', 'shSourceBus__Building" + str(b) + "'), 'flow')",
                         "(('HP__Building" + str(b) + "', 'shSourceBus__Building" + str(b) + "'), 'flow')",
                         "(('ElectricRod__Building" + str(b) + "', 'shSourceBus__Building" + str(b) + "'), 'flow')"]

            dhwGenList = ["(('CHP__Building" + str(b) + "', 'dhwStorageBus__Building" + str(b) + "'), 'flow')",
                          "(('GWHP__Building" + str(b) + "', 'dhwStorageBus__Building" + str(b) + "'), 'flow')",
                          "(('GasBoiler__Building" + str(b) + "', 'dhwStorageBus__Building" + str(b) + "'), 'flow')",
                          "(('HP__Building" + str(b) + "', 'dhwStorageBus__Building" + str(b) + "'), 'flow')",
                          "(('ElectricRod__Building" + str(b) + "', 'dhwStorageBus__Building" + str(b) + "'), 'flow')",
                          "(('solarCollector__Building" + str(b) + "', 'dhwStorageBus__Building" + str(b) + "'), 'flow')"]
            SHUnitIndex = ['CHPsh', "GWHPsh", "GasBoilersh", "HPsh", "ElectricRodsh"]
            DHWUnitIndex = ['CHPdhw', "GWHPdhw", "GasBoilerdhw", "HPdhw", "ElectricRoddhw", "solarCollectordhw"]
        else:
            shGenList = ["(('CHP__Building" + str(b) + "', 'shSourceBus__Building" + str(b) + "'), 'flow')",
                         "(('GWHP__Building" + str(b) + "', 'shSourceBus__Building" + str(b) + "'), 'flow')",
                         "(('GasBoiler__Building" + str(b) + "', 'shSourceBus__Building" + str(b) + "'), 'flow')",
                         "(('HP__Building" + str(b) + "', 'shSourceBus__Building" + str(b) + "'), 'flow')"]

            dhwGenList = ["(('CHP__Building" + str(b) + "', 'dhwStorageBus__Building" + str(b) + "'), 'flow')",
                          "(('GWHP__Building" + str(b) + "', 'dhwStorageBus__Building" + str(b) + "'), 'flow')",
                          "(('GasBoiler__Building" + str(b) + "', 'dhwStorageBus__Building" + str(b) + "'), 'flow')",
                          "(('HP__Building" + str(b) + "', 'dhwStorageBus__Building" + str(b) + "'), 'flow')",
                          "(('solarCollector__Building" + str(b) + "', 'dhwStorageBus__Building" + str(b) + "'), 'flow')"]
            SHUnitIndex = ['CHPsh', "GWHPsh", "GasBoilersh", "HPsh"]
            DHWUnitIndex = ['CHPdhw', "GWHPdhw", "GasBoilerdhw", "HPdhw", "solarCollectordhw"]

        if timeStep == "year":
            shHeatTec[labelDict["shSourceBus__Building"+str(b)]] = dataDict["shSourceBus__Building" + str(b)][shGenList].sum(axis=0).reset_index(drop=True)
            dhwHeatTec[labelDict["domesticHotWaterBus__Building" + str(b)]] = dataDict["domesticHotWaterBus__Building" + str(b)][dhwGenList].sum(axis=0).reset_index(drop=True)

        elif timeStep == "month":
            heatTec = pd.concat([heatTec,
                                    dataDict["shSourceBus__Building" + str(b)][shGenList].groupby(pd.Grouper(freq='M')).sum()
                                    ], axis=1)
            heatTec = pd.concat([heatTec,
                                    dataDict["domesticHotWaterBus__Building" + str(b)][dhwGenList].groupby(pd.Grouper(freq='M')).sum()
                                    ], axis=1)

        elif timeStep == "hour":
            heatTec = pd.concat([heatTec,
                                 dataDict["shSourceBus__Building" + str(b)][shGenList]], axis=1)
            heatTec = pd.concat([heatTec,
                             dataDict["domesticHotWaterBus__Building" + str(b)][dhwGenList]
                             ], axis=1)

    if timeStep == "year":
        shHeatTec = shHeatTec.set_axis(SHUnitIndex, axis=0)
        dhwHeatTec = dhwHeatTec.set_axis(DHWUnitIndex, axis=0)
        shHeatTec['shTotal'] = shHeatTec.sum(axis = 1)
        dhwHeatTec['dhwTotal'] = dhwHeatTec.sum(axis = 1)
        heatTec = pd.concat([shHeatTec, dhwHeatTec], axis =1).fillna(0).astype(float)

        heatTec['total'] = heatTec['shTotal'] + heatTec['dhwTotal']

    elif timeStep == 'month' or timeStep == 'hour':
        heatTec = heatTec.rename(columns=labelDict)

        chpsh_columns = []
        GWHPsh_columns = []
        GasBoilersh_columns = []
        HPsh_columns = []
        ERsh_columns = []
        chpdhw_columns = []
        GWHPdhw_columns = []
        GasBoilerdhw_columns = []
        HPdhw_columns = []
        ERdhw_columns = []
        solarCollectordhw_columns = []
        for b in range(1, buildings+1):
            chpsh_columns.append('CHPshSourceBus_B'+str(b))
            GWHPsh_columns.append('GWHPshSourceBus_B'+str(b))
            GasBoilersh_columns.append('GasBoilershSourceBus_B'+str(b))
            HPsh_columns.append('HPshSourceBus_B'+str(b))
            try:
                ERsh_columns.append('ElectricRodshSourceBus_B'+str(b))
                ERdhw_columns.append('ElectricRoddhwStorageBus_B'+str(b))
            except:
                None
            chpdhw_columns.append('CHPdhwStorageBus_B'+str(b))
            GWHPdhw_columns.append('GWHPdhwStorageBus_B'+str(b))
            GasBoilerdhw_columns.append('GasBoilerdhwStorageBus_B'+str(b))
            HPdhw_columns.append('HPdhwStorageBus_B'+str(b))

            solarCollectordhw_columns.append('solarCollectordhw_B'+str(b))

            for tec in ["CHP", "GWHP", "GasBoiler", "HP", "ElectricRod"]:
                try:
                    heatTec[str(tec) + str(b)] = heatTec[str(tec) + 'shSourceBus_B' +str(b)] + heatTec[str(tec) + 'dhwStorageBus_B'+str(b)]
                except:
                    None
            # try:
            #     heatTec[str("CHP") + str(b)] = heatTec[str("CHP") + 'shSourceBus_B' +str(b)] + heatTec[str("CHP") + 'dhwStorageBus_B'+str(b)]
            # except:
            #     None
            # try:
            #     heatTec[str("ElectricRod") + str(b)] = heatTec[str("ElectricRod") + 'shSourceBus_B' +str(b)] + heatTec[str("ElectricRod") + 'dhwStorageBus_B'+str(b)]
            # except:
            #     None

        heatTec['GWHPshTotal'] = heatTec[GWHPsh_columns].sum(axis=1)
        heatTec['GasBoilershTotal'] = heatTec[GasBoilersh_columns].sum(axis=1)
        heatTec['HPshTotal'] = heatTec[HPsh_columns].sum(axis=1)
        heatTec['GWHPdhwTotal'] = heatTec[GWHPdhw_columns].sum(axis=1)
        heatTec['GasBoilerdhwTotal'] = heatTec[GasBoilerdhw_columns].sum(axis=1)
        heatTec['HPdhwTotal'] = heatTec[HPdhw_columns].sum(axis=1)
        heatTec['solarCollectordhwTotal'] = heatTec[solarCollectordhw_columns].sum(axis=1)
        heatTec['GWHPTotal'] = heatTec['GWHPshTotal'] + heatTec['GWHPdhwTotal']
        heatTec['GasBoilerTotal'] = heatTec['GasBoilershTotal'] + heatTec['GasBoilerdhwTotal']
        heatTec['HPTotal'] = heatTec['HPshTotal'] + heatTec['HPdhwTotal']
        heatTec['solarCollectorTotal'] = heatTec['solarCollectordhwTotal']

        try:
            heatTec['CHPshTotal'] = heatTec[chpsh_columns].sum(axis=1)
            heatTec['CHPdhwTotal'] = heatTec[chpdhw_columns].sum(axis=1)
            heatTec['CHPTotal'] = heatTec['CHPshTotal'] + heatTec['CHPdhwTotal']
        except:
            None

        try:
            heatTec['ElectricRodshTotal'] = heatTec[ERsh_columns].sum(axis=1)
            heatTec['ElectricRoddhwTotal'] = heatTec[ERdhw_columns].sum(axis=1)
            heatTec['ElectricRodTotal'] = heatTec['ElectricRodshTotal'] + heatTec['ElectricRoddhwTotal']
        except:
            None

        if timeStep == 'month':
            heatTec = heatTec.set_axis(['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec'], axis=0)
        elif timeStep == 'hour':
            heatTec.index = heatTec.index #+ pd.offsets.DateOffset(years=3)

    return heatTec


def elec_sell(dataDict, buildings, timeStep):
    """
    Function to calculate electricity sold
    :param :
            dataDict: full result file
            buildings: number of buildings
            timeStep: type str hour, month or year

    :return: elecSell: dataframe with electricity sold for selected time-step
    """

    for b in range(1, buildings+1):
        if b == 1:
            elecSell = dataDict["electricityBus__Building" + str(b)]
        elif b != 1:
            elecSell = pd.concat([elecSell, dataDict["electricityBus__Building" + str(b)]], axis=1)


    elecSell = elecSell.rename(columns=labelDict)
    elecSell = elecSell[["excessElec_B"+str(b) for b in range(1, buildings+1)]]

    if timeStep == 'year':
        elecSell = elecSell .sum(axis=0).reset_index(drop=True)
    elif timeStep == 'month':
        elecSell = elecSell.groupby(pd.Grouper(freq='M')).sum()
        elecSell = elecSell.set_axis(['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec'], axis=0)
        elecSell['excessTotal'] = elecSell.sum(axis=1)
    elif timeStep == 'hour':
        elecSell['excessTotal'] = elecSell.sum(axis=1)
        elecSell.index = elecSell.index

    return elecSell


def elec_gen(dataDict, buildings, timeStep):
    """
    Function to calculate total elec generation per technology
    :param :
            dataDict: full result file
            buildings: number of buildings
            timeStep: type str hour, month or year

    :return: elecTec: dataframe with total amount of heat generated per technology and time step
    """

    elecGenTec = pd.DataFrame()
    elecGrid = []

    for b in range(1, buildings+1):
        elecGenList = ["(('CHP__Building" + str(b) + "', 'electricityProdBus__Building" + str(b) + "'), 'flow')",
                       "(('pv__Building" + str(b) + "', 'electricityProdBus__Building" + str(b) + "'), 'flow')"]

        if timeStep == "year":
            elecGenTec[labelDict["electricityProdBus__Building" + str(b)]] = dataDict["electricityProdBus__Building" + str(b)][elecGenList].sum(axis=0).reset_index(drop=True)
            elecGrid.append(dataDict["gridBus__Building" + str(b)][
                "(('electricityResource__Building" + str(b) + "', 'gridBus__Building" + str(b) + "'), 'flow')"].sum())

        elif timeStep == "month":
            elecGenTec = pd.concat([elecGenTec,
                                    dataDict["electricityProdBus__Building" + str(b)][elecGenList].groupby(pd.Grouper(freq='M')).sum()
                                    ], axis=1)
            elecGenTec = pd.concat([elecGenTec,dataDict["gridBus__Building" + str(b)].loc[:,
                                     "(('electricityResource__Building" + str(b) + "', 'gridBus__Building" + str(b) + "'), 'flow')"
                                                ].groupby(pd.Grouper(freq='M')).sum()
                                    ], axis=1)
        elif timeStep == 'hour':
            elecGenTec = pd.concat([elecGenTec,
                                    dataDict["electricityProdBus__Building" + str(b)][elecGenList]], axis=1)
            elecGenTec = pd.concat([elecGenTec, dataDict["gridBus__Building" + str(b)].loc[:,
                                   "(('electricityResource__Building" + str(b) + "', 'gridBus__Building" + str(b) + "'), 'flow')"
                                   ]], axis=1)

    if timeStep == "year":
        elecGenTec.loc["Grid"] = elecGrid
        elecGenTec = elecGenTec.set_axis(['CHP', "PV", "Grid"
                                          ], axis=0)
        elecGenTec['total'] = elecGenTec.sum(axis=1)

    elif timeStep == "month" or timeStep == "hour":
        if timeStep == "month":
            elecGenTec = elecGenTec.set_axis(['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec'], axis=0)
        else:
            None

        columns = []
        chp_columns = []
        pv_columns = []
        grid_columns = []
        for b in range(1,buildings+1):
            columns = columns + ['CHP_B'+str(b), "PV_B"+str(b), "Grid_B"+str(b)]
            chp_columns.append('CHP_B'+str(b))
            pv_columns.append('PV_B'+str(b))
            grid_columns.append('Grid_B'+str(b))
        elecGenTec.columns = columns

        elecGenTec['CHPTotal'] = elecGenTec[chp_columns].sum(axis=1)
        elecGenTec['PVTotal'] = elecGenTec[pv_columns].sum(axis=1)
        elecGenTec['gridTotal'] = elecGenTec[grid_columns].sum(axis=1)
        elecGenTec['elecTotal'] = elecGenTec['CHPTotal'] + elecGenTec['PVTotal'] + elecGenTec['gridTotal']


        if timeStep == 'hour':
            elecGenTec.index = elecGenTec.index #+ pd.offsets.DateOffset(years=3)

    return elecGenTec


def ElecInfluenceBuilding(dataDict, buildings, PVImpact, CHPImpact, gridImpact, sel_b, parameter):
    """
    Function to create new electrcity price or co2 impact based on a pv and grid mix
    :param :
            dataDict: full result file
            buildings: number of buildings
            PVImpact, CHPImpact: levelized impact of electricity generated
            gridImpact: Impact of the grid (for co2 co2 per hour for the grid
            sel_b: selected building
            parameter: co2 or cost

    :return:
    """
    elecTec = elec_gen(dataDict, buildings, 'hour')
    elecImpact = pd.DataFrame(index = elecTec.index)

    if parameter == 'cost':
        if type(sel_b) is int:
            elecImpact['ElecPrice' + str(sel_b)] = np.where(elecTec['PV_B'+str(sel_b)] + elecTec['Grid_B'+str(sel_b)] > 0,
                                                            (
                                                             elecTec['PV_B'+str(sel_b)] * PVImpact
                                                           + elecTec['CHP_B'+str(sel_b)] * CHPImpact
                                                           + elecTec['Grid_B'+str(sel_b)] * gridImpact ) / (
                                                             elecTec['PV_B'+str(sel_b)] + elecTec['CHP_B'+str(sel_b)]
                                                           + elecTec['Grid_B'+str(sel_b)]),
                                                            gridImpact)

        elif sel_b == 'system':
            pd.set_option('display.max_rows', None)

            elecImpact['ElecPriceSystem'] = np.where(elecTec['PVTotal'] + elecTec['gridTotal'] > 0,
                                                    (elecTec['PVTotal'] * PVImpact
                                                   + elecTec['CHPTotal'] * CHPImpact
                                                   + elecTec['gridTotal'] * gridImpact) / (
                                                     elecTec['PVTotal'] + elecTec['CHPTotal'] + elecTec['gridTotal']),
                                                                                          gridImpact)

        else:
            print('false definition of sel_b')
            raise NotImplementedError

    elif parameter == 'co2':
        # PVimpact : co2 impact per kWh
        # gridImpact : co2 impact per kWh
        if type(sel_b) is int:
            elecImpact['Co2Impact_B' + str(sel_b)] = np.where(elecTec['PV_B'+str(sel_b)] + elecTec['Grid_B'+str(sel_b)] > 0,
                                                     ((elecTec['PV_B'+str(sel_b)] * PVImpact)/(elecTec['PV_B'+str(sel_b)].sum())
                                                    + elecTec['Grid_B'+str(sel_b)] * gridImpact) / (
                                                      elecTec['PV_B'+str(sel_b)]  + elecTec['Grid_B'+str(sel_b)]),
                                                     gridImpact)


        elif sel_b == 'system':

            elecImpact['Co2ImpactSystem'] = np.where(elecTec['PVTotal'] + elecTec['gridTotal'] > 0,
                                                ((elecTec['PVTotal'] * PVImpact)/(elecTec['PVTotal'].sum())
                                                 + elecTec['gridTotal'] * gridImpact) / (
                                                     elecTec['PVTotal'] + elecTec['CHPTotal'] + elecTec['gridTotal']),
                                            gridImpact)

        else:
            print('false definition of sel_b')
            raise NotImplementedError

    return elecImpact


def unit_elecInput(dataDict, buildings):
    """
    Function to read electricity supplied to the heat pumps
    :param :
            dataDict: full result file
            buildings: number of buildings

    :return: elecTec: dataframe with total amount of heat generated per technology
    """

    hpElecInput = pd.DataFrame()

    for b in range(1, buildings+1):
        try:
            ElecList = ["(('electricityInBus__Building" + str(b) + "', '" + str(tec) + "__Building" + str(b) + "'), 'flow')"
                    for tec in ["GWHP", "HP", "ElectricRod", "electricityDemand"]]
            columndict = {ElecList[c]: app_labeldict(labelDict, ElecList[c])[1]
                          for c in range(len(ElecList))}
            hpElecInput = pd.concat([hpElecInput,
                                     dataDict["electricityInBus__Building" + str(b)][ElecList]], axis=1).rename(
                columns=columndict)
        except:
            ElecList = ["(('electricityInBus__Building" + str(b) + "', '" + str(tec) + "__Building" + str(b) + "'), 'flow')"
                        for tec in ["GWHP", "HP", "electricityDemand"]]

            columndict = {ElecList[c]: app_labeldict(labelDict, ElecList[c])[1]
                          for c in range(len(ElecList))}
            hpElecInput = pd.concat([hpElecInput,
                                    dataDict["electricityInBus__Building" + str(b)][ElecList]], axis=1).rename(
                columns=columndict)



        hp_columns = []
        gwhp_columns = []
        ER_columns = []
        BD_columns = []

    for b in range(1, buildings+1):
        gwhp_columns.append('GWHP_B'+str(b))
        hp_columns.append('HP_B'+str(b))
        ER_columns.append('ElectricRod_B'+str(b))
        BD_columns.append('Q_el_B'+str(b))

    hpElecInput['GWHPTotal'] = hpElecInput[gwhp_columns].sum(axis=1)
    hpElecInput['HPTotal'] = hpElecInput[hp_columns].sum(axis=1)
    try:
        hpElecInput['ElectricRodTotal'] = hpElecInput[ER_columns].sum(axis=1)
    except:
        None
    hpElecInput['BuildingDemandElecTotal'] = hpElecInput[BD_columns].sum(axis=1)

    hpElecInput.index = hpElecInput.index

    return hpElecInput


def res_cons(dataDict, buildings):
    """
    Function to determine total amount of ressourses (gas, electricty, etc.) used per technology
    :param :
            dataDict: full result file
            buildings: number of buildings

    :return: resTec: total amount of resources per technology
    """

    gasUse = pd.DataFrame()
    elecUse = pd.DataFrame()

    for b in range(1, buildings+1):

        gasUseList = ["(('naturalGasBus__Building" + str(b) + "', 'CHP__Building" + str(b) + "'), 'flow')",
                     "(('naturalGasBus__Building" + str(b) + "', 'GasBoiler__Building" + str(b) + "'), 'flow')"]
        try:
            elecUseList = ["(('electricityInBus__Building" + str(b) + "', '" + str(tec) + "__Building" + str(b) + "'), 'flow')"
                           for tec in ["GWHP", "HP", "ElectricRod", "solarCollector"]]
            elecAxis = ["GWHP", "HP", "ElectricRod", "solarCollector"]
        except:
            elecUseList = ["(('electricityInBus__Building" + str(b) + "', '" + str(tec) + "__Building" + str(b) + "'), 'flow')"
                           for tec in ["GWHP", "HP", "solarCollector"]]
            elecAxis = ["GWHP", "HP", "solarCollector"]

        gasUse[labelDict["naturalGasBus__Building" + str(b)]] = dataDict["naturalGasBus__Building" + str(b)][gasUseList].sum(axis=0).reset_index(drop=True)
        elecUse[labelDict["electricityInBus__Building" + str(b)]] = dataDict["electricityInBus__Building" + str(b)][elecUseList].sum(axis=0).reset_index(drop=True)

    gasUse = gasUse.set_axis(["CHP", "GasBoiler"], axis=0)
    elecUse = elecUse.set_axis(elecAxis, axis=0)
    gasUse['gasTotal'] = gasUse.sum(axis = 1)
    elecUse['elecTotal'] = elecUse.sum(axis = 1)
    resCons = pd.concat([gasUse, elecUse], axis =1).fillna(0).astype(float)

    return resCons


def cap_technology(dataDict, buildings, energy):
    """
    Function to calculate total capacity per technology selected in the scenario
    :param :
            dataDict: full result file
            buildings: number of buildings
            energy: type str of energy form (heat, elec)

    :return: capTec: capacity per technology
    """

    capTec = pd.DataFrame(index=dataDict['capTransformers__Building'+str(1)].rename(index=labelDict).index)
    capTec["kW in " +str(1)] = pd.DataFrame(dataDict['capTransformers__Building'+str(1)].rename(index=labelDict))

    for b in range(2,buildings+1):
        capTec["kW in " +str(b)] = dataDict['capTransformers__Building'+str(b)].rename(index=labelDict)

    if energy == 'heat':
        try:
            capTec = capTec.loc[['CHP', 'GasBoiler', 'GWHP', 'HP', 'ElectricRod', 'solarConnectBus'], :]
        except:
            capTec = capTec.loc[['CHP', 'GasBoiler', 'GWHP', 'HP', 'solarConnectBus'], :]

    elif energy == 'elec':
        capTec = capTec.loc[['CHP', 'PV'], :]

        for b in range(1,buildings+1):
            capTec.loc['CHP', "kW in " +str(b)] = capTec.loc['CHP', "kW in " +str(b)] * (0.25/0.6) # !!!! should be changed with parametrised efficiency
        capTec = capTec.rename(index={'CHP': 'CHPe'})

    return capTec


def cap_storage(dataDict, buildings):
    """
    Function to calculate total capacity per storage technology selected in the scenario
    :param :
            dataDict: full result file
            buildings: number of buildings
            energy: type str of energy form (heat, elec)

    :return: capSto: capacity per storage
    """

    capSto = pd.DataFrame(index=dataDict['capStorages__Building' + str(1)].rename(index=labelDict).index)
    capSto["kWh in " + str(1)] = pd.DataFrame(dataDict['capStorages__Building' + str(1)].rename(index=labelDict))

    for b in range(2, buildings + 1):
        capSto["kWh in " + str(b)] = dataDict['capStorages__Building' + str(b)].rename(index=labelDict)
    # convert L to kWh
    capSto.loc[capSto.index == 'dhwStorage'] = capSto.loc[capSto.index == 'dhwStorage'] * 4.18 * 45 / 3600
    capSto.loc[capSto.index == 'shStorage'] = capSto.loc[capSto.index == 'shStorage'] * 4.18 * 8 / 3600
    return capSto


def npc_technology(dataDict, inputFileName, buildings, gasCost, elecCost, optMode):
    """
    Function to calculate net preset costs per technology selected in the scenario
    :param :
            dataDict: full result file
            inputFileName: type str of link to input file
            building: number of buildings
            gasCost: cost of gas, timeseries or float
            elecCost: cost of electricity. timeseries or float
            optMode: str group or indiv

    :return: npc per technology
    """

    totalCosts = pd.DataFrame(index = ['naturalGasResource', 'electricityResource', 'Investment', 'Feed-in'])
    for label in [label for label in list(dataDict) if "costs__Building" in label]:
        totalCosts = totalCosts.add(dataDict[label].set_axis(
            ['naturalGasResource', 'electricityResource', 'Investment', 'Feed-in']), fill_value=0)
    print("Net present costs per year are: " + str(float(np.sum(totalCosts))))
    print(totalCosts)

    HeatTec = heat_gen(dataDict, buildings, 'year')
    elecTec = elec_gen(dataDict, buildings, 'year')
    HeatTecHour = heat_gen(dataDict, buildings, 'hour')
    elecTecHour = elec_gen(dataDict, buildings, 'hour')

    if type(elecCost) == str:
        elecPrice = pd.read_csv(elecCost, sep=';', index_col ='timestamp')
        elecPrice.index = pd.to_datetime(elecPrice.index, yearfirst=True)
        elecPrice['cost'] = elecPrice['cost']
    elif type(elecCost) == int or type(elecCost) == float:
        elecPrice = pd.DataFrame(index=rangeToConsider)
        elecPrice["cost"] = elecCost

    levelCosts = {}
    SystemElecPrice = pd.DataFrame()


    ######## Costs input approximation ########
    ## Equivalent Annual Cost for investment costs
    ## input cost as operational costs

    capTecHeat = cap_technology(dataDict, buildings, 'heat')
    capTecElec = cap_technology(dataDict, buildings, 'elec')

    capTec = pd.concat([capTecHeat, capTecElec], axis=0).drop(index='CHPe')
    capTec['total'] = capTec.sum(axis=1)

    sheetTrans = pd.read_excel(inputFileName + str(buildings) + ".xls", sheet_name="transformers")
    sheetSolar = pd.read_excel(inputFileName + str(buildings) + ".xls", sheet_name="solar")

    capTec = capTec.rename(index={'solarConnectBus': 'solarCollector', 'pv': 'PV'})

    CHPData = tec_Data(sheetTrans, 1, 'CHP')
    GasBoilerData = tec_Data(sheetTrans, 1, 'GasBoiler')
    GWHPData = tec_Data(sheetTrans, 1, 'GWHP')
    HPData = tec_Data(sheetTrans, 1, 'HP')
    ERData = tec_Data(sheetTrans, 1, 'ElectricRod')
    SCData = tec_Data(sheetSolar, 1, 'solarCollector')
    PVData = tec_Data(sheetSolar, 1, 'pv')

    elecInput = unit_elecInput(dataDict, buildings)
    capCostsVariable = {'CHP': float(CHPData['invest_cap']),
                        'GasBoiler': float(GasBoilerData['invest_cap']),
                        'GWHP': float(GWHPData['invest_cap']),
                        'HP': float(HPData['invest_cap']),
                        'ElectricRod': float(ERData['invest_cap']),
                        'solarCollector': float(SCData['invest_cap']),
                        'PV': float(PVData['invest_cap'])}
    capCostsBase = {'CHP': float(CHPData['invest_base']),
                    'GasBoiler': float(GasBoilerData['invest_base']),
                    'GWHP': float(GWHPData['invest_base']),
                    'HP': float(HPData['invest_base']),
                    'ElectricRod': float(ERData['invest_base']),
                    'solarCollector': float(SCData['invest_base']),
                    'PV': float(PVData['invest_base'])}
    capCostsMaintenance	 = {'CHP': float(CHPData['maintenance']),
                            'GasBoiler': float(GasBoilerData['maintenance']),
                            'GWHP': float(GWHPData['maintenance']),
                            'HP': float(HPData['maintenance']),
                            'ElectricRod': float(ERData['maintenance']),
                            'solarCollector': float(SCData['maintenance']),
                            'PV': float(PVData['maintenance'])}
    capCostsInstallation = {'CHP': float(CHPData['installation']),
                            'GasBoiler': float(GasBoilerData['installation']),
                            'GWHP': float(GWHPData['installation']),
                            'HP': float(HPData['installation']),
                            'ElectricRod': float(ERData['installation']),
                            'solarCollector': float(SCData['installation']),
                            'PV': float(PVData['installation'])}
    capCostsplanification = {'CHP': float(CHPData['planification']),
                             'GasBoiler': float(GasBoilerData['planification']),
                             'GWHP': float(GWHPData['planification']),
                             'HP': float(HPData['planification']),
                             'ElectricRod': float(ERData['planification']),
                             'solarCollector': float(SCData['planification']),
                             'PV': float(PVData['planification'])}
    gasCostsVariable = gasCost
    tecHeatEff = {'CHP': float(CHPData['efficiency'].split(",")[1]),
                  'GasBoiler': float(GasBoilerData['efficiency'].split(",")[0])}
    timelife = {'CHP': float(CHPData['lifetime']),
                'GasBoiler': float(GasBoilerData['lifetime']),
                'GWHP': float(GWHPData['lifetime']),
                 'HP': float(HPData['lifetime']),
                'ElectricRod': float(ERData['lifetime']),
                'solarCollector': float(SCData['lifetime']),
                'PV': float(PVData['lifetime'])}
    interest = 0.05
    amortisationsValue = {tec: (interest
                                ) / (
                                  1 - (1 + interest)**(-1*timelife[tec])) for tec in capTec.index
                          }

    if "ElectricRodTotal" in HeatTecHour:
        tecInResults = ['CHP', 'GasBoiler', 'GWHP', 'HP', 'ElectricRod', 'solarCollector', 'PV']
    else:
        tecInResults = ['CHP', 'GasBoiler', 'GWHP', 'HP', 'solarCollector', 'PV']

    operationCosts = [
                      gasCostsVariable * sum(HeatTecHour.loc[:, tec + 'Total'])/tecHeatEff[tec] if tec in ['CHP', 'GasBoiler']
                      else sum(elecInput.loc[:,tec + 'Total'] * elecPrice['cost']) if tec in ['GWHP', 'HP']
                      else 0
                      for tec in tecInResults
                      ]

    techInBuilding = {tec: np.count_nonzero(list(capTec.loc[tec, 'kW in 1': 'kW in ' + str(buildings)]))
                     for tec in tecInResults}

    RevenueData = pd.DataFrame({'investmentCosts' :
                                [(techInBuilding[tec] * capCostsBase[tec] + capTec.loc[tec, 'total'] * capCostsVariable[tec]) * amortisationsValue[tec]
                                     if capTec.loc[tec, 'total'] > 0 else 0 for tec in capTec.index ]
                                },
                               index = tecInResults
                                   )

    RevenueData['operationCosts'] = operationCosts
    RevenueData['Maintenance'] = {tec: capCostsMaintenance[tec] * RevenueData.loc[tec, 'investmentCosts']
                                  for tec in tecInResults}


    if elecTec.loc['PV', 'total'] > 0:
        levelCosts['PVElecSystem'] = (RevenueData.loc['PV', 'investmentCosts'] +
                                                RevenueData.loc['PV', 'operationCosts']) / (
                                                elecTec.loc['PV', 'total'])
    elif elecTec.loc['PV', 'total'] == 0:
        levelCosts['PVElecSystem'] = 0

    avoidedCosts = pd.concat([elec_gen(dataDict, buildings, 'hour'), elecPrice['cost']], axis=1)
    avoidedCosts = avoidedCosts.rename(columns={'CHPTotal' : 'CHPElec'})

    avoidedCosts['CHPElecAvoidedCosts'] = avoidedCosts['cost'] * avoidedCosts['CHPElec']

    levelCosts['CHPHeatSystem'] = (RevenueData.loc['CHP', 'investmentCosts'] + RevenueData.loc['CHP', 'operationCosts'] -
                         avoidedCosts['CHPElecAvoidedCosts'].sum()) / (
            HeatTec.loc['CHPsh', 'total'] + HeatTec.loc['CHPdhw', 'total']).sum()

    avoidedCosts = pd.concat([avoidedCosts, HeatTecHour], axis=1)
    avoidedCosts = avoidedCosts.rename(columns={'CHPTotal' : 'CHPHeatTotal'})

    avoidedCosts['CHPHeatAvoidedCosts'] = levelCosts['CHPHeatSystem'] * avoidedCosts['CHPHeatTotal']

    if elecTec.loc['CHP', 'total'] > 0:
        levelCosts['CHPElecSystem'] = (RevenueData.loc['CHP', 'investmentCosts'] + RevenueData.loc['CHP', 'operationCosts'] -
                                   avoidedCosts['CHPHeatAvoidedCosts'].sum()) / (
                                    elecTec.loc['CHP', 'total'])
    else:
        levelCosts['CHPElecSystem'] = 0

    SystemElecPrice['ElecPriceSystem'] = ElecInfluenceBuilding(dataDict, buildings,
                                                               levelCosts['PVElecSystem'],
                                                               levelCosts['CHPElecSystem'],
                                                               elecPrice['cost'],
                                                               'system', 'cost')

    levelCosts['GWHPHeatSystemGrid'] = (RevenueData.loc['GWHP', 'investmentCosts'] + RevenueData.loc['GWHP', 'operationCosts']) / (
            HeatTec.loc['GWHPsh', 'total'] + HeatTec.loc['GWHPdhw', 'total']).sum()

    gwhpOperationCosts = SystemElecPrice['ElecPriceSystem'] * elecInput['GWHPTotal']

    levelCosts['GWHPHeatSystem'] = (RevenueData.loc['GWHP', 'investmentCosts'] + gwhpOperationCosts.sum()) / (
            HeatTec.loc['GWHPsh', 'total'] + HeatTec.loc['GWHPdhw', 'total']).sum()

    levelCosts['GasBoilerHeatSystem'] = (RevenueData.loc['GasBoiler', 'investmentCosts'] +
                               RevenueData.loc['GasBoiler', 'operationCosts']) / (
                HeatTec.loc['GasBoilersh', 'total'] + HeatTec.loc['GasBoilerdhw', 'total']).sum()

    hpOperationCosts = SystemElecPrice['ElecPriceSystem'] * elecInput['HPTotal']

    levelCosts['HPHeatSystem'] = (RevenueData.loc['HP', 'investmentCosts'] + hpOperationCosts.sum()) / (
                HeatTec.loc['HPsh', 'total'] + HeatTec.loc['HPdhw', 'total']).sum()

    levelCosts['HPHeatSystemGrid'] = (RevenueData.loc['HP', 'investmentCosts'] + RevenueData.loc['HP', 'operationCosts']) / (
                HeatTec.loc['HPsh', 'total'] + HeatTec.loc['HPdhw', 'total']).sum()
    try:
        EROperationCosts = SystemElecPrice['ElecPriceSystem'] * elecInput['ElectricRodTotal']

        levelCosts['ElectricRodHeatSystem'] = (RevenueData.loc['ElectricRod', 'investmentCosts'] + EROperationCosts.sum()) / (
                HeatTec.loc['ElectricRodsh', 'total'] + HeatTec.loc['ElectricRoddhw', 'total']).sum()

        levelCosts['ElectricRodHeatSystemGrid'] = (RevenueData.loc['ElectricRod', 'investmentCosts'] + RevenueData.loc[
            'ElectricRod', 'operationCosts']) / (HeatTec.loc['ElectricRodsh', 'total'] + HeatTec.loc['ElectricRoddhw', 'total']).sum()
    except:
        None

    levelCosts['solarCollectorHeatSystem'] = (RevenueData.loc['solarCollector', 'investmentCosts'] +
                                    RevenueData.loc['solarCollector', 'operationCosts']) / (
        HeatTec.loc['solarCollectordhw':,'total']).sum()

    elecLevelCostsSum = pd.DataFrame()
    HeatLevelCostsSum = pd.DataFrame()
    for tec in ['CHP', 'GWHP', 'GasBoiler', 'HP', 'solarCollector']:
        HeatLevelCostsSum['LCSum_' + str(tec)] = HeatTecHour[str(tec) + 'Total'] * levelCosts[str(tec) + 'HeatSystem']
        HeatLevelCostsSum['LCSum_' + str(tec)] = HeatLevelCostsSum['LCSum_' + str(tec)].fillna(0)

    for tec in ['CHP', 'PV']:
        elecLevelCostsSum['LCSum_' + str(tec)] = elecTecHour[str(tec) + 'Total'] * levelCosts[str(tec) + 'ElecSystem']
        elecLevelCostsSum['LCSum_' + str(tec)] = elecLevelCostsSum['LCSum_' + str(tec)].fillna(0)

    HeatLevelCostsSum['LCSystem'] = (HeatLevelCostsSum['LCSum_CHP'].fillna(0)
                                   + HeatLevelCostsSum['LCSum_GWHP'].fillna(0)
                                   + HeatLevelCostsSum['LCSum_GasBoiler'].fillna(0)
                                   + HeatLevelCostsSum['LCSum_HP'].fillna(0)
                                   + HeatLevelCostsSum['LCSum_solarCollector'].fillna(0)
                                     ) / (
                                    HeatTecHour['CHPTotal'].fillna(0) + HeatTecHour['GasBoilerTotal'].fillna(0)
                                    + HeatTecHour['GWHPTotal'].fillna(0) +
                                    HeatTecHour['HPTotal'].fillna(0) + HeatTecHour['solarCollectorTotal']).fillna(0)

    HeatLevelCostsSum['LCSystem'] = HeatLevelCostsSum['LCSystem'].fillna(0)
    elecLevelCostsSum['LCSystem'] = (elecLevelCostsSum['LCSum_CHP'].fillna(0)
                                   + elecLevelCostsSum['LCSum_PV'].fillna(0)
                                    ) / (
                                    elecTecHour['PVTotal'].fillna(0) + elecTecHour['CHPTotal'].fillna(0))

    elecLevelCostsSum['LCSystem'] = elecLevelCostsSum['LCSystem'].fillna(0)

    return HeatLevelCostsSum, elecLevelCostsSum, levelCosts, SystemElecPrice


def levelizedCostsPlot(dataDict, buildings, gasCost, elecCost, iter):
    """
    Function to plot the levelized costs
    :param :
            dataDict: full result file
            building: number of buildings
            gasCost: cost of gas, timeseries or float
            elecCost: cost of electricity. timeseries or float
            iter: optimization iteration

    :return:
    """

    HeatLevelCostsSum, elecLevelCostsSum, levelCosts, SystemElecPrice = npc_technology(dataDict, inputFileName,
                                                                                       buildings, gasCost, elecCost,
                                                                                       optMode)

    ### plot
    fig = plt.figure()
    HeatLevelCostsSum['LCSystem'].plot()
    plt.savefig(outputFileName + "Optimization" + str(iter) + '/heatLevelizedCostsSystem.png', bbox_inches='tight')
    fig.clf()
    plt.close()

    fig = plt.figure()
    elecLevelCostsSum['LCSystem'].plot()
    plt.savefig(outputFileName + "Optimization" + str(iter) + '/elecLevelizedCostsSystem.png', bbox_inches='tight')
    fig.clf()
    plt.close()


def selfsuffisant(dataDict, buildings, outputFileName, selected_days, timeStep, iter, iterRange, results):
    """
    Function to plot electricty supply from self production and the grid
    :param :
            dataDict: full result file
            outputFileName: type str of link to output file
            building: number of buildings
            selected_days: str of datetime for selected day
            timeStep: str: hour, month, year
            iter: optimization iteration
            iterRange: list of ploting order of optimization iterations
            results: dict where iteration results are stored

    :return: bar chart
    """

    if timeStep == 'hour':
        elecGenTec = elec_gen(dataDict, buildings, "hour")
        elecGenTec = pd.concat([elecGenTec, elec_sell(dataDict, buildings, 'hour')*(-1)], axis=1)
        elecGenTecToPlot = elecGenTec[['CHPTotal', 'PVTotal', 'gridTotal', 'excessTotal']]/1000

        fig = plt.figure()
        elecGenTecToPlot.rename(columns=labelDict).plot(
            kind='bar', stacked=True, title="Source of electricity supply", xlabel='Month', ylabel='Electricity per source in MWh')

        plt.legend(loc=(1.04, 0)) #, ['CHP', 'PV', 'Grid', 'Sold']
        ax = plt.gca()
        ax.axes.xaxis.set_ticklabels([])
        plt.tight_layout()
        plt.savefig(outputFileName + "Optimization" + str(iter) + '/hourly_selfsuffisant.png', bbox_inches='tight')
        fig.clf()
        plt.close()

        elecGenTec['ratioSS'] = (elecGenTec['CHPTotal'] + elecGenTec['PVTotal'] + elecGenTec['excessTotal']) / (
                elecGenTec['CHPTotal'] + elecGenTec['PVTotal'] + elecGenTec['excessTotal'] + elecGenTec['gridTotal'])

        fig = plt.figure()
        plt.plot(elecGenTec['ratioSS'].index, elecGenTec['ratioSS'])
        plt.xlabel('')
        plt.ylabel('Ratio of electricity production to consumption')
        plt.tight_layout()
        plt.savefig(outputFileName + "Optimization" + str(iter) + '/hourly_ratio_selfsuffisant.png', bbox_inches='tight')

        fig.clf()
        plt.close()

        for days in selected_days:
            elecGenTecToPlot = elecGenTec[['CHPTotal', 'PVTotal', 'gridTotal', 'excessTotal']]/1000

            fig = plt.figure()
            elecGenTecToPlot.rename(columns=labelDict).loc[days: days + pd.offsets.DateOffset(hour=23)].plot(
                kind='bar', stacked=True, title="Source of electricity supply", xlabel='Month', ylabel='Electricity per source in MWh')

            plt.legend(loc=(1.04, 0)) # , ['CHP', 'PV', 'Grid', 'Sold']
            ax = plt.gca()
            ax.axes.xaxis.set_ticklabels([])
            plt.tight_layout()
            plt.savefig(outputFileName + "Optimization" + str(iter) + '/hourly_selfsuffisant' + str(days)[0:10] + '.png', bbox_inches='tight')
            fig.clf()
            plt.close()

            fig = plt.figure()
            plt.plot(elecGenTec['ratioSS'].loc[days: days + pd.offsets.DateOffset(hour=23)].index,
                     elecGenTec['ratioSS'].loc[days: days + pd.offsets.DateOffset(hour=23)])
            plt.xlabel('')
            plt.ylabel('Ratio of electricity production to consumption')
            plt.tight_layout()
            plt.savefig(outputFileName + "Optimization" + str(iter) + '/hourly_ratio_selfsuffisant' + str(days)[0:10] + '.png', bbox_inches='tight')
            fig.clf()
            plt.close()

    elif timeStep == 'month':
        elecGenTec = elec_gen(dataDict, buildings, "month")
        elecGenTec = pd.concat([elecGenTec, elec_sell(dataDict, buildings, 'month')*(-1)], axis=1)
        elecGenTecToPlot = elecGenTec[['CHPTotal', 'PVTotal', 'gridTotal', 'excessTotal']]/1000

        fig = plt.figure()
        elecGenTecToPlot.rename(columns=labelDict).plot(
            kind='bar', stacked=True, title="Source of electricity supply", xlabel='Month', ylabel='Electricity per source in MWh')
        plt.legend(loc=(1.04, 0))
        plt.tight_layout()
        plt.savefig(outputFileName + "Optimization" + str(iter) + '/selfsuffisant.png', bbox_inches='tight')
        fig.clf()
        plt.close()

        elecGenTec['ratioSS'] = (elecGenTec['CHPTotal'] + elecGenTec['PVTotal'] + elecGenTec['excessTotal']) / (
                elecGenTec['CHPTotal'] + elecGenTec['PVTotal'] + elecGenTec['excessTotal'] + elecGenTec['gridTotal'])

        fig = plt.figure()
        plt.plot(elecGenTec['ratioSS'].index, elecGenTec['ratioSS'])
        plt.title("Ratio of electricity production (excl. injected to grid) to consumption")
        plt.xlabel('')
        plt.ylabel('Ratio')
        plt.tight_layout()
        plt.savefig(outputFileName + "Optimization" + str(iter) + '/ratio_selfsuffisant.png', bbox_inches='tight')
        fig.clf()
        plt.close()


    elif timeStep == 'year':
        if iter == iterRange[0]:
            results['totalselfsuffisant'] = {}
            results['ratioselfsuffisant'] = {}
        else:
            None

        results['totalselfsuffisant'][iter] = elec_gen(dataDict, buildings, "year")['total']
        results['totalselfsuffisant'][iter]['Excess'] = (elec_sell(dataDict, buildings, 'year')*(-1)).sum()

        if iter == iterRange[-1]:
            fig = plt.figure()
            (pd.DataFrame(results['totalselfsuffisant'])/1000).T.plot(kind="bar", stacked=True)
            plt.title("Electricity source of each optimization")
            plt.xlabel('Optimization iteration')
            plt.ylabel('Electricity supply per technology in MWh')
            plt.legend(loc=(1.04, 0) )
            plt.tight_layout()
            plt.savefig(outputFileName + 'selfsuffisant_iterative.png', bbox_inches='tight')
            fig.clf()
            plt.close()

        else:
            None


        results['ratioselfsuffisant'][iter] = (results['totalselfsuffisant'][iter]['CHP']
                                               + results['totalselfsuffisant'][iter]['PV']
                                               + results['totalselfsuffisant'][iter]['Excess']
                                               ) / (results['totalselfsuffisant'][iter]['CHP']
                                                  + results['totalselfsuffisant'][iter]['PV']
                                                  + results['totalselfsuffisant'][iter]['Excess']
                                                  + results['totalselfsuffisant'][iter]['Grid'])

        if iter == iterRange[-1]:
            fig = plt.figure()
            plt.plot([str(x) for x in list(results['ratioselfsuffisant'].keys())], list(results['ratioselfsuffisant'].values()))
            plt.title("Ratio of electricity production (excl. injected to grid) to consumption")
            plt.xlabel('Optimization iteration')
            plt.ylabel('Ratio')
            plt.legend(loc=(1.04, 0) )
            plt.tight_layout()
            plt.savefig(outputFileName + 'ratioselfsuffisant_iterative.png', bbox_inches='tight')
            fig.clf()
            plt.close()

    return results


def heat_distr(dataDict, buildings, iter, outputFileName, timeStep):
    """
    Function to visualize the heat distribution for each technology and selected time step
    :param :
            dataDict: full result file
            outputFileName: type str of link to output file
            building: number of buildings
            timeStep: str: hour, month, year
    :return: plot
    """

    if timeStep == 'year':
        heatTec = heat_gen(dataDict, buildings, 'year')
        totalHeatTec = {'chp': int(heatTec.loc['CHPsh', 'total']) + int(heatTec.loc['CHPdhw', 'total']),
                        'GWHP': int(heatTec.loc['GWHPsh', 'total']) + int(heatTec.loc['GWHPdhw', 'total']),
                        'GasBoiler': int(heatTec.loc['GasBoilersh', 'total'])+ int(heatTec.loc['GasBoilerdhw', 'total']),
                        'HP': int(heatTec.loc['HPsh', 'total']) + int(heatTec.loc['HPdhw', 'total']),
                        'solarCollector': int(heatTec.loc['solarCollectordhw', 'total'])}
        try:
            totalHeatTec['ElectricRod'] = int(heatTec.loc['ElectricRodsh', 'total']) + int(heatTec.loc['ElectricRoddhw', 'total'])
        except:
            None

        fig = plt.figure()
        plt.bar(range(len(totalHeatTec)), [i//1000 for i in list(totalHeatTec.values())], align='center')
        plt.xticks(range(len(totalHeatTec)), list(totalHeatTec.keys()))
        plt.ylabel('Heat generated in MWh')
        plt.savefig(outputFileName + "Optimization" + str(iter) + '/heatDistrYear.png', bbox_inches='tight')
        fig.clf()
        plt.close()


        return totalHeatTec

    elif timeStep == 'month':
        heatTec = heat_gen(dataDict, buildings, 'month')#

        ### divided into sh and dhw
        try:
            heatTecToPlot = heatTec[['CHPshTotal', 'GWHPshTotal', 'GasBoilershTotal', 'HPshTotal', 'ElectricRodshTotal']]/1000
        except:
            heatTecToPlot = heatTec[['CHPshTotal', 'GWHPshTotal', 'GasBoilershTotal', 'HPshTotal']]/1000
        fig = plt.figure()
        heatTecToPlot.rename(columns=labelDict).plot(
            kind='bar', stacked=True, title="Source of heat supply",
            xlabel='Month', ylabel='Space heating per source in MWh')

        plt.legend(loc=(1.04, 0) )
        plt.savefig(outputFileName + "Optimization" + str(iter) + '/shDistrMonth.png', bbox_inches='tight')
        fig.clf()
        plt.close()
        try:
            heatTecToPlot = heatTec[['CHPdhwTotal', 'GWHPdhwTotal', 'GasBoilerdhwTotal', 'HPdhwTotal', 'ElectricRoddhwTotal']]/1000
        except:
            heatTecToPlot = heatTec[['CHPdhwTotal', 'GWHPdhwTotal', 'GasBoilerdhwTotal', 'HPdhwTotal']]/1000
        fig = plt.figure()
        heatTecToPlot.rename(columns=labelDict).plot(
            kind='bar', stacked=True, title="Source of heat supply", xlabel='Month', ylabel='Dhw per source in MWh')
        plt.legend(loc=(1.04, 0))
        plt.savefig(outputFileName + "Optimization" + str(iter) + '/dhwDistrMonth.png', bbox_inches='tight')
        fig.clf()
        plt.close()

        try:
            heatTecToPlot = heatTec[['CHPTotal', 'GWHPTotal', 'GasBoilerTotal', 'HPTotal', 'ElectricRodTotal']]/1000
        except:
            heatTecToPlot = heatTec[['CHPTotal', 'GWHPTotal', 'GasBoilerTotal', 'HPTotal']]/1000

        fig = plt.figure()
        heatTecToPlot.rename(columns=labelDict).plot(
            kind='bar', stacked=True, title="Source of heat supply", xlabel='Month', ylabel='Heat per source in MWh')
        plt.legend(loc=(1.04, 0) )
        plt.savefig(outputFileName + "Optimization" + str(iter) + '/heatDistrMonth.png', bbox_inches='tight')
        fig.clf()
        plt.close()

        return heatTec


def full_load_hour(dataDict, buildings, iter, outputFileName):
    """
    Function to plot part and full load time
    :param :
            dataDict: full result file
            outputFileName: type str of link to output file
            building: number of buildings
            iter: optimization iteration
    :return: bar chart
    """

    shHeatProd = pd.DataFrame()
    dhwHeatProd = pd.DataFrame()

    for b in range(1,buildings+1):
        try:
            shGenList = ["(('"+ str(tec) + "__Building" + str(b) + "', 'shSourceBus__Building" + str(b) + "'), 'flow')"
                         for tec in ["CHP", "GWHP", "GasBoiler", "HP", "ElectricRod"]]
            shColumns = ["CHPsh" + str(b), "GWHPsh" + str(b), "GasBoilersh" + str(b), "HPsh" + str(b), "ElectricRodsh" + str(b)]

            dhwGenList = ["(('"+ str(tec) + "__Building" + str(b) + "', 'dhwStorageBus__Building" + str(b) + "'), 'flow')"
                        for tec in ["CHP", "GWHP", "GasBoiler", "HP", "ElectricRod", "solarCollector"]]

            dhwColumns = ["CHPdhw" + str(b), "GWHPdhw" + str(b), "GasBoilerdhw" + str(b), "HPdhw" + str(b),
                          "ElectricRoddhw" + str(b), "solarCollectordhw" + str(b)]
            shBuilding = dataDict["shSourceBus__Building" + str(b)][shGenList]

            shBuilding.columns = shColumns
            shHeatProd = pd.concat([shHeatProd, shBuilding], axis=1)

            dhwBuilding = dataDict["domesticHotWaterBus__Building" + str(b)][dhwGenList]

            dhwBuilding.columns = dhwColumns
            dhwHeatProd = pd.concat([dhwHeatProd, dhwBuilding], axis=1)

            tecInResults = ['CHP', 'GWHP', 'GasBoiler', 'HP', 'ElectricRod']
            ElectecInResults = ["GWHP", "HP", "ElectricRod"]
        except:
            shGenList = ["(('"+ str(tec) + "__Building" + str(b) + "', 'shSourceBus__Building" + str(b) + "'), 'flow')"
                         for tec in ["CHP", "GWHP", "GasBoiler", "HP"]]
            shColumns = ["CHPsh" + str(b), "GWHPsh" + str(b), "GasBoilersh" + str(b), "HPsh" + str(b)]

            dhwGenList = ["(('"+ str(tec) + "__Building" + str(b) + "', 'dhwStorageBus__Building" + str(b) + "'), 'flow')"
                          for tec in ["CHP", "GWHP", "GasBoiler", "HP", "solarCollector"]]
            dhwColumns = ["CHPdhw" + str(b), "GWHPdhw" + str(b), "GasBoilerdhw" + str(b), "HPdhw" + str(b),
                          "solarCollectordhw" + str(b)]

            shBuilding = dataDict["shSourceBus__Building" + str(b)][shGenList]

            shBuilding.columns = shColumns
            shHeatProd = pd.concat([shHeatProd, shBuilding], axis=1)

            dhwBuilding = dataDict["domesticHotWaterBus__Building" + str(b)][dhwGenList]

            dhwBuilding.columns = dhwColumns
            dhwHeatProd = pd.concat([dhwHeatProd, dhwBuilding], axis=1)

            tecInResults = ['CHP', 'GWHP', 'GasBoiler', 'HP']
            ElectecInResults = ["GWHP", "HP"]

    heatProd = pd.concat([shHeatProd, dhwHeatProd], axis=1)

    capTec = cap_technology(dataDict, buildings, 'heat')
    capTec.loc['HP',:] = capTec.loc['HP',:]/3.5
    capTec.loc['GWHP',:] = capTec.loc['GWHP',:]/4.65

    capTec['total'] = capTec.sum(axis=1)

    if capTec.loc["solarConnectBus", 'total'] > 0:
        fig = plt.figure()
        plt.plot(range(len(heatProd.index)), heatProd[[str("solarCollectordhw") + str(b+1) for b in range(buildings)]].sum(axis=1)
                 .sort_values(ascending=False)/capTec.loc["solarConnectBus", 'total'], label = "Total")

        for b_ in range(buildings):

            plt.plot(range(len(heatProd.index)), heatProd[str("solarCollectordhw") + str(b_ +1)].sort_values(ascending=False)/capTec
                 .loc["solarConnectBus", 'kW in ' + str(b_+1)],
                 label="Building " + str(b_+1) +" (" + str(int(capTec.loc["solarConnectBus", 'kW in ' + str(b_+1)])) + "kW)")

        plt.title("solar Collector for dhw")
        plt.xlabel('Operation hour per year')
        plt.ylabel('Part-load operation in % of installed heat capacity')
        plt.legend(loc=(1.04, 0) )
        plt.savefig(outputFileName + "Optimization" + str(iter) + '/operation_solarCollector_dhw.png', bbox_inches='tight')

        fig.clf()
        plt.close()

    else:
        None

    for tec in tecInResults:
        for h in ['sh', 'dhw']:
            if capTec.loc[str(tec), 'total']>0:
                fig = plt.figure()
                plt.plot(range(len(heatProd.index)), heatProd[[str(tec) + str(h) + str(b+1) for b in range(buildings)]].sum(axis=1)
                         .sort_values(ascending=False)/capTec.loc[str(tec), 'total'], label = "Total")

                for b_ in range(buildings):
                    plt.plot(range(len(heatProd.index)),heatProd[str(tec) + str(h) + str(b_ +1)].sort_values(ascending=False)/capTec
                             .loc[str(tec), 'kW in ' + str(b_+1)],
                             label="Building " + str(b_+1) +" (" + str(int(capTec.loc[str(tec), 'kW in ' + str(b_+1)])) + "kW)")

                plt.title(str(tec) + " for " + str(h))
                plt.xlabel('Operation hour per year')
                plt.ylabel('Part-load operation in % of installed heat capacity')
                plt.legend(loc=(1.04, 0) )
                plt.savefig(outputFileName + "Optimization" + str(iter) + '/operation_' + str(tec) + str(h) + '.png', bbox_inches='tight')
                fig.clf()
                plt.close()

            else:
                None

    hpElecInput = unit_elecInput(dataDict, buildings)
    heatProd = pd.concat([heatProd, hpElecInput], axis=1)


    for tec in ['CHP', 'GasBoiler']:
        if capTec.loc[str(tec), 'total']>0:
            fig = plt.figure()
            plt.plot(range(len(heatProd.index)), heatProd[[str(tec) + str(h) + str(b+1) for h in ['sh', 'dhw']
                                                               for b in range(buildings)]].sum(axis=1)
                             .sort_values(ascending=False)/capTec.loc[str(tec), 'total'], label="Total")

            for b_ in range(buildings):
                    plt.plot(range(len(heatProd.index)),heatProd[[str(tec) + str(h) + str(b_ +1) for h in ['sh', 'dhw']
                                                      ]].sum(axis=1).sort_values(ascending=False)/capTec
                        .loc[str(tec), 'kW in ' + str(b_+1)],
                         label="Building " + str(b_+1) +" (" + str(int(capTec.loc[str(tec), 'kW in ' + str(b_+1)])) + "kW)")

            plt.xlabel('Operation hour per year')
            plt.ylabel('Part-load operation in % of installed capacity')
            plt.title(str(tec) + " for space heating and dhw combined")
            plt.legend(loc=(1.04, 0) )
            plt.savefig(outputFileName + "Optimization" + str(iter) + '/operation_' + str(tec) + '.png', bbox_inches='tight')
            fig.clf()
            plt.close()

        else:
            None

    for tec in ElectecInResults:
        if capTec.loc[str(tec), 'total']>0:
            fig = plt.figure()
            plt.plot(range(len(heatProd.index)), heatProd[[str(tec) + "_B" + str(b+1)
                                                               for b in range(buildings)]].sum(axis=1)
                         .sort_values(ascending=False)/capTec.loc[str(tec), 'total'], label="Total")
            for b_ in range(buildings):
                plt.plot(range(len(heatProd.index)),heatProd[[str(tec) + "_B" + str(b_ +1)
                                                                  ]].sum(axis=1).sort_values(ascending=False)/capTec
                             .loc[str(tec), 'kW in ' + str(b_+1)],
                             label="Building " + str(b_+1) +" (" + str(int(capTec.loc[str(tec), 'kW in ' + str(b_+1)])) + "kW)")

            plt.xlabel('Operation hour per year')
            plt.ylabel('Part-load operation in % of installed capacity (electric for hp)')
            plt.title(str(tec) + " for space heating and dhw combined")
            plt.legend(loc=(1.04, 0) )
            plt.savefig(outputFileName + "Optimization" + str(iter) + '/operation_' + str(tec) + '.png', bbox_inches='tight')
            fig.clf()
            plt.close()

        else:
            None

    fig = plt.figure()
    plt.plot(range(len(heatProd.index)), heatProd[["solarCollectordhw" + str(b+1)
                                               for b in range(buildings)]].sum(axis=1).sort_values(
        ascending=False)/capTec.loc["solarConnectBus", 'total'], label="solar collector (" + str(int(capTec.loc["solarConnectBus", 'total'])) + "kW)")


    for tec in ['CHP',  'GasBoiler']:

        plt.plot(range(len(heatProd.index)), heatProd[[str(tec) + str(h) + str(b+1) for h in ['sh', 'dhw']
                                                           for b in range(buildings)]].sum(axis=1)
                     .sort_values(ascending=False)/capTec.loc[str(tec), 'total'], label=tec + " (" + str(int(capTec.loc[tec,'total'])) + "kW)")

    for tec in ElectecInResults:
        plt.plot(range(len(heatProd.index)), heatProd[[str(tec) + '_B' + str(b+1)
                                                           for b in range(buildings)]].sum(axis=1)
                     .sort_values(ascending=False)/capTec.loc[str(tec), 'total'], label=tec + " (" + str(int(capTec.loc[tec,'total'])) + "kW)")

    plt.xlabel('Operation hour per year')
    plt.ylabel('Part-load operation in % of installed capacity (electric for hp)')
    plt.title("Heating technologies for space heating and dhw combined")
    plt.legend(loc=(1.04, 0) )
    plt.savefig(outputFileName + "Optimization" + str(iter) + '/operation.png', bbox_inches='tight')
    fig.clf()
    plt.close()


def stacked_full_load(dataDict, buildings, iter, outputFileName):
    """
    create plot chart of the heat generation in the system
    :param :
            dataDict: full result file
            outputFileName: type str of link to output file
            building: number of buildings
            iter: optimization iteration

    :return: plot
    """
    heatGen = heat_gen(dataDict, buildings, 'hour')
    try:
        heatGen['heatTotal'] = heatGen['CHPTotal'] + heatGen['GWHPTotal'] + heatGen['GasBoilerTotal'] + \
                               heatGen['HPTotal'] + heatGen['ElectricRod'] + heatGen['solarCollectorTotal']
    except:
        heatGen['heatTotal'] = heatGen['CHPTotal'] + heatGen['GWHPTotal'] + heatGen['GasBoilerTotal'] + \
                               heatGen['HPTotal'] + heatGen['solarCollectorTotal']

    fig = plt.figure()
    plt.plot(range(len(heatGen.index)), heatGen.sum(axis=1).sort_values(ascending=False))
    plt.xlabel('Hours')
    plt.ylabel('Combined heat production in kWh')
    plt.savefig(outputFileName + "Optimization" + str(iter) + '/demand_monoton.png', bbox_inches='tight')
    fig.clf()
    plt.close()

    fig = plt.figure()
    try:
        plt.stackplot(range(len(heatGen)),
                      list(heatGen.sort_values(by='heatTotal', ascending=False)['CHPTotal']/1000),
                      list(heatGen.sort_values(by='heatTotal', ascending=False)['GWHPTotal']/1000),
                      list(heatGen.sort_values(by='heatTotal', ascending=False)['HPTotal']/1000),
                      list(heatGen.sort_values(by='heatTotal', ascending=False)['ElectricRodTotal']/1000),
                      list(heatGen.sort_values(by='heatTotal', ascending=False)['GasBoilerTotal']/1000),
                      list(heatGen.sort_values(by='heatTotal', ascending=False)['solarCollectorTotal']/1000),
                      baseline ='zero',
                      colors =['blue', 'orange', 'green', 'brown', 'red', 'yellow'])
        plt.legend(['CHP', 'GWHP', 'HP', 'ElectricRod', 'Gas Boiler', 'solar Collector'])
    except:
        plt.stackplot(range(len(heatGen)),
                      list(heatGen.sort_values(by='heatTotal', ascending=False)['CHPTotal']/1000),
                      list(heatGen.sort_values(by='heatTotal', ascending=False)['GWHPTotal']/1000),
                      list(heatGen.sort_values(by='heatTotal', ascending=False)['HPTotal']/1000),
                      list(heatGen.sort_values(by='heatTotal', ascending=False)['GasBoilerTotal']/1000),
                      list(heatGen.sort_values(by='heatTotal', ascending=False)['solarCollectorTotal']/1000),
                      baseline ='zero',
                      colors =['blue', 'orange', 'green', 'red', 'yellow'])
        plt.legend(['CHP', 'GWHP', 'HP', 'Gas Boiler', 'solar Collector'])

    plt.legend(loc=(1.04, 0) )

    plt.title('Heat generation operation')
    plt.xlabel('Hour of the year')
    plt.ylabel('Heat generation in MWh')
    plt.savefig(outputFileName + "Optimization" + str(iter) + '/monoton_tec.png', bbox_inches='tight')
    fig.clf()
    plt.close()


def installed_capacity(dataDict, buildings, outputFileName, iter, iterRange, results, system, xlabels=None):
    """
    Function to plot installed capacity of all optimization iterations
    :param :
            dataDict: full result file
            outputFileName: type str of link to output file
            building: number of buildings
            iter: optimization iteration
            iterRange: list of plotting order of optimization iterations
            results: dict where iteration results are stored
            system: how to plot the capacity: type str building or grouped
            xlabels: labels to use on the x axis, iterable or None, if None, xlabels = iter
    :return: plot
    """

    # set xlabels to iterRange if None:
    if not xlabels:
        xlabels = iterRange


    if system == 'building':
        capTecHeat = cap_technology(dataDict, buildings, 'heat')
        capTecHeat = capTecHeat.rename(index={"solarConnectBus":"solarCollector"})

        fig = plt.figure()
        capTecHeat.T.plot(kind='bar', stacked=True)
        plt.title("Installed capacity in each building")
        plt.xlabel('Building')
        plt.ylabel('Capacity in kW')
        plt.legend(loc=(1.04, 0))
        plt.savefig(outputFileName + "Optimization" + str(iter) + '/installedCapHeat.png', bbox_inches='tight')
        fig.clf()
        plt.close()

        capTecElec = cap_technology(dataDict, buildings, 'elec')

        fig = plt.figure()
        capTecElec.T.plot(kind='bar', stacked=True)
        plt.title("Installed capacity in each building")
        plt.xlabel('Building')
        plt.ylabel('Capacity in kW')
        plt.legend(loc=(1.04, 0))
        plt.savefig(outputFileName + "Optimization" + str(iter) + '/installedCapElec.png', bbox_inches='tight')
        fig.clf()
        plt.close()

        capSto = cap_storage(dataDict, buildings)

        fig = plt.figure()
        capSto.T.plot(kind='bar', stacked=True)
        plt.title("Installed capacity in each building")
        plt.xticks(xlabels)
        plt.xlabel('Building')
        plt.ylabel('Capacity in kWh')
        plt.legend(loc=(1.04, 0))
        plt.savefig(outputFileName + "Optimization" + str(iter) + '/installedCapSto.png', bbox_inches='tight',dpi=600)
        fig.clf()
        plt.close()

    elif system == 'grouped':
        if iter == iterRange[0]:
            results['totalCap', 'elec'] = {}
            results['totalCap', 'heat'] = {}
            results['totalCap', 'sto'] = {}
        else:
            None

        results['totalCap', 'elec'][iter] = cap_technology(dataDict, buildings, 'elec').sum(axis=1)
        results['totalCap', 'heat'][iter] = cap_technology(dataDict, buildings, 'heat').sum(axis=1)
        results['totalCap', 'heat'][iter] = results['totalCap', 'heat'][iter].rename(
            index={"solarConnectBus": "solarCollector"})
        results['totalCap', 'sto'][iter] = cap_storage(dataDict, buildings).sum(axis=1)

        if iter == iterRange[-1]:
            fig = plt.figure()
            pd.DataFrame(results['totalCap', 'heat']).T.plot(kind="bar", stacked=True)
            xticks = plt.xticks()
            plt.xticks(xticks[0], xlabels)
            plt.title("Installed capacity of each optimization")
            plt.xlabel('Optimization iteration')
            plt.ylabel('Installed capacity per technology in kW')
            plt.legend(loc=(1.04, 0))
            plt.savefig(outputFileName + 'installedCapHeat_iterative.png', bbox_inches='tight')
            fig.clf()
            plt.close()

            fig = plt.figure()
            pd.DataFrame(results['totalCap','elec']).T.plot(kind="bar", stacked=True)
            xticks = plt.xticks()
            plt.xticks(xticks[0], xlabels)
            plt.title("Installed capacity of each optimization")
            plt.xlabel('Optimization iteration')
            plt.ylabel('Installed capacity per technology in kW')
            plt.legend(loc=(1.04, 0))
            plt.savefig(outputFileName + 'installedCapElec_iterative.png', bbox_inches='tight')
            fig.clf()
            plt.close()

            fig = plt.figure()
            pd.DataFrame(results['totalCap', 'sto']).T.plot(kind="bar", stacked=True)
            xticks = plt.xticks()
            plt.xticks(xticks[0], xlabels)
            plt.title("Installed capacity of each optimization")
            plt.xlabel('Optimization iteration')
            plt.ylabel('Installed capacity per technology in kWh')
            plt.legend(loc=(1.04, 0))
            plt.savefig(outputFileName + 'installedCapSto_iterative.png', bbox_inches='tight')
            fig.clf()
            plt.close()

    else:
        None

    return results


def iter_heat(dataDict, buildings, outputFileName, iter, iterRange, results):
    """
    Plot heat generation of all optimsation iterations
    :param :
            dataDict: full result file
            outputFileName: type str of link to output file
            building: number of buildings
            iter: optimization iteration
            iterRange: list of ploting order of optimization iterations
            results: dict where iteration results are stored

    :return: plot
    """

    if iter == iterRange[0]:
        results['totalHeatTec'] = {}
    else:
        None

    results['totalHeatTec'][iter] = heat_distr(dataDict, buildings, iter, outputFileName, 'year')

    if iter == iterRange[-1]:
        fig = plt.figure()
        (pd.DataFrame(results['totalHeatTec'])/1000).T.plot(kind="bar", stacked=True)
        plt.title("Heat production of each optimization")
        plt.xlabel('Optimization iteration')
        plt.ylabel('Heat generation per technology in MWh')
        plt.legend(loc=(1.04, 0) )
        plt.savefig(outputFileName + 'heatproduction_iterative.png', bbox_inches='tight')
        fig.clf()
        plt.close()

    else:
        None

    return results


def co2_balance(dataDict, inputFileName, buildings, selected_days, gasEmission, elecEmission, rangeToConsider):
    """
    Function to GHG emission impact for different technologies
    :param :
            dataDict: full result file
            inputFileName: type str of link to input file
            building: number of buildings
            selected_days: str of datetime for selected day
            gasEmission: time serie or float with GHG emission impact of gas
            elecEmission: time serie or float with GHG emission impact of grid electricity
            rangeToConsider: time serie

    :return: adding emission impact to heatGenTec and elecGenTec
    """

    print('determine co2 balance')
    heatGenTec = heat_gen(dataDict, buildings, 'hour')
    elecGenTec = elec_gen(dataDict, buildings, 'hour')
    capTecHeat = cap_technology(dataDict, buildings, 'heat')
    capTecElec = cap_technology(dataDict, buildings, 'elec')
    capTec = pd.concat([capTecHeat, capTecElec], axis=0)

    sheetTrans = pd.read_excel(inputFileName + str(buildings) + ".xls", sheet_name="transformers")
    sheetSolar = pd.read_excel(inputFileName + str(buildings) + ".xls", sheet_name="solar")
    sheetSources = pd.read_excel(inputFileName + str(buildings) + ".xls", sheet_name="commodity_sources")


    if type(elecEmission) == str:
        elecImpact = pd.read_csv(elecEmission, sep=';', usecols = ["timestamp", "impact"], index_col ='timestamp')
        elecImpact.index = pd.to_datetime(elecImpact.index, yearfirst=True)

    elif type(elecEmission) == int or type(elecEmission) == float:
        elecImpact = pd.DataFrame(index=pd.date_range(rangeToConsider))
        elecImpact['electricity_impact'] = elecImpact

    elecGenTec['impact'] = list(elecImpact['impact'])

    hpElecInput = unit_elecInput(dataDict, buildings)

    capImpact = {}
    for b in range(1,buildings+1):
        ### multiplying co2 impact with generation
        elecGenTec['GridImpact' + str(b)] = elecGenTec['impact'] * elecGenTec['Grid_B' + str(b)]

        ## CHP co2 impact based on electricity part (gas + cap)

        gasImpactBuilding = gasEmission

        CHPData = tec_Data(sheetTrans, str(b), 'CHP')
        CHPElecShare = 1/(float(CHPData['efficiency'].split(",")[1]) + float(CHPData['efficiency'].split(",")[0]))\
                         * float(CHPData['efficiency'].split(",")[0])
        capImpact[b,"CHPElecCo2"] = capTec.loc["CHP", "kW in " + str(b)] * CHPData['impact_cap'] \
                                    * CHPElecShare / (CHPData['lifetime'].mean())

        if sum(elecGenTec['CHP_B' + str(b)]) > 0:
            elecGenTec['CHPGasImpact' + str(b)] = gasImpactBuilding * CHPElecShare \
                                                  * elecGenTec['CHP_B' + str(b)] / float(CHPData['efficiency'].split(",")[0])
            elecGenTec['CHPImpact' + str(b)] = elecGenTec['CHPGasImpact' + str(b)] + capImpact[b,"CHPElecCo2"]*elecGenTec['CHP_B' + str(b)]/(elecGenTec['CHP_B' + str(b)]).sum()

        else:
            elecGenTec['CHPGasImpact' + str(b)] = 0
            elecGenTec['CHPImpact' + str(b)] = 0

        ## PV
        PVData = tec_Data(sheetSolar, str(b), 'pv')
        capImpact[b,"PVElecCo2"] = capTec.loc["PV", "kW in " + str(b)] * PVData['impact_cap'] / (PVData['lifetime'].mean())
        elecGenTec['PVImpact' + str(b)] = capImpact[b,"PVElecCo2"] * elecGenTec['PV_B' + str(b)] / (elecGenTec['PV_B' + str(b)]).sum()

        ### co2 impact for heat
        ## chp
        capImpact[b,"CHPHeatCo2"] = capTec.loc["CHP", "kW in " + str(b)] * CHPData['impact_cap']\
                                    * (1 - CHPElecShare) / (CHPData['lifetime'].mean())

        if sum(heatGenTec['CHP' + str(b)]) > 0:
            heatGenTec['CHPGasImpact' + str(b)] = gasImpactBuilding * (1 - CHPElecShare) \
                                                  * heatGenTec['CHP' + str(b)] / float(CHPData['efficiency'].split(",")[1])
            heatGenTec['CHPImpact' + str(b)] = heatGenTec['CHPGasImpact' + str(b)] + capImpact[b,"CHPHeatCo2"]*heatGenTec['CHP' + str(b)]/(heatGenTec['CHP' + str(b)]).sum()

        else:
            heatGenTec['CHPGasImpact' + str(b)] = 0
            heatGenTec['CHPImpact' + str(b)] = 0

        ### impact of electricity input
        elecImpactBuilding = ElecInfluenceBuilding(dataDict, buildings, capImpact[b,"PVElecCo2"],
                                                   capImpact[b,"CHPElecCo2"], elecGenTec['impact'], b, 'co2')

        ## GWHP
        GWHPData = tec_Data(sheetTrans, str(b), 'GWHP')
        capImpact[b,"GWHPHeatCo2"] = capTec.loc["GWHP", "kW in " + str(b)] * GWHPData['impact_cap'] /(GWHPData['lifetime'].mean())
        heatGenTec['GWHPHeatElecImpact' + str(b)] = elecImpactBuilding['Co2Impact_B' + str(b)] * hpElecInput['GWHP_B' + str(b)]
        heatGenTec['GWHPImpact' + str(b)] = heatGenTec['GWHPHeatElecImpact' + str(b)] + capImpact[b,"GWHPHeatCo2"]*heatGenTec['GWHP' + str(b)]/(heatGenTec['GWHP' + str(b)]).sum()

        ## HP
        HPData = tec_Data(sheetTrans, str(b), 'HP')
        capImpact[b,"HPHeatCo2"] = capTec.loc["HP", "kW in " + str(b)] * HPData['impact_cap'] /(HPData['lifetime'].mean())
        heatGenTec['HPHeatElecImpact' + str(b)] = elecImpactBuilding['Co2Impact_B' + str(b)] * hpElecInput['HP_B' + str(b)]
        heatGenTec['HPImpact' + str(b)] = heatGenTec['HPHeatElecImpact' + str(b)] + capImpact[b,"HPHeatCo2"]*heatGenTec['HP' + str(b)]/(heatGenTec['HP' + str(b)]).sum()

        ## ER
        try:
            ERData = tec_Data(sheetTrans, str(b), 'ElectricRod')
            capImpact[b,"ElectricRodHeatCo2"] = capTec.loc["ElectricRod", "kW in " + str(b)] * ERData['impact_cap'] /(ERData['lifetime'].mean())
            heatGenTec['ElectricRodHeatElecImpact' + str(b)] = elecImpactBuilding['Co2Impact_B' + str(b)] * hpElecInput['ElectricRod_B' + str(b)]
            heatGenTec['ElectricRodImpact' + str(b)] = heatGenTec['ElectricRodHeatElecImpact' + str(b)] + capImpact[b,"ElectricRodHeatCo2"]*heatGenTec['ElectricRod' + str(b)]/(heatGenTec['ElectricRod' + str(b)]).sum()
        except:
            None

        ## Gas Boiler
        GasBoilerData = tec_Data(sheetTrans, str(b), 'GasBoiler')
        capImpact[b,"GasBoilerHeatCo2"] = capTec.loc["GasBoiler", "kW in " + str(b)] * GasBoilerData['impact_cap'] /(GasBoilerData['lifetime'].mean())
        heatGenTec['GasBoilerHeatGasImpact' + str(b)] = gasImpactBuilding * heatGenTec['GasBoiler' + str(b)] / float(GasBoilerData['efficiency'].split(",")[0])
        heatGenTec['GasBoilerImpact' + str(b)] = heatGenTec['GasBoilerHeatGasImpact' + str(b)] + capImpact[b,"GasBoilerHeatCo2"]*heatGenTec['GasBoiler' + str(b)]/(heatGenTec['GasBoiler' + str(b)]).sum()

        ## Solar Collector
        SolarCollectorData = tec_Data(sheetSolar, str(b), 'solarCollector')
        capImpact[b,"SolarCollectorHeatCo2"] = capTec.loc["solarConnectBus", "kW in " + str(b)] * SolarCollectorData['impact_cap'] /(SolarCollectorData['lifetime'].mean())
        heatGenTec['SolarCollectorImpact' + str(b)] =capImpact[b,"SolarCollectorHeatCo2"]*heatGenTec['solarCollectordhw_B' + str(b)]/(heatGenTec['solarCollectordhw_B' + str(b)]).sum()

    elecGenTec['GridImpactTotal'] = elecGenTec[['GridImpact' + str(b) for b in range(1,buildings+1)]].sum(axis=1)
    elecGenTec['CHPGasImpactTotal'] = elecGenTec[['CHPGasImpact' + str(b) for b in range(1,buildings+1)]].sum(axis=1)
    elecGenTec['CHPElecImpactTotal'] = elecGenTec[['CHPImpact' + str(b) for b in range(1,buildings+1)]].sum(axis=1)
    elecGenTec['PVImpactTotal'] = elecGenTec[['PVImpact' + str(b) for b in range(1,buildings+1)]].sum(axis=1)

    capImpact['system', "PVElecCo2"] = sum([capImpact[b,"PVElecCo2"] for b in range(1,buildings)])
    capImpact['system', "CHPElecCo2"] = sum([capImpact[b, "CHPElecCo2"] for b in range(1,buildings)])

    elecImpactBuilding = ElecInfluenceBuilding(dataDict, buildings, capImpact['system',"CHPElecCo2"],
                                               capImpact['system',"PVElecCo2"], elecGenTec['impact'], 'system', 'co2')

    heatGenTec['GWHPImpactTotal'] = heatGenTec[['GWHPImpact' + str(b) for b in range(1,buildings+1)]].sum(axis=1)
    heatGenTec['CHPGasImpactTotal'] = heatGenTec[['CHPGasImpact' + str(b) for b in range(1,buildings+1)]].sum(axis=1)
    heatGenTec['CHPHeatImpactTotal'] = heatGenTec[['CHPImpact' + str(b) for b in range(1,buildings+1)]].sum(axis=1)
    heatGenTec['GasBoilerImpactTotal'] = heatGenTec[['GasBoilerImpact' + str(b) for b in range(1,buildings+1)]].sum(axis=1)
    heatGenTec['HPImpactTotal'] = heatGenTec[['HPImpact' + str(b) for b in range(1,buildings+1)]].sum(axis=1)
    try:
        heatGenTec['ElectricRodImpactTotal'] = heatGenTec[['ElectricRodImpact' + str(b) for b in range(1,buildings+1)]].sum(axis=1)
    except:
        None
    heatGenTec['SolarCollectorImpactTotal'] = heatGenTec[['SolarCollectorImpact' + str(b) for b in range(1,buildings+1)]].sum(axis=1)

    return heatGenTec, elecGenTec, elecImpactBuilding


def co2_balance_barplot(dataDict, inputFileName, buildings, selected_days, elecEmission, iter, outputFileName):

    """
    Function to GHG emission impact for different technologies
    :param :
            dataDict: full result file
            inputFileName: type str of link to input file
            outputFileName: type str of link to output file
            building: number of buildings
            selected_days: str of datetime for selected day
            iter: optimization iteration
            elecEmission: time serie or float with GHG emission impact of grid electricity

    :return: plotting GHG impact
    """

    heatGenTec, elecGenTec, _ = co2_balance(dataDict, inputFileName, buildings, selected_days, gasEmission, elecEmission, rangeToConsider)

    elecGenTecToPlot = elecGenTec[['CHPElecImpactTotal', 'PVImpactTotal', 'GridImpactTotal']]

    fig = plt.figure()
    (elecGenTecToPlot.rename(columns=labelDict)).plot(
        kind='bar', stacked=True, title="CO2 impact for electricity supply", xlabel='Hours', ylabel='CO2 equivance per source in ')

    plt.legend(loc=(1.04, 0) ) # ['CHP', 'PV', 'Grid']
    ax = plt.gca()
    ax.axes.xaxis.set_ticklabels([])
    plt.savefig(outputFileName + "Optimization" + str(iter) + '/hourly_co1impact_elec.png', bbox_inches='tight')
    fig.clf()
    plt.close()

    heatGenTecToPlot = heatGenTec[['CHPHeatImpactTotal', 'GWHPImpactTotal', 'HPImpactTotal', 'ElectricRodImpactTotal', 'GasBoilerImpactTotal',
            'SolarCollectorImpactTotal']]

    fig = plt.figure()
    (heatGenTecToPlot.rename(columns=labelDict)).plot(
        kind='bar', stacked=True, title="CO2 impact for heat supply", xlabel='Hours', ylabel='CO2 equivance per source in ')

    plt.legend(loc=(1.04, 0) )
    ax = plt.gca()
    ax.axes.xaxis.set_ticklabels([])
    plt.savefig(outputFileName + "Optimization" + str(iter) + '/hourly_co1impact_heat.png', bbox_inches='tight')
    fig.clf()
    plt.close()

    ### Combine heat and elec in one bar chart
    try:
        heatElecGenTec = pd.concat([heatGenTec[['CHPHeatImpactTotal', 'GWHPImpactTotal', 'HPImpactTotal',
                                                'ElectricRodImpactTotal','GasBoilerImpactTotal', 'SolarCollectorImpactTotal']],
                                    elecGenTec[['CHPElecImpactTotal', 'PVImpactTotal', 'GridImpactTotal']]], axis=1)
    except:
        heatElecGenTec = pd.concat([heatGenTec[['CHPHeatImpactTotal', 'GWHPImpactTotal', 'HPImpactTotal',
                                                'GasBoilerImpactTotal', 'SolarCollectorImpactTotal']],
                                    elecGenTec[['CHPElecImpactTotal', 'PVImpactTotal', 'GridImpactTotal']]], axis=1)

    heatElecGenTec['CHPImpactTotal'] = heatElecGenTec['CHPElecImpactTotal'] + heatElecGenTec['CHPHeatImpactTotal']

    fig = plt.figure()
    (heatElecGenTec.rename(columns=labelDict)).plot(
        kind='bar', stacked=True, title="CO2 impact for heat and electricity supply", xlabel='Hours', ylabel='CO2 equivance per source in ')

    plt.legend(loc=(1.04, 0)) #
    ax = plt.gca()
    ax.axes.xaxis.set_ticklabels([])
    plt.savefig(outputFileName + "Optimization" + str(iter) + '/hourly_co1impact.png', bbox_inches='tight')
    fig.clf()
    plt.close()


    for days in selected_days:
        fig = plt.figure()
        (elecGenTec[['CHPElecImpactTotal', 'PVImpactTotal', 'GridImpactTotal']].rename(columns=labelDict).loc[
         days:days + pd.offsets.DateOffset(hour=23)]).reset_index(drop=True).plot(
            kind='bar', stacked=True, title="CO2 impact for electricity supply", xlabel='Hours', ylabel='CO2 equivance per source in ')

        plt.legend(loc=(1.04, 0))
        #ax = plt.gca()
        #ax.axes.xaxis.set_ticklabels([])
        plt.savefig(outputFileName + "Optimization" + str(iter) + '/hourly_co1impact_elec_' + str(days)[0:10] + '.png', bbox_inches='tight')
        fig.clf()
        plt.close()

        fig = plt.figure()
        try:
            (heatGenTec[['CHPHeatImpactTotal', 'GWHPImpactTotal', 'HPImpactTotal', 'ElectricRodImpactTotal', 'GasBoilerImpactTotal',
                         'SolarCollectorImpactTotal']].rename(columns=labelDict).loc[days:days + pd.offsets.DateOffset(hour=23)]).reset_index(drop=True).plot(
                 kind='bar', stacked=True, title="CO2 impact for heat generation", xlabel='Hours', ylabel='CO2 equivance per source in ')
        except:
            heatElecGenTec = pd.concat([heatGenTec[['CHPHeatImpactTotal', 'GWHPImpactTotal', 'HPImpactTotal',
                                                    'GasBoilerImpactTotal', 'SolarCollectorImpactTotal']],
                                        elecGenTec[['CHPElecImpactTotal', 'PVImpactTotal', 'GridImpactTotal']]], axis=1)
        plt.legend(loc=(1.04, 0) ) # ['CHP', 'GWHP', 'HP', 'Gas Boiler', 'Solar collector']

        plt.savefig(outputFileName + "Optimization" + str(iter) + '/hourly_co1impact_heat_' + str(days)[0:10] + '.png', bbox_inches='tight')
        fig.clf()
        plt.close()

        fig = plt.figure()
        try:
            (heatElecGenTec[['CHPImpactTotal', 'GWHPImpactTotal', 'HPImpactTotal', 'ElectricRodImpactTotal',
                             'GasBoilerImpactTotal', 'SolarCollectorImpactTotal',
                             'PVImpactTotal', 'GridImpactTotal']].rename(columns=labelDict).loc[days:days + pd.offsets.DateOffset(hour=23)]).reset_index(drop=True).plot(
                kind='bar', stacked=True, title="CO2 impact for heat and electricity generation", xlabel='Hours',
                ylabel='CO2 equivance per source in ')
        except:
            (heatElecGenTec[['CHPImpactTotal', 'GWHPImpactTotal', 'HPImpactTotal',
                             'GasBoilerImpactTotal', 'SolarCollectorImpactTotal',
                             'PVImpactTotal', 'GridImpactTotal']].rename(columns=labelDict).loc[days:days + pd.offsets.DateOffset(hour=23)]).reset_index(drop=True).plot(
                kind='bar', stacked=True, title="CO2 impact for heat and electricity generation", xlabel='Hours',
                ylabel='CO2 equivance per source in ')
        plt.legend(loc=(1.04, 0) ) # ['CHP', 'GWHP', 'HP', 'Gas Boiler', 'Solar collector', 'PV', 'Grid']
        plt.savefig(outputFileName + "Optimization" + str(iter) + '/hourly_co1impact_' + str(days)[0:10] + '.png', bbox_inches='tight')
        fig.clf()
        plt.close()


# Flexibility KPIs


def grid_periods(dataDict, inputFileName, buildings, gasCost, elecCost, gasEmission, elecEmission, impactDayPeriode,
                 selected_days, k_value, impactPara, optMode, parameter, rangeToConsider):
    """
    Function to GHG emission impact for different technologies
    :param :
        dataDict: full result file
        inputFileName: type str of link to input file
        building: number of buildings
        selected_days: str of datetime for selected day
        gasCost: time serie or float with cost of gas
        elecCost: time serie or float with cost of electricity
        gasEmission: time serie or float with GHG emission impact of gas
        elecEmission: time serie or float with GHG emission impact of grid electricity
        k_value: type float [0,1] to define high and low value periode
        impactPara: electricity price considertion: type str gridData, systemdata
        optMode: str group or indiv
        parameter: flexibility base parameter type str cost or co2
        rangeToConsider: time serie of datetime

    :return: dict with grid periods
    """


    print('determine grid periods')
    gridPeriods = {}

    if parameter == 'cost':
        sheetCol = 'variable costs'
        sheetImpact = "cost"

        print('costs')
    elif parameter == 'co2':
        sheetCol = 'CO2 impact'
        sheetImpact = 'impact'
        print('co2 equivalent impact')

    if impactPara == 'gridData':
        if parameter == 'cost':
            if type(elecCost) == str:
                elecImpact = pd.read_csv(elecCost, sep=';', usecols = ["timestamp", sheetImpact], index_col ='timestamp')
                elecImpact['timestamp'] = elecImpact.index
                elecImpact['timestamp'] = pd.to_datetime(pd.date_range(start=elecImpact['timestamp'].iloc[0],
                                                                       end=elecImpact['timestamp'].iloc[-1],
                                                                       freq='H'), yearfirst=True, dayfirst=False)
                elecImpact = elecImpact.set_index('timestamp')
                elecImpact = elecImpact.loc[rangeToConsider,:]
            elif type(elecCost) == int or type(elecCost) == float:
                elecImpact = pd.DataFrame(index=rangeToConsider)
                elecImpact["cost"] = elecCost
                elecImpact.index.name = "timestamp"
                elecImpact = elecImpact.loc[rangeToConsider,:]

        elif parameter == 'co2':
            if type(elecEmission) == str:
                elecImpact = pd.read_csv(elecEmission, sep=';', usecols = ["timestamp", sheetImpact], index_col ='timestamp')
                elecImpact['timestamp'] = elecImpact.index
                elecImpact['timestamp'] = pd.to_datetime(pd.date_range(start=elecImpact['timestamp'].iloc[0],
                                                                       end=elecImpact['timestamp'].iloc[-1],
                                                                       freq='H'), yearfirst=True, dayfirst=False)
                elecImpact = elecImpact.set_index('timestamp')
                elecImpact = elecImpact.loc[rangeToConsider,:]
            elif type(elecEmission) == int or type(elecEmission) == float:
                elecImpact = pd.DataFrame(index=rangeToConsider)
                elecImpact["impact"] = elecEmission
                elecImpact.index.name = "timestamp"
                elecImpact = elecImpact.loc[rangeToConsider,:]


    elif impactPara == 'systemData':
        if parameter == 'cost':
            _, __, ___, elecImpact = npc_technology(dataDict, inputFileName, buildings, gasCost, elecCost, optMode)

            elecImpact['timestamp'] = pd.to_datetime(elecImpact.index, yearfirst=True, dayfirst=False)
            elecImpact = elecImpact.set_index('timestamp')
            elecImpact = elecImpact.rename(columns={'ElecPriceSystem': 'cost'})
            elecImpact = elecImpact.loc[rangeToConsider,:]

        elif parameter == 'co2':
            _, __, elecImpact = co2_balance(dataDict, inputFileName, buildings, selected_days, gasEmission, elecEmission, rangeToConsider)

            elecImpact['timestamp'] = pd.to_datetime(elecImpact.index, yearfirst=True, dayfirst=False)
            elecImpact = elecImpact.set_index('timestamp')
            elecImpact = elecImpact.rename(columns={'Co2ImpactSystem': 'impact'})
            elecImpact = elecImpact.loc[rangeToConsider,:]

    #iterate days over year
    elecImpact["timestamp"] = elecImpact.index
    for day in elecImpact.timestamp.dt.date.unique():
        startDate = pd.to_datetime(str(day) + ' 00:00:00', yearfirst=True, dayfirst=False)
        endDate = pd.to_datetime(str(day) + ' 23:00:00', yearfirst=True, dayfirst=False)

        if pd.to_datetime(str(day) + ' 23:00:00' , yearfirst=True, dayfirst=False) in elecImpact.index:
            elecImpactPeriode = elecImpact.loc[startDate:endDate, :]

        else:
            print('index is NOT in index')

        upperBoundaryLP = round(elecImpactPeriode.sort_values(by=sheetImpact).iloc[0:int(len(elecImpactPeriode)*k_value),
                                :][sheetImpact].max(),4)
        lowerBoundaryLP = round(elecImpactPeriode[sheetImpact].min(),4)

        lowerBoundaryHP = round(elecImpactPeriode.sort_values(by=sheetImpact).iloc[int(len(elecImpactPeriode)*(1-k_value)):,
                                :][sheetImpact].min(),4)
        upperBoundaryHP = round(elecImpactPeriode[sheetImpact].max(),4)

        elecImpactPeriode = elecImpactPeriode.round(4)

        ## LP
        if not gridPeriods:
            # create dict if not exists
            gridPeriods[parameter, 'lowPeriods'] = elecImpactPeriode[
                (elecImpactPeriode.loc[:, sheetImpact] <= upperBoundaryLP) &
                (elecImpactPeriode.loc[:, sheetImpact] >= lowerBoundaryLP) &
                (elecImpactPeriode.loc[:, sheetImpact] < upperBoundaryHP)
                ]

            ## HP
            gridPeriods[parameter, 'highPeriods'] = elecImpactPeriode[
                (elecImpactPeriode.loc[:, sheetImpact] <= upperBoundaryHP) &
                (elecImpactPeriode.loc[:, sheetImpact] >= lowerBoundaryHP) &
                (elecImpactPeriode.loc[:, sheetImpact] > lowerBoundaryLP)
                ]
        else:
            # expand dict
            gridPeriods[parameter, 'lowPeriods'] = pd.concat([gridPeriods[parameter, 'lowPeriods'], elecImpactPeriode[
                (elecImpactPeriode.loc[:, sheetImpact] <= upperBoundaryLP) &
                (elecImpactPeriode.loc[:, sheetImpact] >= lowerBoundaryLP) &
                (elecImpactPeriode.loc[:, sheetImpact] < upperBoundaryHP)
                ]])

            ## HP
            gridPeriods[parameter, 'highPeriods'] = pd.concat([gridPeriods[parameter, 'highPeriods'], elecImpactPeriode[
                (elecImpactPeriode.loc[:, sheetImpact] <= upperBoundaryHP) &
                (elecImpactPeriode.loc[:, sheetImpact] >= lowerBoundaryHP) &
                (elecImpactPeriode.loc[:, sheetImpact] > lowerBoundaryLP)
                ]])

    return gridPeriods, elecImpact


def energy_flexibility(dataDict, inputFileName, buildings, gasCost, elecCost, gasEmission, elecEmission, optMode,
                       impactPara, selected_days, iter, iterRange, rangeToConsider, results):
    """
    Function to GHG emission impact for different technologies
    :param :
        dataDict: full result file
        inputFileName: type str of link to input file
        building: number of buildings
        selected_days: str of datetime for selected day
        gasCost: time serie or float with cost of gas
        elecCost: time serie or float with cost of electricity
        gasEmission: time serie or float with GHG emission impact of gas
        elecEmission: time serie or float with GHG emission impact of grid electricity
        k_value: type float [0,1] to define high and low value periode
        impactPara: electricity price considertion: type str gridData, systemdata
        iter: optimization iteration
        iterRange: list of ploting order of optimization iterations
        optMode: str group or indiv
        rangeToConsider: time serie of datetime
        results: dict where iteration results are stored

    :return: grid periods
    """

    elecInput = unit_elecInput(dataDict, buildings)
    k_value = 0.4
    impactDayPeriode = 1
    elecImpactPara = {}
    gridPeriodsCost, elecImpactPara['cost'] = grid_periods(dataDict, inputFileName, buildings, gasCost, elecCost, gasEmission,
                                                       elecEmission, impactDayPeriode, selected_days, k_value, impactPara,
                                                       optMode, 'cost', rangeToConsider)

    gridPeriodsCo2, elecImpactPara['CO2'] = grid_periods(dataDict, inputFileName, buildings,  gasCost, elecCost, gasEmission,
                                                            elecEmission, impactDayPeriode, selected_days, k_value, impactPara,
                                                            optMode, 'co2', rangeToConsider)
    gridPeriods = {**gridPeriodsCost, **gridPeriodsCo2}
    gridPeriodPlot(gridPeriods, elecImpactPara, optMode)

    if iter == iterRange[0]:
        results['flexibilityFactor', 'co2'] = {}
        results['flexibilityFactor', 'cost'] = {}
    else:
        None

    if len(gridPeriods['co2', 'lowPeriods']) == 0:
        print('############# co2 -lowPeriods are empty ###########')
    if len(gridPeriods['co2', 'highPeriods']) == 0:
        print('############# co2 -highPeriods are empty###########')
    if len(gridPeriods['cost', 'lowPeriods']) == 0:
        print('############# cost -lowPeriods are empty###########')
    if len(gridPeriods['cost', 'highPeriods']) == 0:
        print('############# cost -highPeriods are empty###########')

    if len(gridPeriods['co2', 'lowPeriods']) == 0 & len(gridPeriods['co2', 'highPeriods']) == 0:
        results['flexibilityFactor', 'co2'][iter] = 0
        print('flexbility Factor for emissions is set to 0, no periods for emissions detected')
    else:

        results['flexibilityFactor', 'co2'][iter] = (float(elecInput.loc[elecInput.index.isin(gridPeriods['co2', 'lowPeriods'].index),
                                                                         'HPTotal'].sum()) \
                                                     - float(elecInput.loc[elecInput.index.isin(gridPeriods['co2', 'highPeriods'].index),
                                                                           'HPTotal'].sum() )
                                                     ) / (
                                                     float(elecInput.loc[elecInput.index.isin(gridPeriods['co2', 'lowPeriods'].index),
                                                                                'HPTotal'].sum()) \
                                                     + float(elecInput.loc[elecInput.index.isin(gridPeriods['co2', 'highPeriods'].index),
                                                                                  'HPTotal'].sum())
                                                    )

    if len(gridPeriods['cost', 'lowPeriods']) == 0 & len(gridPeriods['cost', 'highPeriods']) == 0:
        results['flexibilityFactor', 'cost'][iter] = 0
        print('flexbility Factor for costs is set to 0, no periods for cost detected')
    else:
        results['flexibilityFactor', 'cost'][iter] = (float(elecInput.loc[elecInput.index.isin(gridPeriods['cost', 'lowPeriods'].index),
                                                                          'HPTotal'].sum()) \
                                                      - float(elecInput.loc[elecInput.index.isin(gridPeriods['cost', 'highPeriods'].index),
                                                                            'HPTotal'].sum() )
                                                      ) / (
                                                             float(elecInput.loc[elecInput.index.isin(gridPeriods['cost', 'lowPeriods'].index),
                                                                                 'HPTotal'].sum()) \
                                                             + float(elecInput.loc[elecInput.index.isin(gridPeriods['cost', 'highPeriods'].index),
                                                                                   'HPTotal'].sum())
                                                     )


    return results

def gridPeriodPlot(gridPeriods, elecImpactPara, optMode):
    fig = plt.figure()
    plt.plot(elecImpactPara['cost']['cost'])
    plt.scatter(gridPeriods['cost', 'highPeriods'].index, gridPeriods['cost', 'highPeriods']['cost'])
    plt.scatter(gridPeriods['cost', 'lowPeriods'].index, gridPeriods['cost', 'lowPeriods']['cost'])
    fig.clf()
    plt.close()

def flexibilityBarChart(dataDict, inputFileName, buildings, gasCost, elecCost, gasEmission, elecEmission,
                        optMode, outputFileName, selected_days, iter, iterRange, rangeToConsider, results):

    """
    Function to GHG emission impact for different technologies
    :param :
        dataDict: full result file
        inputFileName: type str of link to input file
        building: number of buildings
        selected_days: str of datetime for selected day
        gasCost: time serie or float with cost of gas
        elecCost: time serie or float with cost of electricity
        gasEmission: time serie or float with GHG emission impact of gas
        elecEmission: time serie or float with GHG emission impact of grid electricity
        k_value: type float [0,1] to define high and low value periode
        impactPara: electricity price considertion: type str gridData, systemdata
        iter: optimization iteration
        iterRange: list of ploting order of optimization iterations
        optMode: str group or indiv
        rangeToConsider: time serie of datetime
        results: dict where iteration results are stored

    :return: grid periods
    """

    for impactPara in ['systemData', 'gridData']:#
        results = energy_flexibility(dataDict, inputFileName, buildings, gasCost, elecCost, gasEmission, elecEmission,
                                     optMode, impactPara, selected_days, iter, iterRange, rangeToConsider, results)

        if iter == iterRange[-1]:

            groups = [[results['flexibilityFactor', 'cost'][i], results['flexibilityFactor', 'co2'][i]] for i in iterRange]
            group_legend = ['cost', 'co2 impact']
            group_labels = results['flexibilityFactor', 'cost'].keys()
            df = pd.DataFrame(groups, columns = group_legend, index=group_labels) #

            fig = plt.figure()
            df.plot.bar()
            plt.title('Energy flexibility (-1: inflexible, +1: flexible ')
            plt.xlabel('Optimization iteration')
            plt.ylabel('Energy flexilibity parameter based on ' + impactPara)
            plt.savefig(outputFileName + '/energy_flexiblity' + impactPara + '.png', bbox_inches='tight')
            fig.clf()
            plt.close()

            fig = plt.figure()
            pd.DataFrame([results['flexibilityFactor', 'cost'][i] for i in iterRange], index=group_labels).plot.bar()
            plt.title('Flexibility factor: electricity costs (-1: inflexible, +1:  flexible ')
            plt.xlabel('Optimization iteration')
            plt.ylabel('Energy flexilibity factor based on ' + impactPara)
            plt.legend().remove()
            plt.savefig(outputFileName + '/cost_flexiblity' + impactPara + '.png', bbox_inches='tight')
            fig.clf()
            plt.close()

            fig = plt.figure()
            pd.DataFrame([results['flexibilityFactor', 'co2'][i] for i in iterRange], index=group_labels).plot.bar()
            plt.title('Flexibility factor: electricity CO2 equ. emissions (-1: inflexible, +1:  flexible ')
            plt.xlabel('Optimization iteration')
            plt.ylabel('Energy flexilibity factor based on ' + impactPara)
            plt.legend().remove()
            plt.savefig(outputFileName + '/c02_flexiblity' + impactPara + '.png', bbox_inches='tight')
            fig.clf()
            plt.close()


def tec_Data(sheet, building, tech):
    tecData = sheet[(sheet['label'] == tech) & (
            sheet['building'] == int(building))].iloc[0]
    return tecData

