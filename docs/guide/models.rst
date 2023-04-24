 .. _Drescription of models :


- An active `Github <https://github.com/>`_ account to clone the repo.
- A solver is installed. `Gurobi solver <https://www.gurobi.com/resource/parallelism-linear-mixed-integer-programming/>`_ is recommended, although other solvers like CBC, GLPK, Cplex could also be used.

Models
===============

Transformers
-------------

As of now, optihood is available as an open source code and needs to be installed from source. Please follow the
Modelsm models modesl....and should be

1. Generic Combined transformers (GDY)


2. Air Source Heat Pumps (SPA)


3. Ground Source Heat pumps (SPA)


4. Combined heat and power (GDY)

5. Boilers (GDY)

6. Electric Rod (GDY)
--> link to generic model

Solar technologies
-------------
1. Solar thermal (SPA)
--> link to OEMOF
+ differential temperature levels

2. Photovoltaics (SPA)
--> link to OEMOF

Energy Storage
-------------
Optihood uses in its code, the modelling of components for energy storage.
The modelling of these storage components is done through the oemof-solph library.
The GenericStorage function is used to model a storage component with an input argument and an output argument.
Two classes were created, one for electrical storage and one for thermal storage.

1. Electric Batteries (GDY)

The modelling of the component for electrical storage is built from a class called "ElectricalStorage" which is a specificity of the "GenericStorage" class of the oemof.solph library. The "ElectricalStorage" class uses the "solph.Flow" object to represent the incoming and outgoing energy flows. This class includes a constructor with several arguments: the definition of the inputs and outputs of the system, the rate of energy loss, the initial storage level, the conversion efficiency of the incoming and outgoing energy, the minimum and maximum capacity of the storage, the investment costs and the distribution mode in particular.
The "solph.Investment" function is used in the code to represent the investment costs associated with an energy system modelled with the oemof.solph package. It is used to calculate the optimal investment costs in the context of solving the optimisation model. Environmental costs are also implemented.


2. Domestic hot water storage & hot water storage(GDY)

The modelling of the thermal storage component is built from a class called "ThermalStorage" which is a specificity of the "GenericStorage" class of the oemof.solph library. Like the "ElectricalStorage" class, this class uses the "solph.Flow" object. This class includes a constructor with several arguments: the definition of the inputs and outputs of the system, the rate of energy loss, the initial storage level, the conversion efficiency of the incoming and outgoing energy, the minimum and maximum capacity of the storage, the investment costs and the distribution mode in particular.
For the calculation of energy losses, a function called "precalculate" has been defined. This function calculates the heat losses and the U-value from the functions "calculate_storage_u_value" and "calculate_losses".
The "solph.Investment" function is used in the code to represent the investment costs associated with an energy system modelled with the oemof.solph package. It is used to calculate the optimal investment costs in the context of solving the optimisation model.

https://oemof-solph.readthedocs.io/en/latest/usage.html#genericstorage-component

Summary
-------------
image of table with nominal performances report WP1

.. image:: ./resources/Summary_Converters.PNG
      :width: 140
      :alt: constraint2
      :align: center

.. image:: ./guide/resources/Summary_Storage.PNG
      :width: 140
      :alt: constraint2
      :align: center