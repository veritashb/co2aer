# co2aer
R algorithm to calculate air exchange rates based on CO2 data

The proc_aer calculates air exchange (ACH) values using the tracer gas decay method. Details of this method are available in ASTM Std E-741 and ASTM Std Guide for using indoor CO2 concentrations to evaluate indoor air quality and ventilation.

The function takes a dataframe with the following fields as input:
- A vector with timestamp (in POSIXct format)
- A vector with CO2 concentration (in ppm)
- A sensor identification number

In addition, the function requires the user to supply a value for outdoor CO2 concentration (*amb* field in the function body). It is recommended to input a constant value, although it is also possible to work with a vector of outdoor CO2 concentration values. 

To ensure the quality of the estimated decay periods, the user can specify whether the detection of the start of each decay period follows a sharp increase in CO2 concentration, or not (*sharp.peak* "yes" or "no"). In case the user wants to limit the decay periods to a minimum delta in CO2 concentration, it can be done by indicating a value to the *min.delta* field in the function body. The user can also specify whether only decays in which CO2 concentration goes back close to outdoors (~50ppm above indicated outdoors) are considered for ACH estimation (*back.to.background* field, "yes" or "no").

The results are offered in two different output types: summary table or time series graph, (*output* = "summary_decay" or "graph_decay")

