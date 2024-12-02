from math import e

print( "    *** ISA calculator ***")
print()
print("1. Calculate ISA for altitude in meters")
print("2. Calculate ISA for altitude in feet")
print("3. Calculate ISA for altitude in FL")
print()

Unit = int(input("Enter your choice:"))
if Unit == 1:
    print()
    h = float(input("Enter altitude [m]:"))
if Unit == 2:
    print()
    h = 0.3048 * float(input("Enter altitude [ft]:"))
if Unit == 3:
    print()
    h = 100.0 * 0.3048 * float(input("Enter altitude [FL]:"))
if Unit < 1:
     print()
     print("Come on, you had one job!")
     print()
if Unit > 3:
     print()
     print("Come on, you had one job!")
     print()

#0<h<=11000.0
T = 288.15 - 0.0065*h
P = 101325.0*(T/288.15) ** (9.80665/(0.0065*287.0))
D = 101325/(287*288.15)*(T/288.15) ** (9.80665/(0.0065*287.0)-1)
Celsius = T - 273.15
Pre = 100*P/101325
Den = 100*D/1.225

#11000.0<h<=20000.0
h1 = 11000.0
T2 = 288.15 - 0.0065*h1
P1 = 101325.0*(T2/288.15) ** (9.80665/(0.0065*287.0))
D1 = 101325/(287*288.15)*(T2/288.15) ** (9.80665/(0.0065*287.0)-1)
P2 = P1 * e ** (-9.80665*(h-h1)/(T2*287.0))
D2 = D1 * e ** (-9.80665*(h-h1)/(T2*287.0))
Celsius2 = T2 - 273.15
Pre2 = 100*P2/101325
Den2 = 100*D2/1.225

#20000.0<h<=32000.0
#here I used the calculation for P2/D2 instead of the values, everywhere else, I copy-pasted the values from previous level.
h2 = 20000.0
T3 = T2 + 0.0010*(h-h2)
P3 = P1 * e ** (-9.80665*(h2-h1)/(T2*287.0))*(T3/T2) ** (-9.80665/(0.0010*287.0))
D3 = D1 * e ** (-9.80665*(h2-h1)/(T2*287.0))*(T3/T2) ** (-1-9.80665/(0.0010*287.0))
Celsius3 = T3 - 273.15
Pre3 = 100*P3/101325
Den3 = 100*D3/1.225

#32000.0<h<=47000.0
h3 = 32000.0
T4 = (T2 + 0.0010*(h3-h2)) + 0.0028*(h-h3)
P4 = 867.255 *(T4/(T2 + 0.0010*(h3-h2))) ** (-9.80665/(0.0028*287.0))
D4 = 0.013216*(T4/(T2 + 0.0010*(h3-h2))) ** (-1-9.80665/(0.0028*287.0))
Celsius4 = T4 - 273.15
Pre4 = 100*P4/101325
Den4 = 100*D4/1.225

#47000.0<h<=51000.0
h4 = 47000.0
T5 = 270.65
P5 = 110.7666 * e ** (-9.80665*(h-h4)/(T5*287.0))
D5 = 0.001426 * e ** (-9.80665*(h-h4)/(T5*287.0))
Celsius5 = T5 - 273.15
Pre5 = 100*P5/101325
Den5 = 100*D5/1.225

#51000.0<h<=71000.0
h5 = 51000.0
T6 = T5 - 0.0028*(h-h5)
P6 = 66.8483 *(T6/T5) ** (9.80665/(0.0028*287.0))
D6 = 0.000861*(T6/T5) ** (-1+9.80665/(0.0028*287.0))
Celsius6 = T6 - 273.15
Pre6 = 100*P6/101325
Den6 = 100*D6/1.225

#71000.0<h<=86000.0
h6 = 71000.0
T7 = 214.65 - 0.0020*(h-h6)
P7 = 3.949 *(T7/214.65) ** (9.80665/(0.0020*287.0))
D7 = 6.4e-05*(T7/214.65) ** (-1+9.80665/(0.0020*287.0))
Celsius7 = T7 - 273.15
Pre7 = 100*P7/101325
Den7 = 100*D7/1.225


if h < 0:
    print()
    print("Sorry, I can only do altitudes above sea level.")
    print()
    print("Ready.")
elif 0 <= h <= 11000.0:
    print()
    print("Temperature:" , round( T , 4) , "K" , "(" , round(Celsius , 1) ,  "'C)" )
    print("Pressure:" , round(P,4) , "Pa" , "(" , round(Pre, 2) , "% SL)")
    print("Density:" , round(D,6), "kg/m3" , "(" , round(Den, 2) , "% SL)")
    print()
    print("Ready.")
elif 11000.0 < h <= 20000:
    print()
    print("Temperature:" , round( T2 , 4) , "K" , "(" , round(Celsius2, 1) ,  "'C)" )
    print("Pressure:" , round(P2,4) , "Pa" , "(" , round(Pre2, 2) , "% SL)")
    print("Density:" , round(D2,6), "kg/m3" , "(" , round(Den2, 2) , "% SL)")
    print()
    print("Ready.")
elif 20000.0 < h <= 32000.0:
    print()
    print("Temperature:" , round( T3 , 4) , "K" , "(" , round(Celsius3 , 1) ,  "'C)" )
    print("Pressure:" , round(P3,4) , "Pa" , "(" , round(Pre3, 2) , "% SL)")
    print("Density:" , round(D3,6), "kg/m3" , "(" , round(Den3, 2) , "% SL)")
    print()
    print("Ready.")
elif 32000.0 < h <= 47000.0:
    print()
    print("Temperature:" , round( T4 , 4) , "K" , "(" , round(Celsius4 , 1) ,  "'C)" )
    print("Pressure:" , round(P4,4) , "Pa" , "(" , round(Pre4, 2) , "% SL)")
    print("Density:" , round(D4,6), "kg/m3" , "(" , round(Den4, 2) , "% SL)")
    print()
    print("Ready.")
elif 47000.0 < h <= 51000.0:
    print()
    print("Temperature:" , round( T5 , 4) , "K" , "(" , round(Celsius5 , 1) ,  "'C)" )
    print("Pressure:" , round(P5,4) , "Pa" , "(" , round(Pre5, 2) , "% SL)")
    print("Density:" , round(D5,6), "kg/m3" , "(" , round(Den5, 2) , "% SL)")
    print()
    print("Ready.")   
elif 51000.0 < h <= 71000.0:
    print()
    print("Temperature:" , round( T6 , 4) , "K" , "(" , round(Celsius6 , 1) ,  "'C)" )
    print("Pressure:" , round(P6,4) , "Pa" , "(" , round(Pre6, 2) , "% SL)")
    print("Density:" , round(D6,6), "kg/m3" , "(" , round(Den6, 2) , "% SL)")
    print()
    print("Ready.")
elif 71000.0 < h <= 86000.0:
    print()
    print("Temperature:" , round( T7 , 4) , "K" , "(" , round(Celsius7 , 1) ,  "'C)" )
    print("Pressure:" , round(P7,4) , "Pa" , "(" , round(Pre7, 2) , "% SL)")
    print("Density:" , round(D7,6), "kg/m3" , "(" , round(Den7, 2) , "% SL)")
    print()
    print("Ready.")
else:
    print()
    print("Sorry, I can only do altitudes up to 86000 m.")
    print()
    print("Ready.")
print()
dummy = input("Press enter to end the ISA calculator.")
