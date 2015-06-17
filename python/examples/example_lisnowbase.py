import sys
import modele
from modele import f90
from pprint import pprint

params = modele.Lisnowbase_Mod.Lisnowparams()
xin = modele.Lisnowbase_Mod.Lisnowin()
xout = modele.Lisnowbase_Mod.Lisnowout()

modele.Lisnowbase_Mod.allocate_snow_adv(params, xin, xout)

print params
print xin
print xout
