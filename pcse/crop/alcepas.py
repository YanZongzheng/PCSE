# -*- coding: utf-8 -*-
# Copyright (c) 2004-2014 Alterra, Wageningen-UR
# Allard de Wit (allard.dewit@wur.nl), October 2017
"""Basic routines for ALCEPAS onion model
"""
from __future__ import print_function
from math import sqrt, exp, cos, pi

from ..traitlets import Instance, Float, AfgenTrait

from ..util import limit, astro, doy, daylength
from ..base_classes import ParamTemplate, StatesTemplate, RatesTemplate, SimulationObject
from ..decorators import prepare_rates, prepare_states
from .. import signals


class ALCEPAS_phenology(SimulationObject):

    class Parameters(ParamTemplate):
        pass

    class StateVariables(StatesTemplate):
        DVS = Float()
        BULBSUM = Float()
        BULB = Float()
        LAI = Float()

    class RateVariables(RatesTemplate):
        TEFF = Float()
        DAYFAC = Float()
        RFR = Float()
        RFRFAC = Float()
        DVR = Float()
        DTSUM = Float()

    def initialize(self, day, kiosk, parvalues):
        self.params = self.Parameters(parvalues)
        self.rates = self.RateVariables(day, kiosk)
        self.states = self.StateVariables(day, kiosk, DVS=0., BULBSUM=0.,
                                          BULB=0., LAI=5.)
        pass

    @prepare_rates
    def calc_rates(self, day, drv):
        r = self.rates
        s = self.states

        r.TEFF = max(0., drv.TEMP - 6.)
        DL = daylength(day, drv.LAT)
        r.DAYFAC = 0 if DL < 12 else (DL - 12.)/12
        # Change v into the right equation
        v = 0.9
        r.RFR = limit(0., 1., v)
        r.RFRFAC = 1 if r.RFR < 0.8 else 5*(1.0 - r.FRF)
        r.DTSUM = r.TEFF * r.DAYFAC * r.RFRFAC
        if s.states.DVS < 1.0:
            r.DVR = r.DTSUM/104.
        else:
            r.DVR = r.DTSUM/303.

    @prepare_states
    def integrate(self, day, delt=1.0):
        s = self.states
        r = self.rates
        s.BULBSUM += r.DTSUM * delt
        s.DVS += r.DVR * delt
        BULB = 0.3 + 100.45 * (exp(-exp(-0.0293*(s.BULBSUM - 91.9))))
        s.BULB = limit(0., 100., BULB)

        if s.DVS > 2:
            self._send_signal(signals.crop_finish)

