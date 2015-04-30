from ..traitlets import Float

from ..base_classes import ParamTemplate, SimulationObject, \
    VariableKiosk, StatesTemplate, prepare_states

class CropWaterBalanceAccumulator(SimulationObject):
    """This class does not do any calculations itself. It only ensures
    that a couple of water balance terms are aggregated over the course
    of the crop cycle instead of over the entire model run.

    To make this aggregation we only store the initial values of these
    variables (the "*I" variables) and subtract then from the same
    state variables in the finalize() call.
    """

    # For computing the totals of the water balance terms only over the
    # crop cycle, we first store the value at the start of the cropping
    # cycle
    _CEVSTI = Float()
    _CEVWTI = Float()
    _CTSRI = Float()
    _CRAINTI = Float()
    _CTOTINFI = Float()
    _CTOTIRRI = Float()
    _CPERCTI = Float()
    _CLOSSTI = Float()

    class StateVariables(StatesTemplate):
        CEVST = Float()
        CEVWT = Float()
        CTSR = Float()
        CRAINT = Float()
        CTOTINF = Float()
        CTOTIRR = Float()
        CPERCT = Float()
        CLOSST = Float()

    def initialize(self, day, kiosk, parvalues):
        self.kiosk = kiosk
        self.states = self.StateVariables(self.kiosk, CEVST=-99., CEVWT=-99., CTSR=-99., CRAINT=-99., CTOTINF=-99.,
                                          CTOTIRR=-99., CPERCT=-99., CLOSST=-99.)
        self._CEVSTI = self.kiosk["EVST"]
        self._CEVWTI = self.kiosk["EVWT"]
        self._CTSRI = self.kiosk["TSR"]
        self._CRAINTI = self.kiosk["RAINT"]
        self._CTOTINFI = self.kiosk["TOTINF"]
        self._CTOTIRRI = self.kiosk["TOTIRR"]
        self._CPERCTI = self.kiosk["PERCT"]
        self._CLOSSTI = self.kiosk["LOSST"]

    def calc_rates(self, day, drv):
        pass

    def integrate(self, day):
        pass

    @prepare_states
    def finalize(self, day):
        s = self.states
        s.CEVST = self.kiosk["EVST"] - self._CEVSTI
        s.CEVWT = self.kiosk["EVWT"] - self._CEVWTI
        s.CTSR = self.kiosk["TSR"] - self._CTSRI
        s.CRAINT = self.kiosk["RAINT"] - self._CRAINTI
        s.CTOTINF = self.kiosk["TOTINF"] - self._CTOTINFI
        s.CTOTIRR = self.kiosk["TOTIRR"] - self._CTOTIRRI
        s.CPERCT = self.kiosk["PERCT"] - self._CPERCTI
        s.CLOSST = self.kiosk["LOSST"] - self._CLOSSTI

