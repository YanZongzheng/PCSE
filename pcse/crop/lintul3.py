# -*- coding: utf-8 -*-
# Copyright (c) 2004-2017 Alterra, Wageningen-UR
# Wim de Winter(Wim.deWinter@wur.nl), April 2015
# Allard de Wit, October 2017
"""
LINTUL3
"""
from math import exp

from ..base_classes import SimulationObject, ParamTemplate, RatesTemplate
from ..base_classes import StatesWithImplicitRatesTemplate as StatesTemplate

from ..traitlets import Float, AfgenTrait, Instance, Bool
from ..decorators import prepare_rates, prepare_states
from ..util import limit
from ..crop.phenology import DVS_Phenology as Phenology
from ..exceptions import CarbonBalanceError, NitrogenBalanceError
from .. import signals


# some functions to make unit conversion explicit.
def cm2mm(x):
    return x*10.


def joule2megajoule(x):
    return x/1e6


def m2mm(x):
    return x*1000.


class Lintul3Partitioning(SimulationObject):
    """Compute partitioning factors for LINTUL3.

    Biomass formed at any time during crop growth is partitioned amongst organs
    (Fig. 1), i.e. roots, stems, leaves and storage organs, with partitioning
    factors defined as a function of development stage (Fig. 2) (Drenth et al.,
    1994), which thus provides the rates of growth of these organs:

    dW/dt[i] = Pc[i] * dW / dt

    where (dW/dt) is the rate of biomass growth (gm-2 d-1); (dW/dt)[i] and Pc[i] are
    the rate of growth (gm-2 d-1) of and the biomass partitioning factor to organ i
    (g organ-i g-1 biomass), respectively. Leaf, stem and root weights of the
    seedlings at the time of transplanting are input parameters for the model. The
    time course of weights of these organs follows from integration of their net
    growth rates, i.e. growth rates minus death rates, the latter being defined as a
    function of physiological age, shading and stress.

    Modification are applied to root/shoot ratio in case of water/nitrogen stress.

    *parameters*
    ======== =============================================== =======  ==========
     Name     Description                                     Type     Unit
    ======== =============================================== =======  ==========
    FLVTB    Partitioning coefficients
    FRTTB    Partitioning coefficients
    FSOTB    Partitioning coefficients
    FSTTB    Partitioning coefficients
    NPART    Coefficient for the effect of N stress on leaf
             biomass reduction                                          -
    ======== =============================================== =======  ==========
    """

    class Parameters(ParamTemplate):
        FLVTB = AfgenTrait()  # Partitioning coefficients
        FRTTB = AfgenTrait()  # Partitioning coefficients
        FSOTB = AfgenTrait()  # Partitioning coefficients
        FSTTB = AfgenTrait()  # Partitioning coefficients
        NPART = Float(-99)  # Coefficient for the effect of N stress on leaf biomass reduction

    def initialize(self, day, kiosk, parvalues):
        self.params = self.Parameters(parvalues)

    def __call__(self, day, drv):
        """Dry matter partitioning fractions to leaves, stem and storage organs.

        Obsolete subroutine name: SUBPAR
        """
        k = self.kiosk
        p = self.params
        FRTWET = p.FRTTB(k.DVS)
        FLVT = p.FLVTB(k.DVS)
        FSTT = p.FSTTB(k.DVS)
        FSOT = p.FSOTB(k.DVS)

        if k.TRANRF < k.NNI:
            #  Water stress is more severe as compared to nitrogen stress and
            #  partitioning will follow the original assumptions of LINTUL2*
            FRTMOD = max(1., 1. / (k.TRANRF + 0.5))
            FRT = FRTWET * FRTMOD
            FSHMOD = (1. - FRT) / (1. - FRT / FRTMOD)
            FLV = FLVT * FSHMOD
            FST = FSTT * FSHMOD
            FSO = FSOT * FSHMOD
        else:
            # Nitrogen stress is more severe as compared to water stress and the
            # less partitioning to leaves will go to the roots*
            FLVMOD = exp(-p.NPART * (1.0 - k.NNI))
            FLV = FLVT * FLVMOD
            MODIF = (1. - FLV) / (1. - (FLV / FLVMOD))
            FST = FSTT * MODIF
            FRT = FRTWET * MODIF
            FSO = FSOT * MODIF

        return FRT, FLV, FST, FSO


class Lintul3NetGrowthRate(SimulationObject):
    """Compute the net growth rate as LINTUL3 does not take into account
    assimilation and respiration separately.

    Monteith (1977), Gallagher and Biscoe (1978) and Monteith (1990) have
    shown that biomass formed per unit intercepted light, LUE (Light
    Use Efficiency, g dry matter MJ-1), is relatively more stable. Hence,
    maximum daily growth rate can be defined as the product of
    intercepted PAR (photosynthetically active radiation, MJm-2 d-1)
    and LUE. Intercepted PAR depends on incident solar radiation,
    the fraction that is photosynthetically active (0.5) (Monteith and
    Unsworth, 1990; Spitters, 1990), and LAI (m2 leaf m-2 soil) according
    to Lambert-Beer's law:

      Q = 0.5Q0[1 - e(-k LAI)]

    where Q is intercepted PAR(MJm-2 d-1), Q0 is daily global radiation
    (MJm-2 d-1), and k is the attenuation coefficient for PAR in the
    canopy.

    Obsolete subroutine name: GROWTH

    *parameters*
    ======== =============================================== =======  ==========
     Name     Description                                     Type     Unit
    ======== =============================================== =======  ==========
    K        Light attenuation coefficient                              m²/m²
    LUE      Light use efficiency                                       g/MJ
    NLAI     Coefficient for the effect of N stress on LAI
             reduction(during juvenile phase)                           -
    NLUE     Coefficient of reduction of LUE under nitrogen
             stress, epsilon                                            -
    ======== =============================================== =======  ==========

    *External variables:*
    =========== ================================================= ==== ===============
     Name        Description                                      Pbl      Unit
    =========== ================================================= ==== ===============
     TRANRF      Transpiration reduction factor calculated         Y      -
     NNI         Nitrogen Nutrition Index                          Y      -
    =========== ================================================= ==== ===============
    """

    class Parameters(ParamTemplate):
        K = Float(-99)  # light extinction coefficient
        LUE = Float(-99)  # Light use efficiency.
        NLUE = Float(-99)  # Extinction coefficient for  Nitrogen distribution down the canopy

    def initialize(self, day, kiosk, parvalues):
        self.params = self.Parameters(parvalues)

    def __call__(self, day, drv):
        p = self.params
        k = self.kiosk

        DTR = joule2megajoule(drv.IRRAD)
        PARINT = 0.5 * DTR * (1. - exp(-p.K * k.LAI))
        RGROWTH = p.LUE * PARINT

        # Reduction factor due to nitrogen stress
        NRF = exp(-p.NLUE * (1.0 - k.NNI))
        # apply strongest reduction factor to net growth rate
        RGROWTH *= min(k.TRANRF, NRF)

        return RGROWTH


class Lintul3NitrogenStress(SimulationObject):
    """Computes nitrogen stress

    A crop is assumed to experience N stress at N concentrations below a critical
    value for unrestricted growth. To quantify crop response to nitrogen shortage, a
    Nitrogen Nutrition Index (NNI) is defined, ranging from 0 (maximum N shortage)
    to 1 (no N shortage):

    NNI = (actual crop [N] - residual [N]) / (critical [N] - residual [N])

    Critical crop nitrogen concentration, the lower limit of canopy nitrogen
    concentration in leaves and stems required for unrestricted growth, has been
    taken as half the maximum nitrogen concentration. An experimental basis for such
    an assumption can be derived from the study of Zhen and Leigh (1990), who
    reported that nitrate accumulation in plant occurs in significant quantity when
    the N needs to reach the maximum growth were fulfilled and the mean value of
    nitrate accumulated beyond the critical N concentration was about 50% for
    different stages.

    *parameters*
    ======== =============================================== =======  ==========
     Name     Description                                     Type     Unit
    ======== =============================================== =======  ==========
    FRNX     Critical N, as a fraction of maximum N
             concentration                                              -
    LRNR     Maximum N concentration of root as fraction of
             that of leaves                                             g/g
    LSNR     Maximum N concentration of stem as fraction of
             that of leaves                                             g/g
    RNFLV    Residual N concentration in leaves                         g/g
    RNFRT    Residual N concentration in roots                          g/g
    RNFST    Residual N concentration in stem                           g/g
    ======== =============================================== =======  ==========
    """

    class Parameters(ParamTemplate):
        FRNX = Float(-99)  # Optimal N concentration as the fraction of maximum N concentration in leaves.
        LRNR = Float(-99)  # Maximum concentration in roots as fraction of maximum N concentration in leaves
        LSNR = Float(-99)  # Maximum N concentration in stems as fraction of maximum N concentration in leaves
        RNFLV = Float(-99)  # Residual N concentration in leaves
        RNFRT = Float(-99)  # Residual N concentration in roots.
        RNFST = Float(-99)  # Residual N concentration in stems
        NMXLV = AfgenTrait()  # Maximum N concentration in the leaves as a function of development stage.

    def initialize(self, day, kiosk, parvalues):
        self.params = self.Parameters(parvalues)

    @prepare_rates
    def __call__(self, day, drv):
        p = self.params
        k = self.kiosk

        # Average residual N concentration.
        NRMR = (k.WLVG * p.RNFLV + k.WST * p.RNFST) / k.TBGMR

        # Maximum N concentration in the leaves, from which the values of the
        # stem and roots are derived, as a function of development stage.
        NMAXLV = p.NMXLV(k.DVS)
        NMAXST = p.LSNR * NMAXLV

        # maximum nitrogen concentration of leaves and stem.
        NOPTLV = p.FRNX * NMAXLV
        NOPTST = p.FRNX * NMAXST

        # Maximum N content in the plant.
        NOPTS = NOPTST * k.WST
        NOPTL = NOPTLV * k.WLVG
        NOPTMR = (NOPTL + NOPTS) / k.TBGMR

        # Total N in green matter of the plant.
        NUPGMR = k.ANLV + k.ANST

        # Nitrogen Nutrition Index.
        NFGMR = NUPGMR / k.TBGMR
        NNI = limit(0.001, 1.0, ((NFGMR - NRMR) / (NOPTMR - NRMR)))

        return NNI


class Lintul3NitrogenDynamics(SimulationObject):
    """Implements crop nitrogen dynamics in LINTUL3.

    The nitrogen dynamics implemented in this module include
    * nitrogen demand and uptake
    * nitrogen translocation from roots, stems & leaves and the uptake in the grains

    Finally, the N balance is computed at each time step. Note that the N dynamics in
    soil is not implemented here.

        *parameters*
        ======== =============================================== =======  ==========
         Name     Description                                     Type     Unit
        ======== =============================================== =======  ==========
        DVSNLT   Development stage after which no nutrients are absorbed    -
        DVSNT    development stage above which N translocation to
                 storage organs does occur                                  -
        FNTRT    Nitrogen translocation from roots to storage
                 organs as a fraction of total amount of
                 nitrogen translocated from leaves and stem to
                 storage organs                                             -
        FRNX     Critical N, as a fraction of maximum N
                 concentration                                              -
        LRNR     Maximum N concentration of root as fraction of
                 that of leaves                                             g/g
        LSNR     Maximum N concentration of stem as fraction of
                 that of leaves                                             g/g
        NMAXSO   Maximum concentration of nitrogen in storage
                 organs                                                     g/g
        RNFLV    Residual N concentration in leaves                         g/g
        RNFRT    Residual N concentration in roots                          g/g
        RNFST    Residual N concentration in stem                           g/g
        TCNT     Time coefficient for N translocation                       d
        ======== =============================================== =======  ==========

        *Tabular parameters*
        ======== =============================================== =======  ==========
         Name     Description                                     Type     Unit
        ======== =============================================== =======  ==========
        NMXLV    Maximum N concentration in the leaves, from
                 which the values of the stem and roots are derived,
                 as a  function of development stage
        ======== =============================================== =======  ==========


        *Initial states*
        ======== =============================================== =======  ==========
         Name     Description                                     Type     Unit
        ======== =============================================== =======  ==========
        NFRLVI   Initial fraction of N in leaves                            gN/gDM
        NFRRTI   Initial fraction of N in roots                             gN/gDM
        NFRSTI   Initial fraction of N in stem                              gN/gDM
        WLVGI    Initial Weight of green leaves                             g/m²
        WSTI     Initial Weight of stem                                     g/m²
        WRTLI    Initial Weight of roots                                    g/m²
        WSOI     Initial Weight of storage organs                           g/m²
        ======== =============================================== =======  ==========


        *State variables*
        =========== ================================================= ==== ===============
         Name        Description                                      Pbl      Unit
        =========== ================================================= ==== ===============
        ANLV        Actual N content in leaves                                 gN/m²
        ANRT        Actual N content in root                                   gN/m²
        ANSO        Actual N content in storage organs                         gN/m²
        ANST        Actual N content in stem                                   gN/m²
        NLOSSL      total N loss by leaves                                     gN/m²
        NLOSSR      total N loss by roots                                      gN/m²
        NUPTT       Total uptake of N over time                                gN/m²
        =========== ================================================= ==== ===============
    """


    # placeholders for initial nitrogen contents for leaves, stems, roots and storage organs
    ANLVI = 0.
    ANSTI = 0.
    ANRTI = 0.
    ANSOI = 0.

    class Parameters(ParamTemplate):
        DVSNLT = Float(-99)  # development stage N-limit
        DVSNT = Float(-99)  # development stage N-threshold
        FNTRT = Float(-99)  # Nitrogen translocation from roots as a fraction of the total amount of nitrogen translocated from leaves and stem.
        FRNX = Float(-99)  # Optimal N concentration as the fraction of maximum N concentration in leaves.
        LRNR = Float(-99)  # Maximum concentration in roots as fraction of maximum N concentration in leaves
        LSNR = Float(-99)  # Maximum N concentration in stems as fraction of maximum N concentration in leaves
        RNFLV = Float(-99)  # Residual N concentration in leaves
        RNFRT = Float(-99)  # Residual N concentration in roots.
        RNFST = Float(-99)  # Residual N concentration in stems
        NMXLV = AfgenTrait()  # Maximum N concentration in the leaves as a function of development stage.
        NMAXSO = Float(-99)  # Maximum N concentration in storage organs
        TCNT = Float(-99)  # Time coefficient(days) for N translocation.
        NFRLVI = Float(-99)  # Initial fraction of N (g N g-1 DM) in leaves.
        NFRRTI = Float(-99)  # Initial fraction of N (g N g-1 DM) in roots.
        NFRSTI = Float(-99)  # Initial fraction of N (g N g-1 DM) in stem.
        WLVGI = Float(-99)  # Initial weight of green leaves
        WSTI = Float(-99)  # Initial weight of stems
        WRTLI = Float(-99)  # Initial weight of roots
        WSOI = Float(-99)  # Initial weight of storage organs
        WCWP = Float(-99)  # Water content at wilting point (1.5 MPa) m3/ m3

    class StateVariables(StatesTemplate):
        ANLV = Float(-99.)  # Actual N content in leaves
        ANST = Float(-99.)  # Actual N content in stem
        ANRT = Float(-99.)  # Actual N content in root
        ANSO = Float(-99.)  # Actual N content in storage organs
        NUPTT = Float(-99.)  # Total uptake of N over time (g N m-2)
        NLOSSL = Float(-99.)  # total N loss by leaves
        NLOSSR = Float(-99.)  # total N loss by roots
        NDEMTO = Float(-99)  # Total nitrogen demand

    def initialize(self, day, kiosk, parvalues):
        self.params = self.Parameters(parvalues)
        p = self.params

        # Initial amount of N (g/m2) in leaves, stem, roots, and storage organs.
        self.ANLVI = p.NFRLVI * p.WLVGI
        self.ANSTI = p.NFRSTI * p.WSTI
        self.ANRTI = p.NFRRTI * p.WRTLI
        self.ANSOI = 0.0
        init = self.StateVariables.initialValues()
        init["NNI"] = 1.0
        init["ANLV"] = self.ANLVI
        init["ANST"] = self.ANSTI
        init["ANRT"] = self.ANRTI
        init["ANSO"] = self.ANSOI

        self.states = self.StateVariables(kiosk, publish=["NDEMTO", "ANLV", "ANST", "ANRT", "ANSO"],
                                          NNI=1.0, ANLV=self.ANLVI, ANST=self.ANSTI, ANRT=self.ANRTI,
                                          ANSO=self.ANSOI, NUPTT=0., NLOSSL=0., NLOSSR=0., NDEMTO=0.,)
        self.states.initialize_rates()

    @prepare_rates
    def calc_rates(self, day, drv):
        p = self.params
        s = self.states
        k = self.kiosk

        # translocatable N in leaves, stem, roots and storage organs.
        ATNLV = max(0., s.ANLV - k.WLVG * p.RNFLV)
        ATNST = max(0., s.ANST - k.WST * p.RNFST)
        ATNRT = min((ATNLV + ATNST) * p.FNTRT, s.ANRT - k.WRT * p.RNFRT)
        ATN = ATNLV + ATNST + ATNRT

        # Loss of nitrogen due to plant death rates
        RNLDLV = p.RNFLV * k.DLV
        RNLDRT = p.RNFRT * k.DRRT

        # Nitrogen demand
        #
        # Total crop nitrogen demand equals the sum of the nitrogen demands of its
        # individual organs (excluding storage organs, for which nitrogen demand is met by
        # translocation from the other organs, i.e. roots, stems and leaves) (Fig. 3).
        # Nitrogen demand of the individual organs is calculated as the difference
        # between maximum and actual organ nitrogen contents. The maximum nitrogen content
        # is defined as a function of canopy development stage (Drenth et al., 1994).
        # Total N demand (TNdem: gm-2 d-1) of the crop is:
        #
        # TNdem = sum(Nmax,i - ANi / dt)
        #
        # where Nmax,i is the maximum nitrogen concentration of organ i (gN/g biomass,
        # with i referring to leaves, stems and roots), Wi is the weight of organ i
        # (g biomass/m2), and ANi is the actual nitrogen content of organ i (gN/m2).
        #
        # N demand of leaves, roots and stem storage organs.
        # Maximum N concentration in the leaves, from which the values of the
        # stem and roots are derived, as a function of development stage.
        NMAXLV = p.NMXLV(k.DVS)
        NMAXST = p.LSNR * NMAXLV
        NMAXRT = p.LRNR * NMAXLV
        NDEMLV = max(NMAXLV * k.WLVG - s.ANLV, 0.)
        NDEMST = max(NMAXST * k.WST - s.ANST, 0.)
        NDEMRT = max(NMAXRT * k.WRT - s.ANRT, 0.)
        NDEMSO = max(p.NMAXSO * k.WSO - s.ANSO, 0.) / p.TCNT

        # Total Nitrogen demand.
        NDEMTO = max(0.0, (NDEMLV + NDEMST + NDEMRT))

        # Nitrogen uptake
        #
        # About 75-90% of the total N uptake at harvest takes place before
        # anthesis and, in conditions of high soil fertility, post-anthesis N uptake
        # may contribute up to 25% but would exclusively end up in the grain
        # as protein. Therefore, this nitrogen would not play any role in the
        # calculation of nitrogen stress that influences the biomass formation.
        # Therefore, nitrogen uptake is assumed to cease at anthesis,
        # as nitrogen content in the vegetative parts hardly increases after
        # anthesis

        # N supply to the storage organs.
        NSUPSO = ATN / p.TCNT if (k.DVS > p.DVSNT) else 0.0

        # Rate of N uptake in grains.
        RNSO = min(NDEMSO, NSUPSO)

        #  Nitrogen uptake limiting factor at low moisture conditions in the
        #  rooted soil layer before anthesis. After anthesis there is no
        #  uptake from the soil anymore.
        NLIMIT = 1.0 if (k.DVS < p.DVSNLT) and (k.WC >= p.WCWP) else 0.0

        DELT = 1.0
        TNSOIL = self.kiosk["TNSOIL"]
        NUPTR = (max(0., min(NDEMTO, k.TNSOIL)) * NLIMIT) / DELT

        # N translocated from leaves, stem, and roots.
        RNTLV = RNSO * ATNLV / ATN
        RNTST = RNSO * ATNST / ATN
        RNTRT = RNSO * ATNRT / ATN

        # compute the partitioning of the total N uptake rate (NUPTR) over the leaves, stem and roots.
        RNULV, RNUST, RNURT = self.N_uptake_rates(NDEMLV, NDEMST, NDEMRT, NUPTR, NDEMTO)

        RNST = RNUST - RNTST
        RNRT = RNURT - RNTRT - RNLDRT

        # Rate of change of N in organs
        RNLV = RNULV - RNTLV - RNLDLV

        NBALAN = (s.NUPTT + (self.ANLVI + self.ANSTI + self.ANRTI + self.ANSOI)
                  - (s.ANLV + s.ANST + s.ANRT + s.ANSO + s.NLOSSL + s.NLOSSR))
        if abs(NBALAN) > 0.0001:
            raise NitrogenBalanceError("Nitrogen un-balance in crop model at day %s" % day)

        # Assign computed rates
        s.rANLV = RNLV
        s.rANST = RNST
        s.rANRT = RNRT
        s.rANSO = RNSO
        s.rNUPTT = NUPTR
        s.rNLOSSL = RNLDLV
        s.rNLOSSR = RNLDRT

    @prepare_states
    def integrate(self, day):

        s = self.states
        p = self.params
        k = self.kiosk

        s.integrate(delta=1.0)

        # N demand of leaves, roots and stem storage organs.
        # Maximum N concentration in the leaves, from which the values of the
        # stem and roots are derived, as a function of development stage.
        NMAXLV = p.NMXLV(k.DVS)
        NMAXST = p.LSNR * NMAXLV
        NMAXRT = p.LRNR * NMAXLV
        NDEMLV = max(NMAXLV * k.WLVG - s.ANLV, 0.)
        NDEMST = max(NMAXST * k.WST - s.ANST, 0.)
        NDEMRT = max(NMAXRT * k.WRT - s.ANRT, 0.)

        # Total Nitrogen demand.
        s.NDEMTO = max(0.0, (NDEMLV + NDEMST + NDEMRT))

    def N_uptake_rates(self, NDEML, NDEMS, NDEMR, NUPTR, NDEMTO):
        """Compute the partitioning of the total N uptake rate (NUPTR) over the leaves, stem, and roots.

        Obsolete subroutine name: RNUSUB
        """
        if NDEMTO > 0:
            RNULV = (NDEML / NDEMTO) * NUPTR
            RNUST = (NDEMS / NDEMTO) * NUPTR
            RNURT = (NDEMR / NDEMTO) * NUPTR

            return RNULV, RNUST, RNURT
        else:
            return 0.0, 0.0, 0.0


class Lintul3(SimulationObject):
    """Implementation of the LINTUL3 crop simulation model.

    * ORIGINAL COPYRGIGHT NOTICE:
    *-------------------------------------------------------------------------*
    * Copyright 2013. Wageningen University, Plant Production Systems group,  *
    * P.O. Box 430, 6700 AK Wageningen, The Netherlands.                      *
    * You may not use this work except in compliance with the Licence.        *
    * You may obtain a copy of the Licence at:                                *
    *                                                                         *
    * http://models.pps.wur.nl/content/licence-agreement                      *
    *                                                                         *
    * Unless required by applicable law or agreed to in writing, software     *
    * distributed under the Licence is distributed on an "AS IS" basis,       *
    * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.*
    *-------------------------------------------------------------------------*

    LINTUL3 is a crop model that calculates biomass production based on intercepted photosynthetically
    active radiation (PAR) and light use efficiency (LUE). It is an adapted version of LINTUL2 (that simulates
    potential and water-limited crop growth), including nitrogen limitation. Nitrogen stress in the model is
    defined through the nitrogen nutrition index (NNI): the ratio of actual nitrogen concentration and critical
    nitrogen concentration in the plant. The effect of nitrogen stress on crop growth is tested in the model
    either through a reduction in LUE or leaf area (LA) or a combination of these two and further evaluated
    with independent datasets. However, water limitation is not considered in the present study as the
    crop is paddy rice. This paper describes the model for the case of rice, test the hypotheses of N stress
    on crop growth and details of model calibration and testing using independent data sets of nitrogen
    treatments (with fertilizer rates of 0 - 400 kgNha-1) under varying environmental conditions in Asia.
    Results of calibration and testing are compared graphically, through Root Mean Square Deviation (RMSD),
    and by Average Absolute Deviation (AAD). Overall average absolute deviation values for calibration and
    testing of total aboveground biomass show less than 26% mean deviation from the observations though
    the values for individual experiments show a higher deviation up to 41%. In general, the model responded
    well to nitrogen stress in all the treatments without fertilizer application as observed, but between
    fertilized treatments the response was varying.

    Nitrogen demand, uptake and stress

    At sub-optimal nitrogen availability in the soil, nitrogen demand of the crop
    cannot be satisfied, which leads to sub-optimal crop nitrogen concentration. The
    crop nitrogen concentration below which a crop experiences nitrogen stress is
    called the critical nitrogen concentration. Nitrogen stress results in reduced
    rates of biomass production and eventually in reduced yields. Actual N content
    is the accumulated N above residual (which forms part of the cell structure).
    The critical N content is the one corresponding to half of the maximum. Nitrogen
    contents of these three reference points include those of leaves and stems,
    whereas roots are not considered since N contents of above-ground (green) parts
    are more important for photosynthesis, because of their chlorophyll content.
    However, calculation of N demand and N uptake also includes the belowground
    part.

    See M.E. Shibu, P.A. Leffelaar, H. van Keulena, P.K. Aggarwal (2010). LINTUL3,
        a simulation model for nitrogen-limited situations: application to rice.
        Eur. J. Agron. (2010) http://dx.doi.org/10.1016/j.eja.2010.01.003

    available from http://models.pps.wur.nl/sites/models.pps.wur.nl/files/LINTUL-N-Shibu-article.pdf


    *parameters*
    ======== =============================================== =======  ==========
     Name     Description                                     Type     Unit
    ======== =============================================== =======  ==========
    DVSI     Initial development stage                                  -
    DVSDR    Development stage above which deathOfLeaves of
             leaves and roots start                                     -
    DVSNLT   Development stage after which no nutrients are absorbed    -
    DVSNT    development stage above which N translocation to
             storage organs does occur                                  -
    FNTRT    Nitrogen translocation from roots to storage
             organs as a fraction of total amount of
             nitrogen translocated from leaves and stem to
             storage organs                                             -
    FRNX     Critical N, as a fraction of maximum N
             concentration                                              -
    LAICR    critical LAI above which mutual shading of
             leaves occurs,                                             °C/d
    LRNR     Maximum N concentration of root as fraction of
             that of leaves                                             g/g
    LSNR     Maximum N concentration of stem as fraction of
             that of leaves                                             g/g
    NLAI     Coefficient for the effect of N stress on LAI
             reduction(during juvenile phase)                           -
    NMAXSO   Maximum concentration of nitrogen in storage
             organs                                                     g/g
    NSLA     Coefficient for the effect of N stress on SLA
             reduction                                                  -
    RDRNS    Relative death rate of leaf weight due to N
             stress                                                     1/d
    RDRRT    Relative death rate of roots                               1/d
    RDRSHM   and the maximum dead rate of leaves due to
             shading                                                    1/d
    RGRL     Relative growth rate of LAI at the exponential
             growth phase                                               °C/d
    RNFLV    Residual N concentration in leaves                         g/g
    RNFRT    Residual N concentration in roots                          g/g
    RNFST    Residual N concentration in stem                           g/g
    ROOTDM   Maximum root depth                                         m
    RRDMAX   Maximum rate of increase in rooting depth                  m/d
    SLAC     Specific leaf area constant                                m²/g
    TBASE    Base temperature for crop development                      °C
    TCNT     Time coefficient for N translocation                       d
    TRANCO   Transpiration constant indicating the level of
             drought tolerance of the wheat crop                        mm/d
    TSUMAG   Temperature sum for ageing of leaves                       °C.d
    WCFC     Water content at field capacity (0.03 MPa)                 m³/m³
    WCST     Water content at full saturation                           m³/m³
    WCWET    Critical Water content for oxygen stress                   m³/m³
    WCWP     Water content at wilting point (1.5 MPa)                   m³/m³
    WMFAC    water management (False = irrigated up to the
             field capacity, true= irrigated up to saturation)          (bool)

    *function tables*
    ======== =============================================== =======  ==========
     Name     Description                                     Type     Unit
    ======== =============================================== =======  ==========
    NMXLV    Maximum N concentration in the leaves, from
             which the values of the stem and roots are derived,
             as a  function of development stage
    RDRT     Relative death rate of leaves as a function of
             Developmental stage                                        1/d
    SLACF    Leaf area correction function as a function of
             development stage, DVS. Reference: Drenth, H.,
             ten Berge, H.F.M. and Riethoven, J.J.M. 1994,
             p.10. (Complete reference under Observed data.)
    ======== =============================================== =======  ==========


    * initial states *
    ======== =============================================== =======  ==========
     Name     Description                                     Type     Unit
    ======== =============================================== =======  ==========
    ROOTDI   Initial rooting depth                                      m
    NFRLVI   Initial fraction of N in leaves                            gN/gDM
    NFRRTI   Initial fraction of N in roots                             gN/gDM
    NFRSTI   Initial fraction of N in stem                              gN/gDM
    WCI      Initial water content in soil                              m³/³
    WLVGI    Initial Weight of green leaves                             g/m²
    WSTI     Initial Weight of stem                                     g/m²
    WRTLI    Initial Weight of roots                                    g/m²
    WSOI     Initial Weight of storage organs                           g/m²
    ======== =============================================== =======  ==========


    **State variables:**
    =========== ================================================= ==== ===============
     Name        Description                                      Pbl      Unit
    =========== ================================================= ==== ===============
    ANLV        Actual N content in leaves
    ANRT        Actual N content in root
    ANSO        Actual N content in storage organs
    ANST        Actual N content in stem
    CUMPAR      PAR accumulator
    LAI         leaf area index                                    *        m²/m²
    NLOSSL      total N loss by leaves
    NLOSSR      total N loss by roots
    NUPTT       Total uptake of N over time                                 gN/m²
    ROOTD       Rooting depth                                      *        m
    TNSOIL      Amount of inorganic N available for crop uptake
    WDRT        dead roots (?)                                              g/m²
    WLVD        Weight of dead leaves                                       g/m²
    WLVG        Weight of green leaves                                      g/m²
    WRT         Weight of roots                                             g/m²
    WSO         Weight of storage organs                                    g/m²
    WST         Weight of stem                                              g/m²
    TAGBM       Total aboveground biomass                                   g/m²
    TGROWTH     Total biomass growth (above and below ground)               g/m²
    =========== ================================================= ==== ===============

    **Rate variables:**
    =========== ================================================= ==== ===============
     Name        Description                                      Pbl      Unit
    =========== ================================================= ==== ===============
     PEVAP       Potential soil evaporation rate                   Y     |mmday-1|
     PTRAN       Potential crop transpiration rate                 Y     |mmday-1|
     TRAN        Actual crop transpiration rate                    N     |mmday-1|
     TRANRF      Transpiration reduction factor calculated         N      -
     RROOTD      Rate of root growth                               Y     |mday-1|
     NNI         Nitrogen Nutrition Index                          N      -
    =========== ================================================= ==== ===============
    """

    # sub-model components for crop simulation
    pheno = Instance(SimulationObject)
    nstress = Instance(SimulationObject)
    n_dynamics = Instance(SimulationObject)
    netgrowthrate = Instance(SimulationObject)
    partitioning = Instance(SimulationObject)
    # placeholder for effective N application rate from the _on_APPLY_N event handler.
    FERTNS = 0.0
    # placeholder for initial leaf area index
    LAII = 0.

    # Parameters, rates and states which are relevant at the main crop simulation level
    class Parameters(ParamTemplate):
        DVSI = Float(-99.)  # Development stage at start of the crop
        DVSDR = Float(-99)  # Development stage above which deathOfLeaves of leaves and roots start.
        DVSNLT = Float(-99)  # development stage N-limit
        LAICR = Float(-99)  # (oC d)-1, critical LAI above which mutual shading of leaves occurs,
        NLAI = Float(-99)  # Coefficient for the effect of N stress on LAI reduction(during juvenile phase)
        NSLA = Float(-99)  # Coefficient for the effect of N stress on SLA reduction
        RDRSHM = Float(-99)  # and the maximum relative deathOfLeaves rate of leaves due to shading.
        RGRL = Float(-99)  # Relative totalGrowthRate rate of LAI at the exponential totalGrowthRate phase
        ROOTDM = Float(-99)  # Maximum root depth for a rice crop.
        RRDMAX = Float(-99)  # Maximum rate of increase in rooting depth (m d-1) for a rice crop.
        SLAC = Float(-99)  # Specific leaf area constant.
        TBASE = Float(-99)  # Base temperature for spring wheat crop.
        TRANCO = Float(-99)  # Transpiration constant (mm/day) indicating the level of drought tolerance of the wheat crop.
        TSUMAG = Float(-99)  # Temperature sum for ageing of leaves
        WCFC = Float(-99)  # Water content at field capacity (0.03 MPa) m3/ m3
        WCI = Float(-99)  # Initial water content in cm3 of water/(cm3 of soil).
        WCST = Float(-99)  # Water content at full saturation m3/ m3
        WCWET = Float(-99)  # Critical Water conten for oxygen stress [m3/m3]
        WCWP = Float(-99)  # Water content at wilting point (1.5 MPa) m3/ m3
        WMFAC = Bool(False)  # water management (0 = irrigated up to the field capacity, 1= irrigated up to saturation)
        RDRNS = Float(-99)  # Relative deathOfLeaves rate of leaves due to N stress.
        RDRRT = Float(-99)  # Relative deathOfLeaves rate of roots.
        # FLVTB = AfgenTrait()  # Partitioning coefficients
        # FRTTB = AfgenTrait()  # Partitioning coefficients
        # FSOTB = AfgenTrait()  # Partitioning coefficients
        # FSTTB = AfgenTrait()  # Partitioning coefficients
        RDRT = AfgenTrait()  #
        SLACF = AfgenTrait()  # Leaf area correction function as a function of development stage, DVS.
        ROOTDI = Float(-99)  # initial rooting depth [m]
        WLVGI = Float(-99)  # Initial weight of green leaves
        WSTI = Float(-99)  # Initial weight of stems
        WRTLI = Float(-99)  # Initial weight of roots
        WSOI = Float(-99)  # Initial weight of storage organs
        
    class Lintul3States(StatesTemplate):
        LAI = Float(-99.)  # leaf area index
        WLVG = Float(-99.)  # Weight of green leaves
        WLVD = Float(-99.)  # Weight of dead leaves
        WST = Float(-99.)  # Weight of stem
        WSO = Float(-99.)  # Weight of storage organs
        WRT = Float(-99.)  # Weight of roots
        ROOTD = Float(-99.)  # Actual root depth [m]
        TGROWTH = Float(-99.)  # Total growth
        WDRT = Float(-99.)  # weight of dead roots
        TNSOIL = Float(-99.)  # Amount of inorganic N available for crop uptake.
        TAGBM = Float(-99.)  # Total aboveground biomass [g /m-2)
        TBGMR = Float(-99.)  # Total aboveground green (living) biomas

    # These are some rates which are not directly connected to a state (PEVAP, TRAN) of which must be published
    # (RROOTD) for the water balance module. Therefore, we explicitly define them here.
    class Lintul3Rates(RatesTemplate):
        PEVAP = Float()
        PTRAN = Float()
        TRAN = Float()
        TRANRF = Float()
        RROOTD = Float()
        DLV = Float()
        DRRT = Float()
        NNI = Float()

    def initialize(self, day, kiosk, parvalues):
        """
        :param day: start date of the simulation
        :param kiosk: variable kiosk of this PCSE  instance
        :param parvalues: `ParameterProvider` object providing parameters as
                key/value pairs
        """
        self.kiosk = kiosk
        self.params = self.Parameters(parvalues)
        self.rates = self.Lintul3Rates(self.kiosk,
                                       publish=["PEVAP", "TRAN", "RROOTD", "DLV", "DRRT", "NNI", "TRANRF"])

        self._connect_signal(self._on_APPLY_N, signals.apply_n)

        # Initialize phenology component of the crop
        self.pheno = Phenology(day, kiosk, parvalues)

        # Calculate initial LAI
        p = self.params
        SLACFI = p.SLACF(p.DVSI)
        ISLA = p.SLAC * SLACFI
        self.LAII = p.WLVGI * ISLA

        # Generate a dict with 'default' initial states (e.g. zero)
        init = self.Lintul3States.initialValues()

        # Initialize state variables
        init["LAI"] = self.LAII
        init["WLVG"] = p.WLVGI
        init["WST"] = p.WSTI
        init["WSO"] = p.WSOI
        init["WRT"] = p.WRTLI
        init["ROOTD"] = p.ROOTDI
        init["TBGMR"] = p.WLVGI + p.WSTI

        # Initialize the states objects
        self.states = self.Lintul3States(kiosk,
                                         publish=["LAI", "ROOTD", "WLVG", "WST",
                                                  "WRT", "WSO", "TBGMR", "TNSOIL"],
                                         **init)
        # Initialize the associated rates of the states
        self.states.initialize_rates()

        # Init sub module
        self.nstress = Lintul3NitrogenStress(day, kiosk, parvalues)
        self.n_dynamics = Lintul3NitrogenDynamics(day, kiosk, parvalues)
        self.netgrowthrate = Lintul3NetGrowthRate(day, kiosk, parvalues)
        self.partitioning = Lintul3Partitioning(day, kiosk, parvalues)

    def _on_APPLY_N(self, amount, recovery):
        """Receive signal for N application with amount the nitrogen amount in g N m-2 and
        recovery the recovery fraction.
        """
        self.FERTNS = amount * recovery

    @prepare_rates
    def calc_rates(self, day, drv):

        p = self.params
        s = self.states
        r = self.rates
        k = self.kiosk

        DELT = 1

        # Crop phenology
        #
        # Crop development, i.e. the order and rate of appearance of vegetative and
        # reproductive organs, is defined in terms of phenological developmental stage
        # (DVS) as a function of heat sum, which is the cumulative daily effective
        # temperature. Daily effective temperature is the average temperature above a
        # crop-specific base temperature (for rice 8C). Some crop or crop varieties are
        # photoperiod sensitive, i.e. flowering depends on the length of the light period
        # during the day in addition to the temperature during the vegetative stage.
        self.pheno.calc_rates(day, drv)

        # if before emergence there is no need to continue
        # because only the phenology is running.
        if self.pheno.get_variable("STAGE") == "emerging":
            return  # no above-ground crop to calculate yet

        # code below is executed only POST-emergence

        # Effective temperature above base temp
        DTEFF = max(0., drv.TEMP - p.TBASE)

        # potential rates of evaporation and transpiration:
        r.PEVAP, r.PTRAN = self._calc_potential_evapotranspiration(drv)

        # actual rate of transpiration:
        r.TRAN = self._calc_actual_transpiration(r.PTRAN, k.WC)

        # Nitrogen Nutrition Index.
        r.NNI = self.nstress(day, drv)

        # Growth reduction function for water stress(actual trans/potential)
        r.TRANRF = r.TRAN / r.PTRAN

        # total totalGrowthRate rate.
        RGROWTH = self.netgrowthrate(day, drv)

        # Biomass partitioning
        FRT, FLV, FST, FSO = self.partitioning(day, drv)

        # Leaf area development
        #
        # The time course of LAI is divided into two stages: an exponential stage during
        # the juvenile phase, where leaf area development is a function of temperature,
        # and a linear stage where it depends on the increase in leaf biomass (Spitters,
        # 1990; Spitters and Schapendonk, 1990). The death of leaves due to senescence
        # that may be enhanced by shading and/or stress leads to a corresponding loss in
        # leaf area. The specific leaf area is used for the conversion of dead leaf
        # biomass to corresponding loss in leaf area. The death of leaves due to
        # senescence occurs only after flowering and the rate depends on crop age
        # (function adopted from ORYZA2000, Bouman et al., 2001). The excessive growth of
        # leaves also result in death of leaves due to mutual shading. The death of leaves
        # due to shading is determined by a maximum death rate and the relative proportion
        # of leaf area above the critical LAI (4.0) (Spitters, 1990; Spitters and
        # Schapendonk, 1990). The net rate of change in leaf area (dLAI/dt) is the
        # difference between its growth rate and its death rate:
        #
        # dLAI/dt = dGLAI / dt - dDLAI / dt
        #
        # where (dGLAI/dt) is the leaf area growth rate and (dDLAI/dt) is the
        # leaf area death rate.

        # Leaf totalGrowthRate and LAI.
        GLV = FLV * RGROWTH

        # Specific Leaf area(m2/g).
        SLA = p.SLAC * p.SLACF(k.DVS) * exp(-p.NSLA * (1. - r.NNI))

        # daily increase of leaf area index.
        GLAI = self._growth_leaf_area(DTEFF, self.LAII, DELT, SLA, GLV, k.WC, k.DVS, r.TRANRF, r.NNI)

        # Relative and actual death rate of leaves due to senescence/ageing.
        RDRTMP = p.RDRT(drv.TEMP)
        DLV, DLAI = self._death_rate_of_leaves(k.TSUM, RDRTMP, r.NNI, SLA)

        # Net rate of change of Leaf area.
        RLAI = GLAI - DLAI

        # Root growth
        #
        # The root system is characterized by its vertical extension in the soil profile.
        # At emergence or at transplanting for transplanted rice. Effect of N stress (NNI)
        # on crop growth through its effect on LA and LUE. rooting depth is initialized.
        # Roots elongate at a constant daily rate, until flowering, provided soil water
        # content is above permanent wilting point (PWP), whereas growth ceases when soil
        # is drier than PWP or when a certain preset maximum rooting depth is reached
        # (Spitters and Schapendonk, 1990; Farr� et al., 2000).
        r.RROOTD = min(p.RRDMAX,  p.ROOTDM - s.ROOTD) if (k.WC > p.WCWP) else 0.0

        # Death of roots.
        DRRT = 0. if (k.DVS < p.DVSDR) else s.WRT * p.RDRRT

        # Actual growth rate of roots, leaves, stem
        # and storage organs.
        RWLVG = GLV - DLV
        RWRT = RGROWTH * FRT - DRRT
        RWST = RGROWTH * FST
        RWSO = RGROWTH * FSO

        # ****************SOIL NITROGEN SUPPLY***********************************
        # Soil–crop nitrogen balance
        #
        # The mineral nitrogen balance of the soil is the difference between nitrogen
        # added through mineralization and/or fertilizer, and nitrogen removed by crop
        # uptake and losses. The net rate of change of N in soil (dN/dt in gm-2 d-1) is:
        #
        # dN/dt[soil]= Nmin + (FERTN * NRF) - dNU/dt
        #
        # where Nmin is the nitrogen supply through mineralization and biological N
        # fixation, FERTN is the fertilizer nitrogen application rate, NRF is the
        # fertilizer nitrogen recovery fraction and dNU/dt is the rate of nitrogen uptake
        # by the crop, which is calculated as the minimum of the N supply from the soil
        # and the N demand from the crop.

        #  Soil N supply (g N m-2 d-1) through mineralization.
        NLIMIT = 1.0 if (k.DVS < p.DVSNLT) and (k.WC >= p.WCWP) else 0.0
        RTMIN = 0.10 * NLIMIT

        #  Change in inorganic N in soil as function of fertilizer
        #  input, soil N mineralization and crop uptake.
        NUPTR = (max(0., min(k.NDEMTO, s.TNSOIL)) * NLIMIT) / DELT
        RNSOIL = self.FERTNS/DELT - NUPTR + RTMIN
        self.FERTNS = 0.0

        # Carbon balance
        WLV = s.WLVG + s.WLVD  # Total leaf weight
        CBALANCE = (s.TGROWTH + (p.WRTLI + p.WLVGI + p.WSTI + p.WSOI)
                    - (WLV + s.WST + s.WSO + s.WRT + s.WDRT))
        if abs(CBALANCE) > 0.0001:
            raise CarbonBalanceError("Carbon un-balance in crop model at day %s" % day)

        # Assign associated rates
        s.rLAI = RLAI
        s.rWLVG = RWLVG
        s.rWLVD = DLV
        s.rWST = RWST
        s.rWSO = RWSO
        s.rWRT = RWRT
        s.rROOTD = r.RROOTD
        s.rTGROWTH = RGROWTH
        s.rWDRT = DRRT
        s.rTNSOIL = RNSOIL

        # Specific rate variable needed by nutrient module
        r.DLV = DLV  # rate of leaf death
        r.DRRT = DRRT  # rate of root death

        # Dynamics of the N module
        self.n_dynamics.calc_rates(day, drv)

    @prepare_states
    def integrate(self, day, delt=1.0):

        self.pheno.integrate(day, delt)

        # if before emergence there is no need to continue
        # because only the phenology is running.
        if self.pheno.get_variable("STAGE") == "emerging":
            # Just run a touch() to to ensure that all state variables are available
            # in the kiosk
            self.touch()
            return

        # run automatic integration on states object.
        s = self.states
        s.integrate(delta=delt)

        # Compute some derived states
        # Total above-ground biomass
        s.TAGBM = s.WLVG + s.WLVD + s.WST + s.WSO
        # Total living vegetative biomass.
        s.TBGMR = s.WLVG + s.WST

        self.n_dynamics.integrate(day)

    def _calc_potential_evapotranspiration(self, drv):
        """Compute potential evaporation and transpiration.
        """
        ES0 = cm2mm(drv.ES0)
        ET0 = cm2mm(drv.ET0)
        pevap = exp(-0.5 * self.states.LAI) * ES0
        pevap = max(0., pevap)
        ptran = (1. - exp(-0.5 * self.states.LAI)) * ET0
        ptran = max(0., ptran)
        return pevap, ptran

    def _calc_actual_transpiration(self, PTRAN, WC):
        """Compute actual rates of evaporation and transpiration.

        Obsolete subroutine name: EVAPTR
        """
        p = self.params
                
        # critical water content    
        WCCR = p.WCWP + max(0.01, PTRAN/(PTRAN + p.TRANCO) * (p.WCFC - p.WCWP))
          
        #  If the soil is flooded for rice, : the soil is assumed to be
        #  permanently saturated and there is no effect of the high water
        #  content on crop totalGrowthRate, because rice has aerenchyma.
        #  Thus FR is formulated as below:
        if p.WMFAC:
            if WC > WCCR:
                FR = 1.
            else:
                FR = limit(0., 1., (WC - p.WCWP) / (WCCR - p.WCWP))
        
        # If soil is irrigated but not flooded, : soil water content is
        # assumed to be at field capacity and the critical water content
        # that affects crop totalGrowthRate (FR) is formulated as below:
        else:
            if WCCR <= WC <= p.WCWET:  # Between critical
                FR = 1.0
            elif WC < WCCR:
                FR = limit(0., 1., (WC - p.WCWP)/(WCCR - p.WCWP))
            else:
                FR = limit(0., 1., (p.WCST - WC)/(p.WCST - p.WCWET))
        
        return PTRAN * FR

    def _growth_leaf_area(self, DTEFF, LAII,  DELT, SLA, GLV, WC, DVS, TRANRF, NNI):
        """Compute daily increase of leaf area index.
        """

        p = self.params
        LAI = self.states.LAI
        
        # Growth during maturation stage:
        GLAI = SLA * GLV

        # Growth during juvenile stage:
        if DVS < 0.2 and LAI < 0.75:
            GLAI = (LAI * (exp(p.RGRL * DTEFF * DELT) - 1.) / DELT) * TRANRF * exp(-p.NLAI * (1.0 - NNI))

        # Growth at day of seedling emergence:
        if LAI == 0. and WC > p.WCWP:
            GLAI = LAII / DELT  
        
        return GLAI

    # def _drymatter_partitioning_fractions(self, NPART, TRANRF, NNI):
        # """Dry matter partitioning fractions to leaves, stem and storage organs.
        #
        # Obsolete subroutine name: SUBPAR
        # """
        # k = self.kiosk
        # p = self.params
        # FRTWET = p.FRTTB(k.DVS)
        # FLVT = p.FLVTB(k.DVS)
        # FSTT = p.FSTTB(k.DVS)
        # FSOT = p.FSOTB(k.DVS)
        #
        # if TRANRF < NNI:
        #     #  Water stress is more severe as compared to nitrogen stress and
        #     #  partitioning will follow the original assumptions of LINTUL2*
        #     FRTMOD = max(1., 1./(TRANRF+0.5))
        #     FRT = FRTWET * FRTMOD
        #     FSHMOD = (1. - FRT) / (1. - FRT/FRTMOD)
        #     FLV = FLVT * FSHMOD
        #     FST = FSTT * FSHMOD
        #     FSO = FSOT * FSHMOD
        # else:
        #     # Nitrogen stress is more severe as compared to water stress and the
        #     # less partitioning to leaves will go to the roots*
        #     FLVMOD = exp(-NPART * (1.0 - NNI))
        #     FLV = FLVT * FLVMOD
        #     MODIF = (1. - FLV)/(1. - (FLV/FLVMOD))
        #     FST = FSTT * MODIF
        #     FRT = FRTWET * MODIF
        #     FSO = FSOT * MODIF
        #
        # return FRT, FLV, FST, FSO

    def _death_rate_of_leaves(self, TSUM, RDRTMP, NNI, SLA):
        """
        compute the relative deathOfLeaves rate of leaves due to age, 
        shading amd due to nitrogen stress.                        
        Obsolete subroutine name: DEATHL                 
        """
      
        p = self.params
        s = self.states
        
        RDRDV = 0. if (TSUM < p.TSUMAG) else RDRTMP
        
        RDRSH = max(0., p.RDRSHM * (s.LAI - p.LAICR) / p.LAICR)
        RDR = max(RDRDV, RDRSH)
        
        if NNI < 1.:
            DLVNS = s.WLVG * p.RDRNS * (1. - NNI)
            DLAINS = DLVNS * SLA
        else:
            DLVNS = 0.
            DLAINS = 0.
        
        DLVS = s.WLVG * RDR
        DLAIS = s.LAI * RDR
        
        DLV = DLVS + DLVNS
        DLAI = DLAIS + DLAINS
    
        return DLV, DLAI