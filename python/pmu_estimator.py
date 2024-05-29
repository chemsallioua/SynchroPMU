from ctypes import CDLL, POINTER, c_bool, c_uint, Structure, byref, c_double, c_int
import numpy as np
import os
import platform
import ctypes.util

class Phasor(Structure):
    _fields_ = [("amp", c_double),
                ("ph", c_double),
                ("freq", c_double)]

class PmuFrame(Structure):
    _fields_ = [("synchrophasor", Phasor),
                ("rocof", c_double)]

    def __init__(self, amp=0.0, ph=0.0, freq=0.0, rocof=0.0):
        self.synchrophasor = Phasor()
        self.synchrophasor.amp = amp
        self.synchrophasor.ph = ph
        self.synchrophasor.freq = freq
        self.rocof = rocof

    def __str__(self):
        return f"[Synchrophasor] amplitude: {self.synchrophasor.amp}, phase: {self.synchrophasor.ph}, frequency: {self.synchrophasor.freq}, rocof: {self.rocof}"


class SynchrophasorEstimatorParams(Structure):
    _fields_ = [("win_len", c_uint),
                ("n_cycles", c_uint),
                ("f0", c_uint),
                ("frame_rate", c_uint),
                ("fs", c_uint),
                ("n_bins", c_uint),
                ("P", c_uint),
                ("Q", c_uint),
                ("iter_eipdft_enabled", c_bool),
                ("interf_trig", c_double),
                ("df", c_double),
                ("norm_factor", c_double),
                ("phasor", Phasor)]

class InternalBuffers(Structure):
    _fields_ = [("Xf", POINTER(c_double)),
                ("Xi", POINTER(c_double)),
                ("dftbins", POINTER(c_double)),
                ("hann_window", POINTER(c_double)),
                ("signal_windows", POINTER(c_double) * 2)]

class HanningTransformConstants(Structure):
    _fields_ = [("C0", c_double),
                ("C1", c_double),
                ("C2", c_double),
                ("C3", c_double),
                ("C4", c_double),
                ("C5", c_double),
                ("C6", c_double),
                ("inv_norm_factor", c_double)]

class RocoFEstimationStates(Structure):
    _fields_ = [("freq_old", c_double * 2),
                ("thresholds", c_double * 3),
                ("low_pass_coeff", c_double * 3),
                ("delay_line", c_double * 2 * 2),
                ("state", c_bool * 2)]
    

class PmuContext(Structure):
    _fields_ = [("synch_params", SynchrophasorEstimatorParams),
                ("buff_params", InternalBuffers),
                ("hann_params", HanningTransformConstants),
                ("rocof_params", RocoFEstimationStates),
                ("pmu_initialized", c_bool)]
class EstimatorConfig(Structure):
    _fields_ = [("n_cycles", c_uint),
                ("f0", c_uint),
                ("frame_rate", c_uint),
                ("fs", c_uint),
                ("n_bins", c_uint),
                ("P", c_uint),
                ("Q", c_uint),
                ("iter_eipdft", c_bool),
                ("interf_trig", c_double),
                ("rocof_thresh", c_double * 3),
                ("rocof_low_pass_coeffs", c_double * 3)]

    def __init__(self, n_cycles, f0, frame_rate, fs, n_bins, P, Q, interf_trig, rocof_thresh, rocof_low_pass_coeffs, iter_eipdft = False):
        self.n_cycles = n_cycles
        self.f0 = f0
        self.frame_rate = frame_rate
        self.fs = fs
        self.n_bins = n_bins
        self.P = P
        self.Q = Q
        self.iter_eipdft = iter_eipdft
        self.interf_trig = c_double(interf_trig)
        self.rocof_thresh = (c_double * 3)(*rocof_thresh)
        self.rocof_low_pass_coeffs = (c_double * 3)(*rocof_low_pass_coeffs)

class PMUEstimator:

    # Values of the configuration modes
    CONFIG_FROM_INI = True
    CONFIG_FROM_STRUCT = False

    def __init__(self, lib_path = None):

        # If lib_path is not given, use default path for each OS
        if lib_path is None:
            
            plat = platform.system()

            if plat == "Linux":
                lib_path = "/usr/local/lib/libpmu_estimator.so"
            elif plat == "Darwin":
                lib_path = "/usr/local/lib/libpmu_estimator.dylib"
            elif plat == "Windows":
                lib_path = "C:\\Program Files (x86)\\PmuEstimator\\lib\\libpmu_estimator.dll"
                if not os.path.exists(str(lib_path)):
                    lib_path = "C:\\Program Files\\PmuEstimator\\lib\\libpmu_estimator.dll" 
            else:
                raise ValueError(f"Unsupported platform: {plat}")
            

        if not os.path.exists(str(lib_path)):
            raise FileNotFoundError(f"Library not found at: {lib_path}")

        self.lib = CDLL(lib_path)

        self.lib.pmu_init.argtypes = [POINTER(PmuContext), POINTER(EstimatorConfig), c_bool]
        self.lib.pmu_init.restype = c_int

        self.lib.pmu_estimate.argtypes = [POINTER(PmuContext), POINTER(c_double), c_double, POINTER(PmuFrame)]
        self.lib.pmu_estimate.restype = c_int

        self.lib.pmu_deinit.argtypes = [POINTER(PmuContext)]
        self.lib.pmu_deinit.restype = c_int

        self.ctx = PmuContext()
        self.ctx.pmu_initialized = 0

    def __del__(self):
        return self.lib.pmu_deinit(byref(self.ctx))
    
    def deinit(self):
        return self.lib.pmu_deinit(byref(self.ctx))

    def configure_from_ini(self, ini_file_path):
        return self.lib.pmu_init(byref(self.ctx), ini_file_path, self.CONFIG_FROM_INI)

    def configure_from_class(self, config):
        return self.lib.pmu_init(byref(self.ctx), byref(config), self.CONFIG_FROM_STRUCT)

    def estimate(self, input_signal_window, mid_window_fracsec):
        frame = PmuFrame()
        input_signal = (c_double * len(input_signal_window))(*input_signal_window)
        result = self.lib.pmu_estimate(byref(self.ctx), input_signal, c_double(mid_window_fracsec), byref(frame))

        framedict = {
            "amp": frame.synchrophasor.amp,
            "ph": frame.synchrophasor.ph,
            "freq": frame.synchrophasor.freq,
            "rocof": frame.rocof
        }

        if result != 0:
            return None
        else:
            return framedict

