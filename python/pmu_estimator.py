from ctypes import CDLL, POINTER, c_bool, c_uint, Structure, byref, c_double
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

class EstimatorConfig(Structure):
    #_pack_ = 1  # Set structure alignment to 1 byte
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

            # Find the installed library
            library_name = 'pmu_estimator'
            lib_path = ctypes.util.find_library(library_name)

            if lib_path is None:
                plat = platform.system()

                if plat == "Linux":
                    lib_path = "/usr/local/lib/libpmu_estimator.so"
                elif plat == "Darwin":
                    lib_path = "/usr/local/lib/libpmu_estimator.dylib"
                elif plat == "Windows":
                    lib_path = "C:\Program Files (x86)\PmuEstimator\lib\libpmu_estimator.dll"
                    if not os.path.exists(str(lib_path)):
                        lib_path = "C:\Program Files\PmuEstimator\lib\libpmu_estimator.dll" 
                else:
                    raise ValueError(f"Unsupported platform: {plat}")
            

        if not os.path.exists(str(lib_path)):
            raise FileNotFoundError(f"Library not found at: {lib_path}")

        self.lib = CDLL(lib_path)

        self.lib.pmu_init.argtypes = [POINTER(EstimatorConfig), c_bool]
        self.lib.pmu_estimate.argtypes = [POINTER(c_double), c_double, POINTER(PmuFrame)]
        self.lib.pmu_deinit.argtypes = []

    def __del__(self):
        self.lib.pmu_deinit()

    def configure_from_ini(self, ini_file_path):
        self.lib.pmu_init(ini_file_path, self.CONFIG_FROM_INI)

    def configure_from_class(self, config):
        self.lib.pmu_init(byref(config), self.CONFIG_FROM_STRUCT)

    def estimate(self, input_signal_window, mid_window_fracsec):
        frame = PmuFrame()
        input_signal = (c_double * len(input_signal_window))(*input_signal_window)
        self.lib.pmu_estimate(input_signal, mid_window_fracsec , byref(frame))
        return frame
