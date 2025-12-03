"""
SynchroPMU Python API
=====================

Python bindings for the SynchroPMU C library, providing access to the Phasor Measurement
Unit (PMU) estimator based on the Iterative Enhanced Interpolated DFT algorithm.

This module provides a Pythonic interface to the C library, allowing users to:
- Initialize PMU estimator instances
- Configure from INI files or Python objects
- Estimate synchrophasors, frequency, and ROCOF from signal samples
- Manage multiple independent PMU instances

Example:
    Basic usage with configuration from structure::

        from pmu_estimator import PMUEstimator, EstimatorConfig
        
        # Create configuration
        config = EstimatorConfig(
            n_cycles=2, f0=50, frame_rate=50, fs=2000,
            n_bins=200, P=3, Q=2, interf_trig=0.1,
            rocof_thresh=[0.1, 0.5, 1.0],
            rocof_low_pass_coeffs=[0.1, 0.2, 0.3]
        )
        
        # Initialize estimator
        pmu = PMUEstimator()
        pmu.configure_from_class(config)
        
        # Process signal
        signal_window = [...]  # Your signal samples
        result = pmu.estimate(signal_window, mid_fracsec=0.5)
        print(f"Frequency: {result['freq']} Hz")

Author:
    Chemseddine Allioua, Brahim Mazighi

Copyright:
    Copyright (c) 2023. All Rights Reserved.
    Confidential and Proprietary - University of Bologna.
"""

from ctypes import CDLL, POINTER, c_bool, c_uint, Structure, byref, c_double, c_int, c_void_p, c_char_p, c_float
import numpy as np
import os
import platform
import ctypes.util

class Phasor(Structure):
    """
    Synchrophasor representation structure.
    
    Represents a complex voltage or current measurement as amplitude, phase, and frequency.
    
    Attributes:
        amp (float): Phasor amplitude (magnitude)
        ph (float): Phasor phase angle (radians)
        freq (float): Frequency (Hz)
    """
    _fields_ = [("amp", c_float),
                ("ph", c_float),
                ("freq", c_float)]

class PmuFrame(Structure):
    """
    PMU output frame structure.
    
    Contains the complete output of a single PMU estimation, including the
    synchrophasor and Rate of Change of Frequency (ROCOF).
    
    Attributes:
        synchrophasor (Phasor): Estimated synchrophasor
        rocof (float): Rate of Change of Frequency (Hz/s)
    """
    _fields_ = [("synchrophasor", Phasor),
                ("rocof", c_float)]

    def __init__(self, amp=0.0, ph=0.0, freq=0.0, rocof=0.0):
        """
        Initialize a PMU frame.
        
        Args:
            amp (float): Phasor amplitude
            ph (float): Phasor phase (radians)
            freq (float): Frequency (Hz)
            rocof (float): Rate of Change of Frequency (Hz/s)
        """
        self.synchrophasor = Phasor()
        self.synchrophasor.amp = amp
        self.synchrophasor.ph = ph
        self.synchrophasor.freq = freq
        self.rocof = rocof

    def __str__(self):
        """String representation of the PMU frame."""
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
                ("interf_trig", c_float),
                ("df", c_float),
                ("norm_factor", c_float),
                ("phasor", Phasor)]

class InternalBuffers(Structure):
    _fields_ = [("Xf", POINTER(c_float)),
                ("Xi", POINTER(c_float)),
                ("dftbins", POINTER(c_float)),
                ("hann_window", POINTER(c_float)),
                ("signal_windows", POINTER(c_float) * 2)]

class HanningTransformConstants(Structure):
    _fields_ = [("C0", c_float),
                ("C1", c_float),
                ("C2", c_float),
                ("C3", c_float),
                ("C4", c_float),
                ("C5", c_float),
                ("C6", c_float),
                ("inv_norm_factor", c_float)]

class RocoFEstimationStates(Structure):
    _fields_ = [("freq_old", c_float * 2),
                ("thresholds", c_float * 3),
                ("low_pass_coeff", c_float * 3),
                ("delay_line", c_float * 2 * 2),
                ("state", c_bool * 2)]
    

class PmuContext(Structure):
    _fields_ = [("synch_params", SynchrophasorEstimatorParams),
                ("buff_params", InternalBuffers),
                ("hann_params", HanningTransformConstants),
                ("rocof_params", RocoFEstimationStates),
                ("pmu_initialized", c_bool)]
class EstimatorConfig(Structure):
    """
    PMU estimator configuration structure.
    
    Contains all parameters needed to configure the PMU estimator algorithm.
    
    Attributes:
        n_cycles (int): Number of cycles of the fundamental frequency in the observation window
        f0 (int): Nominal fundamental frequency (Hz)
        frame_rate (int): Output frame rate (frames per second)
        fs (int): Sampling frequency (Hz)
        n_bins (int): Number of DFT bins to compute
        P (int): Number of iterations for the iterative algorithm
        Q (int): Number of interference tones for enhanced algorithm
        iter_eipdft (bool): Enable/disable iterative enhanced interpolated DFT
        interf_trig (float): Interference detection trigger threshold
        rocof_thresh (array): ROCOF thresholds for different conditions (3 values)
        rocof_low_pass_coeffs (array): Low-pass filter coefficients for ROCOF (3 values)
    """
    _fields_ = [("n_cycles", c_uint),
                ("f0", c_uint),
                ("frame_rate", c_uint),
                ("fs", c_uint),
                ("n_bins", c_uint),
                ("P", c_uint),
                ("Q", c_uint),
                ("iter_eipdft", c_bool),
                ("interf_trig", c_float),
                ("rocof_thresh", c_float * 3),
                ("rocof_low_pass_coeffs", c_float * 3)]

    def __init__(self, n_cycles, f0, frame_rate, fs, n_bins, P, Q, interf_trig, rocof_thresh, rocof_low_pass_coeffs, iter_eipdft = False):
        """
        Initialize estimator configuration.
        
        Args:
            n_cycles (int): Number of fundamental frequency cycles in window
            f0 (int): Nominal frequency (Hz)
            frame_rate (int): Frame rate (fps)
            fs (int): Sampling frequency (Hz)
            n_bins (int): Number of DFT bins
            P (int): Number of iterations
            Q (int): Number of interference tones
            interf_trig (float): Interference threshold
            rocof_thresh (list): ROCOF thresholds [3 values]
            rocof_low_pass_coeffs (list): Filter coefficients [3 values]
            iter_eipdft (bool): Enable iterative enhanced ipDFT (default: False)
        """
        self.n_cycles = n_cycles
        self.f0 = f0
        self.frame_rate = frame_rate
        self.fs = fs
        self.n_bins = n_bins
        self.P = P
        self.Q = Q
        self.iter_eipdft = iter_eipdft
        self.interf_trig = interf_trig
        for i in range(3):
            self.rocof_thresh[i] = rocof_thresh[i]
        for i in range(3):
            self.rocof_low_pass_coeffs[i] = rocof_low_pass_coeffs[i]

class PMUEstimator:
    """
    Python interface to the SynchroPMU C library.
    
    This class provides a high-level interface to the PMU estimator, handling
    library loading, initialization, and estimation operations.
    
    Multiple independent instances can be created and used simultaneously.
    
    Attributes:
        CONFIG_FROM_INI (bool): Configuration mode flag for INI file
        CONFIG_FROM_STRUCT (bool): Configuration mode flag for structure
    
    Example:
        >>> pmu = PMUEstimator()
        >>> pmu.configure_from_ini("config.ini")
        >>> result = pmu.estimate(signal_samples, 0.5)
        >>> print(f"Frequency: {result['freq']} Hz")
    """

    # Values of the configuration modes
    CONFIG_FROM_INI = True
    CONFIG_FROM_STRUCT = False

    def __init__(self, lib_path = None):
        """
        Initialize PMU estimator and load the C library.
        
        Args:
            lib_path (str, optional): Path to the shared library. If None, uses
                platform-specific default paths:
                - Linux: /usr/local/lib/libpmu_estimator.so
                - macOS: /usr/local/lib/libpmu_estimator.dylib
                - Windows: C:\\Program Files\\PmuEstimator\\lib\\libpmu_estimator.dll
        
        Raises:
            ValueError: If the platform is not supported
            FileNotFoundError: If the library is not found at the specified path
        """

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

        self.lib.pmu_init.argtypes = [POINTER(PmuContext), c_void_p, c_bool]
        self.lib.pmu_init.restype = c_int

        self.lib.pmu_estimate.argtypes = [POINTER(PmuContext), POINTER(c_float), c_float, POINTER(PmuFrame)]
        self.lib.pmu_estimate.restype = c_int

        self.lib.pmu_deinit.argtypes = [POINTER(PmuContext)]
        self.lib.pmu_deinit.restype = c_int

        self.ctx = PmuContext()
        self.ctx.pmu_initialized = 0

    def __del__(self):
        """Destructor - deinitializes the PMU estimator."""
        return self.lib.pmu_deinit(byref(self.ctx))
    
    def deinit(self):
        """
        Deinitialize the PMU estimator and free resources.
        
        Returns:
            int: 0 on success, -1 on error
        """
        return self.lib.pmu_deinit(byref(self.ctx))

    def configure_from_ini(self, ini_file_path):
        """
        Configure the PMU estimator from an INI file.
        
        Args:
            ini_file_path (str): Path to the configuration INI file
        
        Returns:
            int: 0 on success, -1 on error
        """
        ini_path_bytes = ini_file_path.encode('utf-8') if isinstance(ini_file_path, str) else ini_file_path
        return self.lib.pmu_init(byref(self.ctx), ini_path_bytes, self.CONFIG_FROM_INI)

    def configure_from_class(self, config):
        """
        Configure the PMU estimator from an EstimatorConfig object.
        
        Args:
            config (EstimatorConfig): Configuration object with all parameters
        
        Returns:
            int: 0 on success, -1 on error
        """
        return self.lib.pmu_init(byref(self.ctx), byref(config), self.CONFIG_FROM_STRUCT)

    def estimate(self, input_signal_window, mid_window_fracsec):
        """
        Estimate synchrophasor, frequency, and ROCOF from input signal.
        
        This is the main estimation function that processes a window of samples
        and returns the estimated PMU frame.
        
        Args:
            input_signal_window (list or array): Input signal samples. Length must
                match the configured window size (n_cycles * fs / f0).
            mid_window_fracsec (float): Fractional second timestamp (relative to PPS)
                of the midpoint of the observation window, used for phase correction.
        
        Returns:
            dict or None: Dictionary with keys 'amp', 'ph', 'freq', 'rocof' on success,
                None on error. Dictionary values:
                - 'amp' (float): Phasor amplitude
                - 'ph' (float): Phasor phase (radians)
                - 'freq' (float): Frequency (Hz)
                - 'rocof' (float): Rate of Change of Frequency (Hz/s)
        
        Example:
            >>> result = pmu.estimate(signal_window, 0.5)
            >>> if result:
            ...     print(f"Frequency: {result['freq']} Hz")
            ...     print(f"ROCOF: {result['rocof']} Hz/s")
        """
        frame = PmuFrame()
        input_signal = (c_float * len(input_signal_window))(*input_signal_window)
        result = self.lib.pmu_estimate(byref(self.ctx), input_signal, c_float(mid_window_fracsec), byref(frame))

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

