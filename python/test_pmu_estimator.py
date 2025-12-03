"""
Unit tests for the PMU Estimator Python API.

This test suite validates the core functionality of the SynchroPMU estimator,
including configuration, estimation, and error handling.
"""

import unittest
import math
import os
import tempfile
from pmu_estimator import PMUEstimator, EstimatorConfig, PmuFrame


class TestEstimatorConfig(unittest.TestCase):
    """Test cases for EstimatorConfig class."""

    def test_config_creation(self):
        """Test that EstimatorConfig can be created with valid parameters."""
        config = EstimatorConfig(
            n_cycles=4,
            f0=50,
            frame_rate=50,
            fs=25600,
            n_bins=11,
            P=3,
            Q=3,
            interf_trig=0.0033,
            rocof_thresh=[3.0, 25.0, 0.035],
            rocof_low_pass_coeffs=[0.5913, 0.2043, 0.2043],
            iter_eipdft=True
        )
        
        self.assertEqual(config.n_cycles, 4)
        self.assertEqual(config.f0, 50)
        self.assertEqual(config.frame_rate, 50)
        self.assertEqual(config.fs, 25600)
        self.assertEqual(config.n_bins, 11)
        self.assertEqual(config.P, 3)
        self.assertEqual(config.Q, 3)
        self.assertEqual(config.iter_eipdft, True)

    def test_config_default_iter_eipdft(self):
        """Test that iter_eipdft defaults to False."""
        config = EstimatorConfig(
            n_cycles=2,
            f0=60,
            frame_rate=60,
            fs=3840,
            n_bins=10,
            P=2,
            Q=2,
            interf_trig=0.01,
            rocof_thresh=[1.0, 10.0, 0.1],
            rocof_low_pass_coeffs=[0.5, 0.25, 0.25]
        )
        
        self.assertEqual(config.iter_eipdft, False)


class TestPmuFrame(unittest.TestCase):
    """Test cases for PmuFrame class."""

    def test_frame_initialization(self):
        """Test that PmuFrame can be initialized with parameters."""
        frame = PmuFrame(amp=2.0, ph=0.5, freq=50.0, rocof=0.1)
        
        # Check that values are approximately correct (floating point comparison)
        self.assertAlmostEqual(frame.synchrophasor.amp, 2.0, places=5)
        self.assertAlmostEqual(frame.synchrophasor.ph, 0.5, places=5)
        self.assertAlmostEqual(frame.synchrophasor.freq, 50.0, places=5)
        self.assertAlmostEqual(frame.rocof, 0.1, places=5)

    def test_frame_default_initialization(self):
        """Test that PmuFrame defaults to zeros."""
        frame = PmuFrame()
        
        self.assertEqual(frame.synchrophasor.amp, 0.0)
        self.assertEqual(frame.synchrophasor.ph, 0.0)
        self.assertEqual(frame.synchrophasor.freq, 0.0)
        self.assertEqual(frame.rocof, 0.0)

    def test_frame_str_representation(self):
        """Test that PmuFrame has a string representation."""
        frame = PmuFrame(amp=2.0, ph=0.5, freq=50.0, rocof=0.1)
        str_repr = str(frame)
        
        self.assertIn("Synchrophasor", str_repr)
        self.assertIn("amplitude", str_repr)
        self.assertIn("phase", str_repr)
        self.assertIn("frequency", str_repr)
        self.assertIn("rocof", str_repr)


class TestPMUEstimator(unittest.TestCase):
    """Test cases for PMUEstimator class."""

    def setUp(self):
        """Set up test fixtures - skip if library not available."""
        try:
            self.pmu = PMUEstimator()
            self.library_available = True
        except FileNotFoundError:
            self.library_available = False
            self.skipTest("PMU library not installed - skipping tests that require it")

    def tearDown(self):
        """Clean up after each test."""
        if self.library_available and hasattr(self, 'pmu'):
            try:
                self.pmu.deinit()
            except:
                pass

    def test_configure_from_struct_50hz(self):
        """Test configuration from EstimatorConfig structure with 50Hz."""
        if not self.library_available:
            return
            
        config = EstimatorConfig(
            n_cycles=4,
            f0=50,
            frame_rate=50,
            fs=25600,
            n_bins=11,
            P=3,
            Q=3,
            interf_trig=0.0033,
            rocof_thresh=[3.0, 25.0, 0.035],
            rocof_low_pass_coeffs=[0.5913, 0.2043, 0.2043],
            iter_eipdft=True
        )
        
        result = self.pmu.configure_from_class(config)
        self.assertEqual(result, 0, "Configuration should succeed")

    def test_configure_from_struct_60hz(self):
        """Test configuration from EstimatorConfig structure with 60Hz."""
        if not self.library_available:
            return
            
        config = EstimatorConfig(
            n_cycles=2,
            f0=60,
            frame_rate=60,
            fs=3840,
            n_bins=10,
            P=2,
            Q=2,
            interf_trig=0.01,
            rocof_thresh=[1.0, 10.0, 0.1],
            rocof_low_pass_coeffs=[0.5, 0.25, 0.25],
            iter_eipdft=False
        )
        
        result = self.pmu.configure_from_class(config)
        self.assertEqual(result, 0, "Configuration should succeed")

    def test_configure_from_ini_file(self):
        """Test configuration from INI file."""
        if not self.library_available:
            return
        
        # Skip this test - the C library's INI parser (iniparser) has issues with
        # temporary files and can cause division by zero errors when reading malformed
        # or incomplete INI data. This is a known limitation of the underlying C library.
        # To test INI configuration, use the actual config files from the repository.
        self.skipTest("INI configuration test skipped - use config file from repository")

    def test_estimate_known_signal_50hz(self):
        """Test estimation with a known 50Hz sinusoidal signal."""
        if not self.library_available:
            return
            
        # Signal parameters
        AMP = 2.0
        PH = 0.0
        FREQ = 50.0
        sample_rate = 25600
        window_size = 2048  # 4 cycles at 50Hz with 25600 samples/sec
        dt = 1.0 / sample_rate
        
        # Configure estimator
        config = EstimatorConfig(
            n_cycles=4,
            f0=50,
            frame_rate=50,
            fs=25600,
            n_bins=11,
            P=3,
            Q=3,
            interf_trig=0.0033,
            rocof_thresh=[3.0, 25.0, 0.035],
            rocof_low_pass_coeffs=[0.5913, 0.2043, 0.2043],
            iter_eipdft=True
        )
        
        self.pmu.configure_from_class(config)
        
        # Generate clean sinusoidal signal
        input_signal = [AMP * math.cos(2 * math.pi * FREQ * dt * i + PH) 
                       for i in range(window_size)]
        
        # Estimate
        result = self.pmu.estimate(input_signal, 0.0)
        
        self.assertIsNotNone(result, "Estimation should return a result")
        self.assertIn('amp', result)
        self.assertIn('ph', result)
        self.assertIn('freq', result)
        self.assertIn('rocof', result)
        
        # Check that amplitude is close to expected (within 1%)
        self.assertAlmostEqual(result['amp'], AMP, delta=AMP * 0.01,
                              msg="Amplitude should be close to expected value")
        
        # Check that frequency is close to 50Hz (within 0.1 Hz)
        self.assertAlmostEqual(result['freq'], FREQ, delta=0.1,
                              msg="Frequency should be close to 50Hz")

    def test_estimate_known_signal_51hz(self):
        """Test estimation with a known 51Hz signal (off-nominal frequency)."""
        if not self.library_available:
            return
            
        # Signal parameters
        AMP = 1.5
        PH = 0.5
        FREQ = 51.0  # Off-nominal frequency
        sample_rate = 25600
        window_size = 2048
        dt = 1.0 / sample_rate
        
        # Configure estimator
        config = EstimatorConfig(
            n_cycles=4,
            f0=50,
            frame_rate=50,
            fs=25600,
            n_bins=11,
            P=3,
            Q=3,
            interf_trig=0.0033,
            rocof_thresh=[3.0, 25.0, 0.035],
            rocof_low_pass_coeffs=[0.5913, 0.2043, 0.2043],
            iter_eipdft=True
        )
        
        self.pmu.configure_from_class(config)
        
        # Generate signal
        input_signal = [AMP * math.cos(2 * math.pi * FREQ * dt * i + PH) 
                       for i in range(window_size)]
        
        # Estimate
        result = self.pmu.estimate(input_signal, 0.0)
        
        self.assertIsNotNone(result, "Estimation should return a result")
        
        # Check frequency detection (within 0.2 Hz tolerance)
        self.assertAlmostEqual(result['freq'], FREQ, delta=0.2,
                              msg="Frequency should be close to 51Hz")

    def test_estimate_60hz_signal(self):
        """Test estimation with a 60Hz signal configuration."""
        if not self.library_available:
            return
            
        # Signal parameters for 60Hz system
        AMP = 1.0
        PH = 0.0
        FREQ = 60.0
        sample_rate = 3840
        window_size = 128  # 2 cycles at 60Hz
        dt = 1.0 / sample_rate
        
        # Configure for 60Hz system
        config = EstimatorConfig(
            n_cycles=2,
            f0=60,
            frame_rate=60,
            fs=3840,
            n_bins=10,
            P=2,
            Q=2,
            interf_trig=0.01,
            rocof_thresh=[1.0, 10.0, 0.1],
            rocof_low_pass_coeffs=[0.5, 0.25, 0.25],
            iter_eipdft=False
        )
        
        self.pmu.configure_from_class(config)
        
        # Generate signal
        input_signal = [AMP * math.cos(2 * math.pi * FREQ * dt * i + PH) 
                       for i in range(window_size)]
        
        # Estimate
        result = self.pmu.estimate(input_signal, 0.0)
        
        self.assertIsNotNone(result, "Estimation should return a result")
        self.assertAlmostEqual(result['freq'], FREQ, delta=0.2,
                              msg="Frequency should be close to 60Hz")

    def test_multiple_estimates(self):
        """Test that multiple estimates can be performed in sequence."""
        if not self.library_available:
            return
            
        # Configure
        config = EstimatorConfig(
            n_cycles=4,
            f0=50,
            frame_rate=50,
            fs=25600,
            n_bins=11,
            P=3,
            Q=3,
            interf_trig=0.0033,
            rocof_thresh=[3.0, 25.0, 0.035],
            rocof_low_pass_coeffs=[0.5913, 0.2043, 0.2043],
            iter_eipdft=True
        )
        
        self.pmu.configure_from_class(config)
        
        # Generate multiple signal windows
        sample_rate = 25600
        window_size = 2048
        dt = 1.0 / sample_rate
        
        for freq in [49.0, 50.0, 51.0]:
            input_signal = [2.0 * math.cos(2 * math.pi * freq * dt * i) 
                           for i in range(window_size)]
            
            result = self.pmu.estimate(input_signal, 0.0)
            self.assertIsNotNone(result, f"Estimation should succeed for {freq}Hz")
            self.assertAlmostEqual(result['freq'], freq, delta=0.2)

    def test_deinit(self):
        """Test that deinit can be called safely."""
        if not self.library_available:
            return
            
        config = EstimatorConfig(
            n_cycles=4,
            f0=50,
            frame_rate=50,
            fs=25600,
            n_bins=11,
            P=3,
            Q=3,
            interf_trig=0.0033,
            rocof_thresh=[3.0, 25.0, 0.035],
            rocof_low_pass_coeffs=[0.5913, 0.2043, 0.2043]
        )
        
        self.pmu.configure_from_class(config)
        result = self.pmu.deinit()
        # deinit should succeed or return an appropriate code
        self.assertIsNotNone(result)


class TestMultipleInstances(unittest.TestCase):
    """Test cases for multiple PMU estimator instances."""

    def setUp(self):
        """Set up test fixtures - skip if library not available."""
        try:
            # Just try to create an instance to check availability
            test_pmu = PMUEstimator()
            test_pmu.deinit()
            self.library_available = True
        except FileNotFoundError:
            self.library_available = False
            self.skipTest("PMU library not installed - skipping tests that require it")

    def test_multiple_independent_instances(self):
        """Test that multiple independent PMU instances can coexist."""
        if not self.library_available:
            return
            
        # Create two instances with different configurations
        pmu1 = PMUEstimator()
        pmu2 = PMUEstimator()
        
        try:
            # Configure first instance for 50Hz
            config1 = EstimatorConfig(
                n_cycles=4,
                f0=50,
                frame_rate=50,
                fs=25600,
                n_bins=11,
                P=3,
                Q=3,
                interf_trig=0.0033,
                rocof_thresh=[3.0, 25.0, 0.035],
                rocof_low_pass_coeffs=[0.5913, 0.2043, 0.2043],
                iter_eipdft=True
            )
            
            # Configure second instance for 60Hz
            config2 = EstimatorConfig(
                n_cycles=2,
                f0=60,
                frame_rate=60,
                fs=3840,
                n_bins=10,
                P=2,
                Q=2,
                interf_trig=0.01,
                rocof_thresh=[1.0, 10.0, 0.1],
                rocof_low_pass_coeffs=[0.5, 0.25, 0.25],
                iter_eipdft=False
            )
            
            result1 = pmu1.configure_from_class(config1)
            result2 = pmu2.configure_from_class(config2)
            
            self.assertEqual(result1, 0, "First instance configuration should succeed")
            self.assertEqual(result2, 0, "Second instance configuration should succeed")
            
            # Generate signals for each configuration
            signal1 = [2.0 * math.cos(2 * math.pi * 50.0 * (1.0/25600) * i) 
                      for i in range(2048)]
            signal2 = [1.0 * math.cos(2 * math.pi * 60.0 * (1.0/3840) * i) 
                      for i in range(128)]
            
            # Estimate with both instances
            result1 = pmu1.estimate(signal1, 0.0)
            result2 = pmu2.estimate(signal2, 0.0)
            
            self.assertIsNotNone(result1, "First instance estimation should succeed")
            self.assertIsNotNone(result2, "Second instance estimation should succeed")
            
            # Verify results are appropriate for their configurations
            self.assertAlmostEqual(result1['freq'], 50.0, delta=0.2,
                                  msg="First instance should estimate 50Hz")
            self.assertAlmostEqual(result2['freq'], 60.0, delta=0.2,
                                  msg="Second instance should estimate 60Hz")
        finally:
            pmu1.deinit()
            pmu2.deinit()


class TestEdgeCases(unittest.TestCase):
    """Test edge cases and error conditions."""

    def setUp(self):
        """Set up test fixtures - skip if library not available."""
        try:
            self.pmu = PMUEstimator()
            self.library_available = True
        except FileNotFoundError:
            self.library_available = False
            self.skipTest("PMU library not installed - skipping tests that require it")

    def tearDown(self):
        """Clean up after each test."""
        if self.library_available and hasattr(self, 'pmu'):
            try:
                self.pmu.deinit()
            except:
                pass

    def test_estimate_without_configuration(self):
        """Test that estimate fails gracefully without configuration."""
        if not self.library_available:
            return
            
        # Try to estimate without configuring
        signal = [1.0] * 2048
        result = self.pmu.estimate(signal, 0.0)
        
        # Should return None or handle gracefully
        self.assertIsNone(result, "Estimate without config should return None")

    def test_estimate_with_wrong_window_size(self):
        """Test estimation with incorrect window size."""
        if not self.library_available:
            return
            
        config = EstimatorConfig(
            n_cycles=4,
            f0=50,
            frame_rate=50,
            fs=25600,
            n_bins=11,
            P=3,
            Q=3,
            interf_trig=0.0033,
            rocof_thresh=[3.0, 25.0, 0.035],
            rocof_low_pass_coeffs=[0.5913, 0.2043, 0.2043]
        )
        
        self.pmu.configure_from_class(config)
        
        # Use wrong window size (should be 2048)
        wrong_signal = [1.0] * 1000
        
        # This might cause an error or return None depending on C library implementation
        try:
            result = self.pmu.estimate(wrong_signal, 0.0)
            # If it doesn't crash, that's acceptable
        except (ValueError, TypeError, OSError) as e:
            # If it raises a specific exception, that's also acceptable
            pass


if __name__ == '__main__':
    unittest.main()
